/*
  Copyright (C) 2011 - 2021 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include <aspect/simulator.h>
#include <aspect/global.h>

#include <deal.II/lac/solver_gmres.h>

#include <deal.II/lac/trilinos_solver.h>

#include <deal.II/base/tensor_function.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

namespace aspect
{
  namespace internal
  {

    /**
     * A class we use when setting up the data structures for nullspace removal
     * of the rotations in spherical or annular shells.
     */
    template <int dim>
    class Rotation : public TensorFunction<1,dim>
    {
      private:
        const Tensor<1,dim> axis;

      public:
        // Constructor for TensorFunction that takes cartesian direction (1,2, or 3)
        // and creates a solid body rotation around that axis.
        Rotation(const unsigned int a)
          :
          axis(Tensor<1,dim>(Point<dim>::unit_vector(a)))
        {}

        // Constructor for TensorFunction that takes an axis
        // and creates a solid body rotation around that axis.
        Rotation(const Tensor<1,dim> &rotation_axis)
          :
          axis(rotation_axis)
        {}

        Tensor<1,dim> value (const Point<dim> &p) const override
        {
          if ( dim == 2)
            return cross_product_2d(p);
          else
            return cross_product_3d(axis, p);
        }
    };


    /**
     * A class we use when setting up the data structures for nullspace removal
     * of the translations in box-like geometries.
     */
    template <int dim>
    class Translation : public TensorFunction<1,dim>
    {
      private:
        const Tensor<1,dim> translation;

      public:
        // Constructor for TensorFunction that takes a Cartesian direction (1,2, or 3)
        // and creates a translation along that axis
        Translation(const unsigned int d)
          :
          translation( Point<dim>::unit_vector(d) )
        {}

        // Constructor for TensorFunction that takes a vector
        // and creates a translation along that vector
        Translation(const Tensor<1,dim> &t)
          :
          translation(t)
        {}

        Tensor<1,dim> value(const Point<dim> &) const override
        {
          return translation;
        }
    };
  }



  template <int dim>
  void Simulator<dim>::setup_nullspace_constraints(AffineConstraints<double> &constraints)
  {
    if (!(parameters.nullspace_removal & (NullspaceRemoval::linear_momentum
                                          | NullspaceRemoval::net_translation)))
      return;

    // Note: We want to add a single Dirichlet zero constraint for each
    // translation direction. This is complicated by the fact that we need to
    // find a DoF that is not already constrained. In parallel the constraint
    // needs to be added on all processors where it is locally_relevant and
    // all processors need to agree on the index.

    // First find candidates for DoF indices to constrain for each velocity component.
    types::global_dof_index vel_idx[dim];
    {
      for (unsigned int d=0; d<dim; ++d)
        vel_idx[d] = numbers::invalid_dof_index;

      unsigned int n_left_to_find = dim;

      std::vector<types::global_dof_index> local_dof_indices (finite_element.dofs_per_cell);
      typename DoFHandler<dim>::active_cell_iterator cell;
      for (const auto &cell : dof_handler.active_cell_iterators())
        if (cell->is_locally_owned())
          {
            cell->get_dof_indices (local_dof_indices);

            for (unsigned int i=0; i<finite_element.dofs_per_cell; ++i)
              {
                const unsigned int component = finite_element.system_to_component_index(i).first;

                if (component < introspection.component_indices.velocities[0]
                    || component > introspection.component_indices.velocities[dim-1])
                  continue; // only look at velocity

                const unsigned int velocity_component = component - introspection.component_indices.velocities[0];

                if (vel_idx[velocity_component] != numbers::invalid_dof_index)
                  continue; // already found one

                const types::global_dof_index idx = local_dof_indices[i];

                if (constraints.can_store_line(idx) && !constraints.is_constrained(idx))
                  {
                    vel_idx[velocity_component] = idx;
                    --n_left_to_find;
                  }

                // are we done searching?
                if (n_left_to_find == 0)
                  break; // exit inner loop
              }

            if (n_left_to_find == 0)
              break; // exit outer loop
          }

    }


    const unsigned int flags[] = {(NullspaceRemoval::linear_momentum_x
                                   |NullspaceRemoval::net_translation_x),
                                  (NullspaceRemoval::linear_momentum_y
                                   |NullspaceRemoval::net_translation_y),
                                  (NullspaceRemoval::linear_momentum_z
                                   |NullspaceRemoval::net_translation_z)
                                 };

    for (unsigned int d=0; d<dim; ++d)
      if (parameters.nullspace_removal & flags[d])
        {
          // Make a reduction to find the smallest index (processors that
          // found a larger candidate just happened to not be able to store
          // that index with the minimum value). Note that it is possible that
          // some processors might not be able to find a potential DoF, for
          // example because they don't own any DoFs. On those processors we
          // will use dof_handler.n_dofs() when building the minimum (larger
          // than any valid DoF index).
          const types::global_dof_index global_idx = dealii::Utilities::MPI::min(
                                                       (vel_idx[d] != numbers::invalid_dof_index)
                                                       ?
                                                       vel_idx[d]
                                                       :
                                                       dof_handler.n_dofs(),
                                                       mpi_communicator);

          Assert(global_idx < dof_handler.n_dofs(),
                 ExcMessage("Error, couldn't find a velocity DoF to constrain."));

          // Finally set this DoF to zero (if the current MPI process
          // cares about it):
          if (constraints.can_store_line(global_idx))
            {
              Assert(!constraints.is_constrained((global_idx)),
                     ExcInternalError());
              constraints.add_line(global_idx);
            }
        }
  }


  template <int dim>
  void Simulator<dim>::remove_nullspace(LinearAlgebra::BlockVector &relevant_dst,
                                        LinearAlgebra::BlockVector &tmp_distributed_stokes)
  {
    if (parameters.nullspace_removal & NullspaceRemoval::angular_momentum)
      {
        // use_constant_density = false, remove net angular momentum
        remove_net_angular_momentum( false, relevant_dst, tmp_distributed_stokes);
      }
    if (parameters.nullspace_removal & NullspaceRemoval::linear_momentum)
      {
        // use_constant_density = false, remove net momentum
        remove_net_linear_momentum( false, relevant_dst, tmp_distributed_stokes);
      }
    if (parameters.nullspace_removal & NullspaceRemoval::net_rotation)
      {
        // use_constant_density = true, remove net rotation
        remove_net_angular_momentum( true, relevant_dst, tmp_distributed_stokes);
      }
    if (parameters.nullspace_removal & NullspaceRemoval::net_surface_rotation)
      {
        // use_constant_density = true, remove net rotation
        // limit_to_top_faces = true, to limit to surface motions
        remove_net_angular_momentum( true, relevant_dst, tmp_distributed_stokes, true);
      }
    if (parameters.nullspace_removal & NullspaceRemoval::net_translation)
      {
        // use_constant_density = true, remove net translation
        remove_net_linear_momentum( true, relevant_dst, tmp_distributed_stokes);
      }
  }

  template <int dim>
  void Simulator<dim>::remove_net_linear_momentum( const bool use_constant_density,
                                                   LinearAlgebra::BlockVector &relevant_dst,
                                                   LinearAlgebra::BlockVector &tmp_distributed_stokes )
  {
    Assert(introspection.block_indices.velocities != introspection.block_indices.pressure,
           ExcNotImplemented());

    // compute and remove net linear momentum from velocity field, by computing
    // \int \rho (v + v_const) = 0

    QGauss<dim> quadrature(parameters.stokes_velocity_degree+1);
    const unsigned int n_q_points = quadrature.size();
    FEValues<dim> fe(*mapping, finite_element, quadrature,
                     UpdateFlags(update_quadrature_points | update_JxW_values | update_values | update_gradients));

    Tensor<1,dim> local_momentum;
    double local_mass = 0.0;


    // Vectors for evaluating the finite element solution
    std::vector<std::vector<double>> composition_values (introspection.n_compositional_fields,
                                                         std::vector<double> (n_q_points));
    std::vector< Tensor<1,dim>> velocities( n_q_points );

    typename DoFHandler<dim>::active_cell_iterator cell;
    // loop over all local cells
    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        {
          fe.reinit (cell);

          // get the velocity at each quadrature point
          fe[introspection.extractors.velocities].get_function_values (relevant_dst, velocities);

          // get the density at each quadrature point if necessary
          MaterialModel::MaterialModelInputs<dim> in(n_q_points,
                                                     introspection.n_compositional_fields);
          MaterialModel::MaterialModelOutputs<dim> out(n_q_points,
                                                       introspection.n_compositional_fields);
          if ( ! use_constant_density)
            {
              fe[introspection.extractors.pressure].get_function_values (relevant_dst, in.pressure);
              fe[introspection.extractors.temperature].get_function_values (relevant_dst, in.temperature);
              in.velocity = velocities;
              fe[introspection.extractors.pressure].get_function_gradients (relevant_dst, in.pressure_gradient);
              for (unsigned int c=0; c<introspection.n_compositional_fields; ++c)
                fe[introspection.extractors.compositional_fields[c]].get_function_values(relevant_dst,
                                                                                         composition_values[c]);

              for (unsigned int i=0; i<n_q_points; ++i)
                {
                  in.position[i] = fe.quadrature_point(i);
                  for (unsigned int c=0; c<introspection.n_compositional_fields; ++c)
                    in.composition[i][c] = composition_values[c][i];
                }
              material_model->evaluate(in, out);
            }

          // actually compute the momentum and mass
          for (unsigned int k=0; k<n_q_points; ++k)
            {
              // get the density at this quadrature point
              const double rho = (use_constant_density ? 1.0 : out.densities[k]);

              local_momentum += velocities[k] * rho * fe.JxW(k);
              local_mass += rho * fe.JxW(k);
            }
        }

    // Calculate the total mass and velocity correction
    const double mass = Utilities::MPI::sum( local_mass, mpi_communicator);
    Tensor<1,dim> velocity_correction = Utilities::MPI::sum(local_momentum, mpi_communicator)/mass;

    // We may only want to remove the nullspace for a single component, so zero out
    // the velocity correction if it is not selected by the NullspaceRemoval flag
    if (use_constant_density) // disable translation correction
      {
        if ( !(parameters.nullspace_removal & NullspaceRemoval::net_translation_x) )
          velocity_correction[0] = 0.0;  // don't correct x translation
        if ( !(parameters.nullspace_removal & NullspaceRemoval::net_translation_y) )
          velocity_correction[1] = 0.0;  // don't correct y translation
        if ( !(parameters.nullspace_removal & NullspaceRemoval::net_translation_z) && dim == 3 )
          velocity_correction[2] = 0.0;  // don't correct z translation
      }
    else // disable momentum correction
      {
        if ( !(parameters.nullspace_removal & NullspaceRemoval::linear_momentum_x) )
          velocity_correction[0] = 0.0;  // don't correct x translation
        if ( !(parameters.nullspace_removal & NullspaceRemoval::linear_momentum_y) )
          velocity_correction[1] = 0.0;  // don't correct y translation
        if ( !(parameters.nullspace_removal & NullspaceRemoval::linear_momentum_z) && dim == 3 )
          velocity_correction[2] = 0.0;  // don't correct z translation
      }

    // vector for storing the correction to the velocity field
    LinearAlgebra::Vector correction(tmp_distributed_stokes.block(introspection.block_indices.velocities));

    // Now construct a translation vector with the desired rate and subtract it from our vector
    internal::Translation<dim> translation( velocity_correction );
    interpolate_onto_velocity_system(translation, correction);
    tmp_distributed_stokes.block(introspection.block_indices.velocities).add(-1.0,correction);

    // copy into the locally relevant vector
    relevant_dst.block(introspection.block_indices.velocities) =
      tmp_distributed_stokes.block(introspection.block_indices.velocities);
  }



  template <int dim>
  RotationProperties<dim>
  Simulator<dim>::compute_net_angular_momentum(const bool use_constant_density,
                                               const LinearAlgebra::BlockVector &solution,
                                               const bool limit_to_top_faces) const
  {
    // compute the momentum from velocity field, by computing
    // \int \rho u \cdot r_orth = \omega  * \int \rho x^2    ( 2 dimensions)
    // \int \rho r \times u =  I \cdot \omega  (3 dimensions)

    QGauss<dim> quadrature(parameters.stokes_velocity_degree+1);
    QGauss<dim-1> surface_quadrature(parameters.stokes_velocity_degree+1);

    const unsigned int n_q_points = (limit_to_top_faces == false) ? quadrature.size() : surface_quadrature.size();
    UpdateFlags flags = update_quadrature_points | update_JxW_values | update_values | update_gradients;

    FEValues<dim> fe_values (*mapping, finite_element, quadrature, flags);
    FEFaceValues<dim> fe_face_values (*mapping, finite_element, surface_quadrature, flags);
    FEValuesBase<dim> &fe((limit_to_top_faces == false)
                          ?
                          dynamic_cast<FEValuesBase<dim> &>(fe_values)
                          :
                          dynamic_cast<FEValuesBase<dim> &>(fe_face_values));

    // moment of inertia and angular momentum for 3D
    SymmetricTensor<2,dim> local_moment_of_inertia;
    Tensor<1,dim> local_angular_momentum;

    // analogues to the moment of inertia and angular momentum for 2D
    double local_scalar_moment_of_inertia = 0.0;
    double local_scalar_angular_momentum = 0.0;

    // Structures for evaluating the velocities and the material model
    MaterialModel::MaterialModelInputs<dim> in(n_q_points,
                                               introspection.n_compositional_fields);
    MaterialModel::MaterialModelOutputs<dim> out(n_q_points,
                                                 introspection.n_compositional_fields);

    // loop over all local cells
    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        {
          if (limit_to_top_faces == false)
            fe_values.reinit (cell);
          else
            {
              // We only want the output at the top boundary, so only compute it if the current cell
              // has a face at the top boundary.
              bool cell_at_top_boundary = false;
              for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
                if (cell->at_boundary(f) &&
                    (geometry_model->translate_id_to_symbol_name (cell->face(f)->boundary_id()) == "top"))
                  {
                    Assert(cell_at_top_boundary == false,
                           ExcInternalError("Error, more than one top surface found in a cell."));

                    cell_at_top_boundary = true;
                    fe_face_values.reinit(cell,f);
                  }

              if (cell_at_top_boundary == false)
                continue;
            }

          const std::vector<Point<dim>> &q_points = fe.get_quadrature_points();

          if (use_constant_density == false)
            {
              // Set use_strain_rates to false since we don't need viscosity
              in.reinit(fe, cell, introspection, solution, false);
              material_model->evaluate(in, out);
            }
          else
            {
              // Get the velocity at each quadrature point
              fe[introspection.extractors.velocities].get_function_values (solution, in.velocity);
            }

          // actually compute the moment of inertia and angular momentum
          for (unsigned int k=0; k<n_q_points; ++k)
            {
              // get the position and density at this quadrature point
              const Point<dim> r_vec = q_points[k];
              const double rho = (use_constant_density ? 1.0 : out.densities[k]);
              const double JxW = fe.JxW(k);

              if (dim == 2)
                {
                  // Get the velocity perpendicular to the position vector
                  const Tensor<1,dim> r_perp = cross_product_2d(r_vec);

                  local_scalar_angular_momentum += in.velocity[k] * r_perp * rho * JxW;
                  local_scalar_moment_of_inertia += r_vec.norm_square() * rho * JxW;
                }
              else
                {
                  const Tensor<1,dim> r_cross_v = cross_product_3d(r_vec, in.velocity[k]);

                  local_angular_momentum += r_cross_v * rho * JxW;
                  local_moment_of_inertia += (r_vec.norm_square() * unit_symmetric_tensor<dim>() - symmetrize(outer_product(r_vec,r_vec))) * rho * JxW;
                }
            }
        }

    RotationProperties<dim> properties;

    // Sum up the local contributions and solve for the overall rotation
    if (dim == 2)
      {
        properties.scalar_moment_of_inertia = Utilities::MPI::sum(local_scalar_moment_of_inertia, mpi_communicator);
        properties.scalar_angular_momentum = Utilities::MPI::sum(local_scalar_angular_momentum, mpi_communicator);
        properties.scalar_rotation = properties.scalar_angular_momentum / properties.scalar_moment_of_inertia;
      }
    else
      {
        properties.tensor_moment_of_inertia = Utilities::MPI::sum(local_moment_of_inertia,
                                                                  mpi_communicator);
        properties.tensor_angular_momentum = Utilities::MPI::sum(local_angular_momentum, mpi_communicator );
        const SymmetricTensor<2,dim> inverse_moment (invert(Tensor<2,dim>(properties.tensor_moment_of_inertia)));
        properties.tensor_rotation = - inverse_moment * properties.tensor_angular_momentum;
      }

    return properties;
  }



  template <int dim>
  void Simulator<dim>::remove_net_angular_momentum( const bool use_constant_density,
                                                    LinearAlgebra::BlockVector &relevant_dst,
                                                    LinearAlgebra::BlockVector &tmp_distributed_stokes,
                                                    const bool limit_to_top_faces)
  {
    Assert(introspection.block_indices.velocities != introspection.block_indices.pressure,
           ExcNotImplemented());

    RotationProperties<dim> rotation_properties = compute_net_angular_momentum(use_constant_density,
                                                                               relevant_dst,
                                                                               limit_to_top_faces);

    // vector for storing the correction to the velocity field
    LinearAlgebra::Vector correction(tmp_distributed_stokes.block(introspection.block_indices.velocities));

    if (dim == 2)
      {
        // Now construct a rotation vector with the desired rate and subtract it from our vector
        const internal::Rotation<dim> rot(0);
        interpolate_onto_velocity_system(rot, correction);
        tmp_distributed_stokes.block(introspection.block_indices.velocities).add(-1.0*rotation_properties.scalar_rotation,correction);
      }
    else
      {
        // Remove that rotation from the solution vector
        const internal::Rotation<dim> rot(rotation_properties.tensor_rotation);
        interpolate_onto_velocity_system(rot, correction);
        tmp_distributed_stokes.block(introspection.block_indices.velocities).add(1.0,correction);
      }

    // copy into the locally relevant vector
    relevant_dst.block(introspection.block_indices.velocities) =
      tmp_distributed_stokes.block(introspection.block_indices.velocities);
  }

}





// explicit instantiation of the functions we implement in this file
namespace aspect
{
#define INSTANTIATE(dim) \
  template struct RotationProperties<dim>; \
  template void Simulator<dim>::remove_nullspace (LinearAlgebra::BlockVector &,LinearAlgebra::BlockVector &vector); \
  template void Simulator<dim>::setup_nullspace_constraints (AffineConstraints<double> &);

  ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
}
