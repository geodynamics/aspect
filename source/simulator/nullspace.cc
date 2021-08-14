/*
  Copyright (C) 2011 - 2022 by the authors of the ASPECT code.

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


#include <aspect/global.h>
#include <aspect/simulator.h>
#include <aspect/geometry_model/sphere.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/mesh_deformation/interface.h>

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
     * Add constraints to @p constraints for the active components
     * specified by @p mask for the DoFs in the support point given by
     * @p location assumed on the boundary of the domain. This will set their
     * value to zero.
     *
     * If successful (a support point at the position @p location was found),
     * this function returns @p true.
     *
     * The DoFs are found by iterating over all faces on the boundary
     * and computing the distance to @p location. We consider all locally relevant
     * DoFs here, so that the resulting AffineConstraints object is consistent between
     * processors.
     */
    template <int dim, int spacedim=dim>
    bool
    constrain_point_on_boundary_to_zero_active(AffineConstraints<double> &constraints,
                                               const DoFHandler<dim,spacedim> &dof_handler,
                                               const Mapping<dim> &mapping,
                                               const Point<dim> &location,
                                               const ComponentMask &mask)
    {
      const auto &fe = dof_handler.get_fe();
      const std::vector<Point<dim - 1>> &unit_support_points = fe.get_unit_face_support_points();
      const Quadrature<dim - 1> quadrature(unit_support_points);
      const unsigned int dofs_per_face = fe.dofs_per_face;
      std::vector<types::global_dof_index> face_dofs(dofs_per_face);

      FEFaceValues<dim, spacedim> fe_face_values(mapping,
                                                 fe,
                                                 quadrature,
                                                 update_quadrature_points);

      for (const auto &cell : dof_handler.active_cell_iterators())
        if (cell->at_boundary()
            &&
            (cell->is_locally_owned() || cell->is_ghost()))
          for (unsigned int face_no = 0;
               face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
            if (cell->at_boundary(face_no))
              {
                const typename DoFHandler<dim, spacedim>::face_iterator face = cell->face(face_no);
                face->get_dof_indices(face_dofs);
                fe_face_values.reinit(cell, face_no);

                bool found = false;
                for (unsigned int i = 0; i < face_dofs.size(); ++i)
                  {
                    const unsigned int component = fe.face_system_to_component_index(i).first;
                    if (mask[component])
                      {
                        const Point<dim> position = fe_face_values.quadrature_point(i);
                        if (position.distance(location) < 1e-6*cell->diameter())
                          {
                            found = true;
                            if (!constraints.is_constrained(face_dofs[i]) &&
                                constraints.can_store_line(face_dofs[i]))
                              constraints.add_line(face_dofs[i]);
                          }
                      }
                  }

                // Success! No reason to look at any other faces:
                if (found)
                  return true;
              }

      return false;
    }

    /**
     * Like above, but on the given multigrid level @p level.
     */
    template <int dim, int spacedim=dim>
    bool
    constrain_point_on_boundary_to_zero_on_level(AffineConstraints<double> &constraints,
                                                 const DoFHandler<dim,spacedim> &dof_handler,
                                                 const Mapping<dim> &mapping,
                                                 const Point<dim> &location,
                                                 const ComponentMask &mask,
                                                 const int level)
    {
      const auto &fe = dof_handler.get_fe();
      const std::vector<Point<dim - 1>> &unit_support_points = fe.get_unit_face_support_points();
      const Quadrature<dim - 1> quadrature(unit_support_points);
      const unsigned int dofs_per_face = fe.dofs_per_face;
      std::vector<types::global_dof_index> face_dofs(dofs_per_face);

      FEFaceValues<dim, spacedim> fe_face_values(mapping,
                                                 fe,
                                                 quadrature,
                                                 update_quadrature_points);

      std::set<types::boundary_id>::iterator b_id;

      for (const auto &cell : dof_handler.cell_iterators_on_level(level))
        if (
#if DEAL_II_VERSION_GTE(9,4,0)
          cell->is_locally_owned_on_level() ||  cell->is_ghost_on_level()
#else
          cell->level_subdomain_id() != numbers::artificial_subdomain_id
#endif
        )
          for (unsigned int face_no = 0;
               face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
            if (cell->at_boundary(face_no))
              {
                const typename DoFHandler<dim, spacedim>::level_face_iterator face = cell->face(face_no);
                face->get_mg_dof_indices(level, face_dofs);
                fe_face_values.reinit(cell, face_no);

                bool found = false;
                for (unsigned int i = 0; i < face_dofs.size(); ++i)
                  {
                    const unsigned int component = fe.face_system_to_component_index(i).first;
                    if (mask[component])
                      {
                        const Point<dim> position = fe_face_values.quadrature_point(i);
                        if (position.distance(location) < 1e-6*cell->diameter())
                          {
                            found = true;
                            if (!constraints.is_constrained(face_dofs[i]) &&
                                constraints.can_store_line(face_dofs[i]))
                              constraints.add_line(face_dofs[i]);
                          }
                      }
                  }
                // Success! No reason to look at any other faces:
                if (found)
                  return true;

              }

      return false;
    }



    /**
     * Calls either the active or the multilevel version of constrain_point_on_boundary_to_zero() above.
     */
    template <int dim, int spacedim=dim>
    bool
    constrain_point_on_boundary_to_zero(
      AffineConstraints<double> &constraints,
      const DoFHandler<dim,spacedim> &dof_handler,
      const Mapping<dim> &mapping,
      const Point<dim> &location,
      const ComponentMask &mask,
      const bool multilevel_constraints,
      const unsigned int level = numbers::invalid_unsigned_int)
    {
      if (multilevel_constraints)
        {
          Assert(level!=numbers::invalid_unsigned_int, ExcInternalError());
          return constrain_point_on_boundary_to_zero_on_level(constraints, dof_handler, mapping, location, mask, level);
        }
      else
        {
          Assert(level==numbers::invalid_unsigned_int, ExcInternalError());
          return constrain_point_on_boundary_to_zero_active(constraints, dof_handler, mapping, location, mask);
        }
    }



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

  /**
   * Add a number of constraints to @p constraints to remove the
   * rotational nullspace of the problem either for the active DoFs
   * or on the specified multigrid level @p level.
   *
   * This code currently only works for spherical shells and spheres
   * in 2d and 3d. A spherical object has 1 rotational mode in 2d and
   * 3 in 3d (rotations around x, y, and z axis). We will fix an
   * appropriate number of velocity components to get a unique
   * solution. We will need to identify 1 point in 2d and 2 points in
   * 3d to do that. Mesh deformation complicates the situation, as we
   * can not easily identify points on that surface.
   */
  template <int dim>
  void setup_rotational_constraints(AffineConstraints<double> &constraints,
                                    const GeometryModel::Interface<dim> &geometry_model,
                                    const MeshDeformation::MeshDeformationHandler<dim> *mesh_deformation,
                                    const DoFHandler<dim> &dof_handler,
                                    const Mapping<dim> &mapping,
                                    const MPI_Comm &mpi_communicator,
                                    const bool multilevel_constraints,
                                    const unsigned int level = numbers::invalid_unsigned_int)
  {
    // We assert() that we find the DoFs we are trying to fix, but
    // this only works for active DoFs (level==-1) and not level
    // DoFs.  When we have adaptive refinement, where the finer
    // levels do not contain cells touching the chosen point, we
    // won't be able to find a DoF on these finer levels. This
    // works okay, though.

    const bool is_spherical_shell = Plugins::plugin_type_matches<const GeometryModel::SphericalShell<dim>>(geometry_model);
    const bool is_sphere = Plugins::plugin_type_matches<const GeometryModel::Sphere<dim>>(geometry_model);
    Assert(is_spherical_shell || is_sphere,
           ExcNotImplemented("Nullspace constraints for rotation are currently only implemented for the "
                             "Sphere and SphericalShell geometries."));

    const bool mesh_deformation_enabled = (mesh_deformation != nullptr);
    bool mesh_deformation_at_top = false;
    bool mesh_deformation_at_bottom = false;
    if (mesh_deformation_enabled)
      {
        // Check if mesh deformation happens at the top boundary, bottom boundary, of both:
        const std::set<types::boundary_id> ids = mesh_deformation->get_active_mesh_deformation_boundary_indicators();

        if (ids.find(geometry_model.get_symbolic_boundary_names_map()["top"]) != ids.end())
          mesh_deformation_at_top = true;
        if (is_spherical_shell && ids.find(geometry_model.get_symbolic_boundary_names_map()["bottom"]) != ids.end())
          mesh_deformation_at_bottom = true;
      }

    if (mesh_deformation_at_top && is_sphere)
      {
        // Not sure what to do in this situation. For now, let's not add a constraint and hope that it works.
        return;
      }
    if (mesh_deformation_at_top && mesh_deformation_at_bottom && is_spherical_shell)
      {
        // Again, not much we can do here...
        return;
      }

    // Use the bottom surface unless we have mesh deformation at the bottom:
    const double point_depth = ((is_sphere || (is_spherical_shell && mesh_deformation_at_bottom)) ? 0.0 : geometry_model.maximal_depth());

    if (dim==2)
      {
        // Pick a point at the desired depth. We assume that the
        // representative point is in positive y direction, so we can
        // fix the x component of velocity. This is true for shell and
        // spherical shell implementations.
        const Point<dim> location = geometry_model.representative_point(point_depth);
        Assert(location[0] == 0. && location[1]>0., ExcInternalError());

        ComponentMask x_velocity_mask(dof_handler.get_fe().n_components(), false);
        x_velocity_mask.set(0, true);

        const bool success = internal::constrain_point_on_boundary_to_zero(constraints,
                                                                           dof_handler,
                                                                           mapping,
                                                                           location,
                                                                           x_velocity_mask,
                                                                           multilevel_constraints,
                                                                           level);

        const bool global_success = (Utilities::MPI::max(success?1:0, mpi_communicator) == 1);
        AssertThrow(multilevel_constraints || global_success,
                    ExcInternalError("Could not find the specified support point for adding a nullspace constraint."));
      }
    else if (dim==3)
      {
        // Pick a point at the desired depth. We assume that the
        // representative point is on the positive z axis. This is
        // true for shell and spherical shell implementations. This
        // way we can fix the tangential components (x and y
        // velocity).
        const Point<dim> location1 = geometry_model.representative_point(point_depth);
        Assert(location1[0] == 0. && location1[1] == 0. && location1[dim-1]>0., ExcInternalError());

        {
          ComponentMask x_and_y_velocity_mask(dof_handler.get_fe().n_components(), false);
          x_and_y_velocity_mask.set(0, true);
          x_and_y_velocity_mask.set(1, true);
          bool success = internal::constrain_point_on_boundary_to_zero(constraints,
                                                                       dof_handler,
                                                                       mapping,
                                                                       location1,
                                                                       x_and_y_velocity_mask,
                                                                       multilevel_constraints,
                                                                       level);
          const bool global_success = (Utilities::MPI::max(success?1:0, mpi_communicator) == 1);
          AssertThrow(multilevel_constraints || global_success,
                      ExcInternalError("Could not find the specified support point for adding a nullspace constraint."));
        }

        {
          // Construct point 2 by rotating location1 so it is on the positive x axis. Then constrain y velocity.
          Point<dim> location2;
          location2[0] = location1[dim-1];
          ComponentMask y_velocity_mask(dof_handler.get_fe().n_components(), false);
          y_velocity_mask.set(1, true);
          bool success = internal::constrain_point_on_boundary_to_zero(constraints,
                                                                       dof_handler,
                                                                       mapping,
                                                                       location2,
                                                                       y_velocity_mask,
                                                                       multilevel_constraints,
                                                                       level);
          const bool global_success = (Utilities::MPI::max(success?1:0, mpi_communicator) == 1);
          AssertThrow(multilevel_constraints || global_success,
                      ExcInternalError("Could not find the specified support point for adding a nullspace constraint."));
        }
      }


  }

  template <int dim>
  void Simulator<dim>::setup_nullspace_constraints(AffineConstraints<double> &constraints)
  {

    if ((parameters.nullspace_removal & Parameters<dim>::NullspaceRemoval::any_rotation)
        &&
        parameters.constrain_rotational_nullspace)
      {
        setup_rotational_constraints(constraints,
                                     *geometry_model,
                                     mesh_deformation.get(),
                                     dof_handler,
                                     *mapping,
                                     mpi_communicator,
                                     /* multilevel_constraints = */ false);
      }

    if (parameters.nullspace_removal & NullspaceRemoval::any_translation)
      {
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
  }


  template <int dim>
  void Simulator<dim>::setup_nullspace_constraints(AffineConstraints<double> &constraints,
                                                   const DoFHandler<dim> &dof_handler,
                                                   const int level)
  {
    if ((parameters.nullspace_removal & Parameters<dim>::NullspaceRemoval::any_rotation)
        &&
        parameters.constrain_rotational_nullspace)
      {
        const Mapping<dim> &current_mapping =
          (mesh_deformation) ? mesh_deformation->get_level_mapping(level) : *mapping;

        setup_rotational_constraints(constraints,
                                     *geometry_model,
                                     mesh_deformation.get(),
                                     dof_handler,
                                     current_mapping,
                                     mpi_communicator,
                                     /* multilevel_constraints = */ true,
                                     level);
      }
  }



  template <int dim>
  void Simulator<dim>::remove_nullspace(LinearAlgebra::BlockVector &relevant_dst,
                                        LinearAlgebra::BlockVector &tmp_distributed_stokes)
  {
    if (parameters.nullspace_removal & NullspaceRemoval::angular_momentum)
      {
        remove_net_angular_momentum( /* use_constant_density = */ false, // remove net momentum
                                                                  relevant_dst,
                                                                  tmp_distributed_stokes);
      }
    if (parameters.nullspace_removal & NullspaceRemoval::net_rotation)
      {
        remove_net_angular_momentum( /* use_constant_density = */ true, // remove net rotation
                                                                  relevant_dst,
                                                                  tmp_distributed_stokes);
      }
    if (parameters.nullspace_removal & NullspaceRemoval::net_surface_rotation)
      {
        remove_net_angular_momentum( /* use_constant_density = */ true, // remove net rotation
                                                                  relevant_dst,
                                                                  tmp_distributed_stokes,
                                                                  /* limit_to_top_faces = */ true);
      }

    if (parameters.nullspace_removal & NullspaceRemoval::linear_momentum)
      {
        remove_net_linear_momentum( /* use_constant_density = */ false, // remove net momentum
                                                                 relevant_dst,
                                                                 tmp_distributed_stokes);
      }
    if (parameters.nullspace_removal & NullspaceRemoval::net_translation)
      {
        remove_net_linear_momentum( /* use_constant_density = */ true, // remove net translation
                                                                 relevant_dst,
                                                                 tmp_distributed_stokes);
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
    const Quadrature<dim> &quadrature = introspection.quadratures.velocities;
    const unsigned int n_q_points = quadrature.size();
    FEValues<dim> fe(*mapping, finite_element, quadrature,
                     UpdateFlags(update_quadrature_points | update_JxW_values | update_values | update_gradients));

    Tensor<1,dim> local_momentum;
    double local_mass = 0.0;


    // Vectors for evaluating the finite element solution
    std::vector<std::vector<double>> composition_values (introspection.n_compositional_fields,
                                                          std::vector<double> (n_q_points));
    std::vector<Tensor<1,dim>> velocities( n_q_points );

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
          in.requested_properties = MaterialModel::MaterialProperties::density;

          if (!use_constant_density)
            {
              fe[introspection.extractors.pressure].get_function_values(relevant_dst, in.pressure);
              fe[introspection.extractors.temperature].get_function_values(relevant_dst, in.temperature);
              in.velocity = velocities;
              fe[introspection.extractors.pressure].get_function_gradients(relevant_dst, in.pressure_gradient);
              for (unsigned int c = 0; c < introspection.n_compositional_fields; ++c)
                fe[introspection.extractors.compositional_fields[c]].get_function_values(relevant_dst,
                                                                                         composition_values[c]);

              for (unsigned int i = 0; i < n_q_points; ++i)
                {
                  in.position[i] = fe.quadrature_point(i);
                  for (unsigned int c = 0; c < introspection.n_compositional_fields; ++c)
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

    const Quadrature<dim> &quadrature = introspection.quadratures.velocities;
    const Quadrature<dim-1> &surface_quadrature = introspection.face_quadratures.velocities;

    const unsigned int n_q_points = (limit_to_top_faces == false) ? quadrature.size() : surface_quadrature.size();
    UpdateFlags flags = update_quadrature_points | update_JxW_values | update_values | update_gradients;

    FEValues<dim> fe_values (*mapping, finite_element, quadrature, flags);
    FEFaceValues<dim> fe_face_values (*mapping, finite_element, surface_quadrature, flags);
    FEValuesBase<dim> &fe((limit_to_top_faces == false)
                          ?
                          dynamic_cast<FEValuesBase<dim> &>(fe_values)
                          :
                          dynamic_cast<FEValuesBase<dim> &>(fe_face_values));

    // moment of inertia and angular momentum for 3d
    SymmetricTensor<2,dim> local_moment_of_inertia;
    Tensor<1,dim> local_angular_momentum;

    // analogues to the moment of inertia and angular momentum for 2d
    double local_scalar_moment_of_inertia = 0.0;
    double local_scalar_angular_momentum = 0.0;

    // Structures for evaluating the velocities and the material model
    MaterialModel::MaterialModelInputs<dim> in(n_q_points,
                                               introspection.n_compositional_fields);
    MaterialModel::MaterialModelOutputs<dim> out(n_q_points,
                                                 introspection.n_compositional_fields);
    in.requested_properties = MaterialModel::MaterialProperties::density;

    // loop over all local cells
    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        {
          if (limit_to_top_faces == false)
            fe_values.reinit(cell);
          else
            {
              // We only want the output at the top boundary, so only compute it if the current cell
              // has a face at the top boundary.
              bool cell_at_top_boundary = false;
              for (const unsigned int f : cell->face_indices())
                if (cell->at_boundary(f) &&
                    (geometry_model->translate_id_to_symbol_name(cell->face(f)->boundary_id()) == "top"))
                  {
                    Assert(cell_at_top_boundary == false,
                           ExcInternalError("Error, more than one top surface found in a cell."));

                    cell_at_top_boundary = true;
                    fe_face_values.reinit(cell, f);
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
              fe[introspection.extractors.velocities].get_function_values(solution, in.velocity);
            }

          // actually compute the moment of inertia and angular momentum
          for (unsigned int k = 0; k < n_q_points; ++k)
            {
              // get the position and density at this quadrature point
              const Point<dim> r_vec = q_points[k];
              const double rho = (use_constant_density ? 1.0 : out.densities[k]);
              const double JxW = fe.JxW(k);

              if (dim == 2)
                {
                  // Get the velocity perpendicular to the position vector
                  const Tensor<1, dim> r_perp = cross_product_2d(r_vec);

                  local_scalar_angular_momentum += in.velocity[k] * r_perp * rho * JxW;
                  local_scalar_moment_of_inertia += r_vec.norm_square() * rho * JxW;
                }
              else
                {
                  const Tensor<1, dim> r_cross_v = cross_product_3d(r_vec, in.velocity[k]);

                  local_angular_momentum += r_cross_v * rho * JxW;
                  local_moment_of_inertia += (r_vec.norm_square() * unit_symmetric_tensor<dim>() - symmetrize(outer_product(r_vec, r_vec))) * rho * JxW;
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
  template void Simulator<dim>::remove_nullspace (LinearAlgebra::BlockVector &,LinearAlgebra::BlockVector &); \
  template void Simulator<dim>::setup_nullspace_constraints (AffineConstraints<double> &); \
  template void Simulator<dim>::setup_nullspace_constraints (AffineConstraints<double> &, const DoFHandler<dim>&, const int);

  ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
}
