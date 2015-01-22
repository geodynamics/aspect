/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include <aspect/simulator.h>
#include <aspect/global.h>

#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/constraint_matrix.h>

#ifdef ASPECT_USE_PETSC
#include <deal.II/lac/solver_cg.h>
#else
#include <deal.II/lac/trilinos_solver.h>
#endif

#include <deal.II/lac/pointer_matrix.h>
#include <deal.II/base/tensor_function.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

namespace aspect
{
  namespace internal
  {
    using namespace dealii;



    /**
     * A class we use when setting up the data structures for nullspace removal
     * of the rotations in spherical or annular shells.
     */
    template<int dim>
    class Rotation : public TensorFunction<1,dim>
    {
      private:
        Tensor<1,dim> axis;

      public:
        Rotation(const unsigned int a)
          :
          axis(Tensor<1,dim>(Point<dim>::unit_vector(a)))
        {}

        virtual Tensor<1,dim> value (const Point<dim> &p) const
        {
          Tensor<1,dim> vel;
          if ( dim == 2)
            cross_product(vel, p);
          else
            cross_product(vel, axis, p);
          return vel;
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
        const unsigned int direction;

      public:
        Translation(const unsigned int d)
          :
          direction(d)
        {}

        virtual Tensor<1,dim> value(const Point<dim> &) const
        {
          return Point<dim>::unit_vector(direction);
        }
    };

  }

  template <int dim>
  void Simulator<dim>::setup_nullspace_removal()
  {
    if (parameters.nullspace_removal & NullspaceRemoval::linear_momentum)
      AssertThrow(false, ExcNotImplemented());

    std::vector<std_cxx1x::shared_ptr<TensorFunction<1,dim> > > funcs;

    if (parameters.nullspace_removal & NullspaceRemoval::net_rotation)
      {
        if (dim==2)
          funcs.push_back(std_cxx1x::shared_ptr<TensorFunction<1,dim> >(new internal::Rotation<dim>(0)));
        if (dim==3)
          for (unsigned int a=0; a<dim; ++a)
            funcs.push_back(std_cxx1x::shared_ptr<TensorFunction<1,dim> >(new internal::Rotation<dim>(a)));
      }

    if (parameters.nullspace_removal & NullspaceRemoval::net_translation_x)
      funcs.push_back(std_cxx1x::shared_ptr<TensorFunction<1,dim> >(new internal::Translation<dim>(0)));

    if (parameters.nullspace_removal & NullspaceRemoval::net_translation_y)
      funcs.push_back(std_cxx1x::shared_ptr<TensorFunction<1,dim> >(new internal::Translation<dim>(1)));

    if (parameters.nullspace_removal & NullspaceRemoval::net_translation_z)
      {
        //Only do z direction if dim == 3
        AssertThrow( dim == 3, ExcMessage("Can't remove z translational mode in 2 dimensions"));
        funcs.push_back(std_cxx1x::shared_ptr<TensorFunction<1,dim> >(new internal::Translation<dim>(2)));
      }

    if (funcs.size()>0)
      {
        Assert(introspection.block_indices.velocities != introspection.block_indices.pressure,
               ExcNotImplemented());

        net_rotations_translations.resize(funcs.size());
        for (unsigned int i=0; i<funcs.size(); ++i)
          {
            // for each of the null space dimensions, set up
            // a vector, fill it with an element of the null space,
            // and normalize it
            net_rotations_translations[i].reinit(
              introspection.index_sets.system_partitioning[introspection.block_indices.velocities],
              mpi_communicator);
            interpolate_onto_velocity_system(*funcs[i],
                                             net_rotations_translations[i]);
            net_rotations_translations[i] /= net_rotations_translations[i].l2_norm();
          }

      }
  }



  template <int dim>
  void Simulator<dim>::remove_nullspace(LinearAlgebra::BlockVector &relevant_dst,
                                        LinearAlgebra::BlockVector &tmp_distributed_stokes)
  {
    if (parameters.nullspace_removal & NullspaceRemoval::net_rotation ||
        parameters.nullspace_removal & NullspaceRemoval::net_translation )
      {
        Assert(introspection.block_indices.velocities != introspection.block_indices.pressure,
               ExcNotImplemented());

        for (unsigned int i=0; i<net_rotations_translations.size(); ++i)
          {
            // compute the magnitude of the solution vector in direction
            // of this null space vector and subtract the corresponding multiple
            const double power = net_rotations_translations[i]
                                 * tmp_distributed_stokes.block(introspection.block_indices.velocities);
            tmp_distributed_stokes.block(introspection.block_indices.velocities).sadd(1.0,
                                                                                      -1.0*power,
                                                                                      net_rotations_translations[i]);
          }
        relevant_dst.block(0) = tmp_distributed_stokes.block(0);
      }
    if (parameters.nullspace_removal & NullspaceRemoval::angular_momentum)
      {
        remove_net_angular_momentum( relevant_dst, tmp_distributed_stokes);
      }
  }



  template <>
  void
  Simulator<3>::remove_net_angular_momentum( LinearAlgebra::BlockVector &relevant_dst,
                                             LinearAlgebra::BlockVector &tmp_distributed_stokes )
  {
    AssertThrow(false, ExcNotImplemented());
  }


  template <>
  void
  Simulator<2>::remove_net_angular_momentum( LinearAlgebra::BlockVector &relevant_dst,
                                             LinearAlgebra::BlockVector &tmp_distributed_stokes )
  {
    Assert(introspection.block_indices.velocities != introspection.block_indices.pressure,
           ExcNotImplemented());

    // compute and remove angular momentum from velocity field, by computing
    // int rho V \cdot r_orth = omega  * int rho x^2
    const unsigned int dim=2;

    QGauss<dim> quadrature(parameters.stokes_velocity_degree+1);
    FEValues<dim> fe(mapping, finite_element, quadrature,
                     UpdateFlags(update_quadrature_points | update_JxW_values | update_values));

    DoFHandler<dim>::active_cell_iterator cell;
    std::vector<Point<dim> > q_points(quadrature.size());
    std::vector<Vector<double> > fe_vals(quadrature.size(), Vector<double>(finite_element.n_components()));

    // analogues to the moment of inertia and angular momentum for 2D
    double local_scalar_moment = 0;
    double local_scalar_angular_momentum = 0;

    //loop over all local cells
    for (cell = dof_handler.begin_active(); cell != dof_handler.end(); ++cell)
      if (cell->is_locally_owned())
        {
          fe.reinit (cell);
          q_points = fe.get_quadrature_points();
          fe.get_function_values(relevant_dst, fe_vals);

          // get the density at each quadrature point
          MaterialModel::Interface<dim>::MaterialModelInputs in(q_points.size(),
                                                                parameters.n_compositional_fields);
          MaterialModel::Interface<dim>::MaterialModelOutputs out(q_points.size(),
                                                                  parameters.n_compositional_fields);
          for (unsigned int i=0; i< q_points.size(); i++)
            {
              in.pressure[i] = fe_vals[i][introspection.component_indices.pressure];
              in.temperature[i] = fe_vals[i][introspection.component_indices.temperature];
              for (unsigned int c=0; c<parameters.n_compositional_fields; ++c)
                in.composition[i][c] = fe_vals[i][introspection.component_indices.compositional_fields[c]];
              in.position[i] = q_points[i];

            }
          material_model->evaluate(in, out);

          Point<dim> r_vec;
          Tensor<1,dim> velocity;

          // actually compute the moment of inertia and angular momentum
          for (unsigned int k=0; k< quadrature.size(); ++k)
            {
              // get the position and velocity at this quadrature point
              r_vec = q_points[k];
              for (unsigned int i=0; i<dim; ++i) velocity[i] = fe_vals[k][i];
              // Get the velocity perpendicular to the position vector
              Tensor<1,dim> r_perp;
              cross_product(r_perp, r_vec);
              Tensor<1,dim> v_perp = velocity - (velocity*r_vec)*r_vec/(r_vec.norm_square());

              // calculate a signed scalar angular momentum
              local_scalar_angular_momentum += fe.JxW(k) * velocity*r_perp * out.densities[k];
              // calculate a scalar moment of inertia
              local_scalar_moment += fe.JxW(k) * r_vec.norm_square() * out.densities[k];
            }
        }

    const double scalar_moment = Utilities::MPI::sum( local_scalar_moment, mpi_communicator);
    const double scalar_angular_momentum = Utilities::MPI::sum( local_scalar_angular_momentum, mpi_communicator);

    const double rotation_rate = scalar_angular_momentum/scalar_moment;  //solve for the rotation rate to cancel the angular momentum

    // Now construct a rotation vector with the desired rate and subtract it from our vector
    LinearAlgebra::Vector correction(tmp_distributed_stokes.block(introspection.block_indices.velocities));
    internal::Rotation<dim> rot(0);
    interpolate_onto_velocity_system(rot, correction);
    tmp_distributed_stokes.block(introspection.block_indices.velocities).sadd(1.0, -1.0*rotation_rate,correction);

    relevant_dst.block(0) = tmp_distributed_stokes.block(0);
  }

}





// explicit instantiation of the functions we implement in this file
namespace aspect
{
#define INSTANTIATE(dim) \
  template void Simulator<dim>::setup_nullspace_removal (); \
  template void Simulator<dim>::remove_nullspace (LinearAlgebra::BlockVector &,LinearAlgebra::BlockVector &vector);

  ASPECT_INSTANTIATE(INSTANTIATE)
}
