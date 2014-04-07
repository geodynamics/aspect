/*
  Copyright (C) 2011, 2012, 2013 by the authors of the ASPECT code.

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
/*  $Id: solver.cc 2325 2014-02-27 23:21:30Z bangerth $  */


#include <aspect/simulator.h>
#include <aspect/global.h>

#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/constraint_matrix.h>

#ifdef USE_PETSC
#include <deal.II/lac/solver_cg.h>
#else
#include <deal.II/lac/trilinos_solver.h>
#endif

#include <deal.II/lac/pointer_matrix.h>
#include <deal.II/base/tensor_function.h>

namespace aspect
{
  namespace internal
  {
    using namespace dealii;



    template<int dim>
    class Rotation : public TensorFunction<1,dim>
      {
        private:
          Tensor<1,dim> axis;
        public:
          Rotation(const unsigned int a) : axis(Tensor<1,dim>(Point<dim>::unit_vector(a))) {}
          virtual Tensor<1,dim> value (const Point<dim> &p) const { Tensor<1,dim> vel; if( dim == 2) cross_product(vel, p); else cross_product(vel, axis, p); return vel;}
      };

    template <int dim>
    class Translation : public TensorFunction<1,dim>
    {
      private:
        const unsigned int direction;
      public:
        Translation(const unsigned int d) : direction(d) {}
        virtual Tensor<1,dim> value(const Point<dim> &) const { return Tensor<1,dim>(Point<dim>::unit_vector(direction)); }
    };

  }

  template <int dim>
  void Simulator<dim>::setup_nullspace_removal()
  {
    if (parameters.nullspace_removal & NullspaceRemoval::angular_momentum)
        AssertThrow(false, ExcNotImplemented());
    if (parameters.nullspace_removal & NullspaceRemoval::translational_momentum)
        AssertThrow(false, ExcNotImplemented());

    std::vector<std_cxx1x::shared_ptr<TensorFunction<1,dim> > > funcs;

    if (parameters.nullspace_removal & NullspaceRemoval::net_rotation)
      {
        if (dim==2)
          funcs.push_back(std_cxx1x::shared_ptr<TensorFunction<1,dim> >(new internal::Rotation<dim>(0)));
        if (dim==3)
          for(unsigned int a=0; a<dim; ++a)
            funcs.push_back(std_cxx1x::shared_ptr<TensorFunction<1,dim> >(new internal::Rotation<dim>(a)));
      }

    if (parameters.nullspace_removal & NullspaceRemoval::net_translation)
      {
          for(unsigned int a=0; a<dim; ++a)
            funcs.push_back(std_cxx1x::shared_ptr<TensorFunction<1,dim> >(new internal::Translation<dim>(a)));
      }

    if (funcs.size()>0)
      {
        net_rotations_translations.resize(funcs.size());
        for (unsigned int i=0;i<funcs.size();++i)
          net_rotations_translations[i].reinit(
              introspection.index_sets.system_partitioning[introspection.block_indices.velocities],
              mpi_communicator);

        unsigned int idx = 0;
        for (unsigned int i=0;i<funcs.size();++i)
          {
            interpolate_onto_velocity_system(*funcs[i],
                net_rotations_translations[i]);
            net_rotations_translations[i] /= net_rotations_translations[i].l2_norm();
          }

      }
  }

  template <int dim>
  void Simulator<dim>::remove_nullspace(LinearAlgebra::BlockVector &vector)
  {
    if (parameters.nullspace_removal & NullspaceRemoval::net_rotation ||
        parameters.nullspace_removal & NullspaceRemoval::net_translation)
      {
        for(unsigned int i=0; i<net_rotations_translations.size(); ++i)
        {
           double power = net_rotations_translations[i]
                          * vector.block(introspection.block_indices.velocities);
           vector.block(introspection.block_indices.velocities).sadd(1.0,
                     -1.0*power,
                     net_rotations_translations[i]);
           pcout << "removing NULLSPACE " << i << " power: " << power << std::endl;
        }
      }
  }

}





// explicit instantiation of the functions we implement in this file
namespace aspect
{
#define INSTANTIATE(dim) \
  template void Simulator<dim>::setup_nullspace_removal (); \
  template void Simulator<dim>::remove_nullspace (LinearAlgebra::BlockVector &vector);

  ASPECT_INSTANTIATE(INSTANTIATE)
}
