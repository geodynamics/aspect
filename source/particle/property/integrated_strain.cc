/*
  Copyright (C) 2015 - 2017 by the authors of the ASPECT code.

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

#include <aspect/particle/property/integrated_strain.h>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      template <int dim>
      void
      IntegratedStrain<dim>::initialize_one_particle_property(const Point<dim> &,
                                                              std::vector<double> &data) const
      {
        const static Tensor<2,dim> identity = unit_symmetric_tensor<dim>();
        for (unsigned int i = 0; i < Tensor<2,dim>::n_independent_components ; ++i)
          data.push_back(identity[Tensor<2,dim>::unrolled_to_component_indices(i)]);
      }

      template <int dim>
      void
      IntegratedStrain<dim>::update_one_particle_property(const unsigned int data_position,
                                                          const Point<dim> &,
                                                          const Vector<double> &,
                                                          const std::vector<Tensor<1,dim> > &gradients,
                                                          const ArrayView<double> &data) const
      {
        Tensor<2,dim> old_strain;
        for (unsigned int i = 0; i < Tensor<2,dim>::n_independent_components ; ++i)
          old_strain[Tensor<2,dim>::unrolled_to_component_indices(i)] = data[data_position + i];

        Tensor<2,dim> grad_u;
        for (unsigned int d=0; d<dim; ++d)
          grad_u[d] = gradients[d];

        const double dt = this->get_timestep();

        Tensor<2,dim> new_strain;

        // here we integrate the equation
        // new_deformation_gradient = velocity_gradient * old_deformation_gradient
        // using a RK4 integration scheme.
        const Tensor<2,dim> k1 = grad_u * old_strain * dt;
        new_strain = old_strain + 0.5*k1;

        const Tensor<2,dim> k2 = grad_u * new_strain * dt;
        new_strain = old_strain + 0.5*k2;

        const Tensor<2,dim> k3 = grad_u * new_strain * dt;
        new_strain = old_strain + k3;

        const Tensor<2,dim> k4 = grad_u * new_strain * dt;

        // the new strain is the rotated old strain plus the
        // strain of the current time step
        new_strain = old_strain + (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;

        for (unsigned int i = 0; i < Tensor<2,dim>::n_independent_components ; ++i)
          data[data_position + i] = new_strain[Tensor<2,dim>::unrolled_to_component_indices(i)];
      }

      template <int dim>
      UpdateTimeFlags
      IntegratedStrain<dim>::need_update() const
      {
        return update_time_step;
      }

      template <int dim>
      UpdateFlags
      IntegratedStrain<dim>::get_needed_update_flags () const
      {
        return update_gradients;
      }

      template <int dim>
      std::vector<std::pair<std::string, unsigned int> >
      IntegratedStrain<dim>::get_property_information() const
      {
        const unsigned int n_components = Tensor<2,dim>::n_independent_components;
        const std::vector<std::pair<std::string,unsigned int> > property_information (1,std::make_pair("integrated strain",n_components));
        return property_information;
      }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      ASPECT_REGISTER_PARTICLE_PROPERTY(IntegratedStrain,
                                        "integrated strain",
                                        "A plugin in which the particle property tensor is "
                                        "defined as the deformation gradient tensor "
                                        "$\\mathbf F$ this particle has experienced. "
                                        "$\\mathbf F$ can be polar-decomposed into the left stretching tensor "
                                        "$\\mathbf L$ (the finite strain we are interested in), and the "
                                        "rotation tensor $\\mathbf Q$. See the corresponding cookbook in "
                                        "the manual for more detailed information.")
    }
  }
}

