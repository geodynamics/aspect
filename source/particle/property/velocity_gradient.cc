/*
  Copyright (C) 2015 - 2025 by the authors of the ASPECT code.

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

#include <aspect/particle/property/velocity_gradient.h>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      
      template <int dim>
      void
      VelocityGradient<dim>::initialize_one_particle_property(const Point<dim> &,
                                                              std::vector<double> &data) const
      {
        const static Tensor<2,dim> identity = unit_symmetric_tensor<dim>();
        for (unsigned int i = 0; i < Tensor<2,dim>::n_independent_components ; ++i)
          data.push_back(identity[Tensor<2,dim>::unrolled_to_component_indices(i)]);
      }



      template <int dim>
      void
      VelocityGradient<dim>::update_particle_property(const unsigned int data_position,
                                                      const Vector<double> &/*solution*/,
                                                      const std::vector<Tensor<1,dim>> &gradients,
                                                      typename ParticleHandler<dim>::particle_iterator &particle) const
      {
        ArrayView<double> data = particle->get_properties();
        // Velocity gradients
        Tensor<2,dim> grad_u;
        for (unsigned int d=0; d<dim; ++d)
          grad_u[d] = gradients[d];
        for (unsigned int i = 0; i < Tensor<2,dim>::n_independent_components ; ++i)
          data[data_position + i] = grad_u[Tensor<2,dim>::unrolled_to_component_indices(i)];
      }



      template <int dim>
      UpdateTimeFlags
      VelocityGradient<dim>::need_update() const
      {
        return update_time_step;
      }



      template <int dim>
      UpdateFlags
      VelocityGradient<dim>::get_needed_update_flags () const
      {
        return update_gradients;
      }


      
      template <int dim>
      std::vector<std::pair<std::string, unsigned int>>
      VelocityGradient<dim>::get_property_information() const
      {
        const unsigned int n_components = Tensor<2,dim>::n_independent_components;
        const std::vector<std::pair<std::string,unsigned int>> property_information (1,std::make_pair("velocity gradient",n_components));
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
      ASPECT_REGISTER_PARTICLE_PROPERTY(VelocityGradient,
                                        "velocity gradient",
                                        "Implementation of a plugin in which the particle "
                                        "property is defined as the recent velocity gradient "
                                        "at this position.")
    }
  }
}
