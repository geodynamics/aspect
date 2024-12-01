/*
  Copyright (C) 2015 - 2024 by the authors of the ASPECT code.

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

#include <aspect/particle/property/strain_rate.h>
#include <aspect/particle/property/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      template <int dim>
      void
      StrainRate<dim>::initialize_one_particle_property(const Point<dim> &,
                                                        std::vector<double> &data) const
      {
        const static Tensor<2,dim> identity = unit_symmetric_tensor<dim>();
        for (unsigned int i = 0; i < Tensor<2,dim>::n_independent_components ; ++i)
          data.push_back(identity[Tensor<2,dim>::unrolled_to_component_indices(i)]);

      }

      template <int dim>
      void
      StrainRate<dim>::update_particle_properties(const ParticleUpdateInputs<dim> &inputs,
                                                  typename ParticleHandler<dim>::particle_iterator_range &particles) const
      {
        unsigned int p = 0;
        for (auto &particle: particles)
          {
            const auto data = particle.get_properties();
            // Velocity gradients
            Tensor<2,dim> grad_u;
            for (unsigned int d=0; d<dim; ++d)
              grad_u[d] = inputs.gradients[p][d];

            // Calculate strain rate from velocity gradients
            const SymmetricTensor<2,dim> strain_rate = symmetrize (grad_u);

            for (unsigned int i = 0; i < Tensor<2,dim>::n_independent_components ; ++i)
              data[this->data_position + i] = strain_rate[Tensor<2,dim>::unrolled_to_component_indices(i)];
            ++p;
          }

      }

      template <int dim>
      UpdateTimeFlags
      StrainRate<dim>::need_update() const
      {
        return update_time_step;
      }

      template <int dim>
      UpdateFlags
      StrainRate<dim>::get_update_flags (const unsigned int component) const
      {
        if (this->introspection().component_masks.velocities[component] == true)
          return update_gradients;

        return update_default;
      }

      template <int dim>
      std::vector<std::pair<std::string, unsigned int>>
      StrainRate<dim>::get_property_information() const
      {
        const unsigned int n_components = Tensor<2,dim>::n_independent_components;
        const std::vector<std::pair<std::string,unsigned int>> property_information (1,std::make_pair("strainrate",n_components));
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
      ASPECT_REGISTER_PARTICLE_PROPERTY(StrainRate,
                                        "strain rate",
                                        "Implementation of a plugin in which the time evolution of "
                                        "strain rate is saved and stored on the particles.")
    }
  }
}
