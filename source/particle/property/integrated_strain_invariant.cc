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

#include <aspect/particle/property/integrated_strain_invariant.h>
#include <aspect/utilities.h>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      template <int dim>
      void
      IntegratedStrainInvariant<dim>::initialize_one_particle_property(const Point<dim> &,
                                                                       std::vector<double> &data) const
      {
        data.push_back(0.0);
      }



      template <int dim>
      void
      IntegratedStrainInvariant<dim>::update_particle_properties(const ParticleUpdateInputs<dim> &inputs,
                                                                 typename ParticleHandler<dim>::particle_iterator_range &particles) const
      {

        unsigned int p = 0;
        for (auto &particle: particles)
          {
            // Integrated strain invariant from prior time step
            const auto data = particle.get_properties();
            const double old_strain = data[this->data_position];

            // Current timestep
            const double dt = this->get_timestep();

            // Velocity gradients
            Tensor<2,dim> grad_u;
            for (unsigned int d=0; d<dim; ++d)
              grad_u[d] = inputs.gradients[p][this->introspection().component_indices.velocities[d]];

            // Calculate strain rate from velocity gradients
            const SymmetricTensor<2,dim> strain_rate = symmetrize (grad_u);

            // Calculate strain rate second invariant
            const double edot_ii = std::sqrt(std::max(-Utilities::Tensors::consistent_second_invariant_of_deviatoric_tensor(Utilities::Tensors::consistent_deviator(strain_rate)), 0.));

            // New strain is the old strain plus dt*edot_ii
            const double new_strain = old_strain + dt*edot_ii;
            data[this->data_position] = new_strain;
            ++p;
          }
      }



      template <int dim>
      UpdateTimeFlags
      IntegratedStrainInvariant<dim>::need_update() const
      {
        return update_time_step;
      }



      template <int dim>
      UpdateFlags
      IntegratedStrainInvariant<dim>::get_update_flags (const unsigned int component) const
      {
        if (this->introspection().component_masks.velocities[component] == true)
          return update_gradients;

        return update_default;
      }



      template <int dim>
      std::vector<std::pair<std::string, unsigned int>>
      IntegratedStrainInvariant<dim>::get_property_information() const
      {
        const std::vector<std::pair<std::string,unsigned int>> property_information (1,std::make_pair("integrated strain invariant",1));
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
      ASPECT_REGISTER_PARTICLE_PROPERTY(IntegratedStrainInvariant,
                                        "integrated strain invariant",
                                        "A plugin in which the particle property is defined as "
                                        "the finite strain invariant ($\\varepsilon_{ii}$). "
                                        "This property is calculated with the timestep ($dt$) and "
                                        "the second invariant of the deviatoric strain rate tensor "
                                        "($\\dot{\\varepsilon}_{ii}$), where the value at time step $n$ is "
                                        "$\\varepsilon_{ii}^{n} = \\varepsilon_{ii}^{n-1} + "
                                        "dt\\dot{\\varepsilon}_{ii}$.")
    }
  }
}
