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

#include <aspect/particle/property/pT_path.h>
#include <aspect/simulator_signals.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/initial_temperature/interface.h>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      template <int dim>
      void
      PTPath<dim>::initialize()
      {
        // Make sure we keep track of the initial temperature manager and
        // that it continues to live beyond the time when the simulator
        // class releases its pointer to it.
        initial_temperature = this->get_initial_temperature_manager_pointer();
      }




      template <int dim>
      void
      PTPath<dim>::initialize_one_particle_property(const Point<dim> &position,
                                                    std::vector<double> &data) const
      {
        // The following is strictly only correct whenever a particle is
        // created in the first time step. After that, taking pressure and
        // temperature from adiabatic and initial temperature objects is
        // not quite correct -- instead, we should be initializing from
        // the current pressure and temperature, but that is substantially
        // more complicated since we are not passed this information.
        //
        // The issue is probably not terribly important because at least for
        // all following time steps, we set temperature and pressure to
        // their correct then-current values.
        data.push_back(this->get_adiabatic_conditions().pressure(position));
        data.push_back(initial_temperature->initial_temperature(position));
      }



      template <int dim>
      void
      PTPath<dim>::update_particle_properties(const ParticleUpdateInputs<dim> &inputs,
                                              typename ParticleHandler<dim>::particle_iterator_range &particles) const
      {
        unsigned int p = 0;
        for (auto &particle: particles)
          {
            particle.get_properties()[this->data_position]   = inputs.solution[p][this->introspection().component_indices.pressure];
            particle.get_properties()[this->data_position+1] = inputs.solution[p][this->introspection().component_indices.temperature];
            ++p;
          }
      }



      template <int dim>
      UpdateTimeFlags
      PTPath<dim>::need_update() const
      {
        return update_output_step;
      }



      template <int dim>
      UpdateFlags
      PTPath<dim>::get_update_flags (const unsigned int component) const
      {
        if (this->introspection().component_masks.pressure[component] == true ||
            this->introspection().component_masks.temperature[component] == true)
          return update_values;

        return update_default;
      }



      template <int dim>
      std::vector<std::pair<std::string, unsigned int>>
      PTPath<dim>::get_property_information() const
      {
        return {{"p", 1}, {"T", 1}};
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
      ASPECT_REGISTER_PARTICLE_PROPERTY(PTPath,
                                        "pT path",
                                        "Implementation of a plugin in which the particle "
                                        "property is defined as the current pressure and "
                                        "temperature at this position. This can be used "
                                        "to generate pressure-temperature paths of "
                                        "material points over time.")
    }
  }
}
