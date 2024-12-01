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


#include <aspect/particle/property/melt_particle.h>
#include <aspect/simulator.h>
#include <aspect/melt.h>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      template <int dim>
      void
      MeltParticle<dim>::initialize_one_particle_property(const Point<dim> &/*position*/,
                                                          std::vector<double> &data) const
      {
        data.push_back(0.0);
      }

      template <int dim>
      void
      MeltParticle<dim>::update_particle_properties(const ParticleUpdateInputs<dim> &inputs,
                                                    typename ParticleHandler<dim>::particle_iterator_range &particles) const
      {
        const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");

        unsigned int p = 0;
        for (auto &particle: particles)
          {
            if (inputs.solution[p][this->introspection().component_indices.compositional_fields[porosity_idx]] > threshold_for_melt_presence)
              particle.get_properties()[this->data_position] = 1.0;
            else
              particle.get_properties()[this->data_position] = 0.0;
            ++p;
          }
      }

      template <int dim>
      UpdateTimeFlags
      MeltParticle<dim>::need_update() const
      {
        return update_time_step;
      }

      template <int dim>
      UpdateFlags
      MeltParticle<dim>::get_update_flags (const unsigned int component) const
      {
        if (this->introspection().component_masks.compositions[component] == true)
          return update_values;

        return update_default;
      }

      template <int dim>
      std::vector<std::pair<std::string, unsigned int>>
      MeltParticle<dim>::get_property_information() const
      {
        std::vector<std::pair<std::string,unsigned int>> property_information (1,std::make_pair("melt_presence",1));
        return property_information;
      }

      template <int dim>
      void
      MeltParticle<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Melt particle");
        {
          prm.declare_entry ("Threshold for melt presence", "1e-3",
                             Patterns::Double (0., 1.),
                             "The minimum porosity that has to be present at the position of a particle "
                             "for it to be considered a melt particle (in the sense that the melt presence "
                             "property is set to 1).");
        }
        prm.leave_subsection();
      }


      template <int dim>
      void
      MeltParticle<dim>::parse_parameters (ParameterHandler &prm)
      {
        AssertThrow(this->introspection().compositional_name_exists("porosity"),
                    ExcMessage("Particle property melt particle only works if"
                               "there is a compositional field called porosity."));

        prm.enter_subsection("Melt particle");
        {
          threshold_for_melt_presence = prm.get_double ("Threshold for melt presence");
        }
        prm.leave_subsection();
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
      ASPECT_REGISTER_PARTICLE_PROPERTY(MeltParticle,
                                        "melt particle",
                                        "Implementation of a plugin in which the particle "
                                        "property is defined as presence of melt above a "
                                        "threshold, which can be set as an input parameter. "
                                        "This property is set to 0 if melt is not present and "
                                        "set to 1 if melt is present.")
    }
  }
}
