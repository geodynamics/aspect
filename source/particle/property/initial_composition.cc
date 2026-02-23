/*
  Copyright (C) 2015 - 2023 by the authors of the ASPECT code.

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

#include <aspect/particle/property/initial_composition.h>
#include <aspect/initial_composition/interface.h>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      template <int dim>
      void
      InitialComposition<dim>::initialize_one_particle_property(const Point<dim> &position,
                                                                std::vector<double> &data) const
      {
        const std::string prefix = "initial ";
        for (const auto &key_and_value : this->get_parameters().mapped_particle_properties)
          {
            const std::string &property_name = key_and_value.second.first;
            if (property_name.size() > prefix.size() &&
                property_name.compare(0, prefix.size(), prefix) == 0)
              data.push_back(this->get_initial_composition_manager().initial_composition(position, key_and_value.first));
          }
      }



      template <int dim>
      InitializationModeForLateParticles
      InitialComposition<dim>::late_initialization_mode () const
      {
        return interpolate_respect_boundary;
      }



      template <int dim>
      AdvectionField
      InitialComposition<dim>::advection_field_for_boundary_initialization (const unsigned int property_component) const
      {
        Assert (property_component < this->n_compositional_fields(),
                ExcInternalError());

        return AdvectionField::composition(property_component);
      }



      template <int dim>
      std::vector<std::pair<std::string, unsigned int>>
      InitialComposition<dim>::get_property_information() const
      {
        std::vector<std::pair<std::string,unsigned int>> property_information;

        // Find the compositional fields that are mapped to this particle property
        const std::string prefix = "initial ";
        for (const auto &key_and_value : this->get_parameters().mapped_particle_properties)
          {
            const std::string &property_name = key_and_value.second.first;
            if (property_name.size() > prefix.size() &&
                property_name.compare(0, prefix.size(), prefix) == 0)
              {
                AssertThrow(key_and_value.second.second == 0, ExcNotImplemented());
                property_information.emplace_back(key_and_value.second.first, 1);
              }
          }

        AssertThrow(property_information.size() > 0,
                    ExcMessage("You have requested the particle property <initial "
                               "composition>, but there is no compositional field "
                               "mapped to this particle property. Please prefix the "
                               "property names of the corresponding compositional fields "
                               "by 'initial ', or remove this particle property."));

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
      ASPECT_REGISTER_PARTICLE_PROPERTY(InitialComposition,
                                        "initial composition",
                                        "Implementation of a plugin in which the particle "
                                        "property is given as the initial composition "
                                        "at the particle's initial position. The particle "
                                        "gets as many properties as there are "
                                        "compositional fields.")
    }
  }
}
