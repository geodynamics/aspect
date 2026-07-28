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
#include <aspect/prescribed_solution/interface.h>

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
        for (unsigned int i = 0; i < this->n_compositional_fields(); ++i)
          data.push_back(this->get_initial_composition_manager().initial_composition(position,i));
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
        AssertThrow(this->n_compositional_fields() > 0,
                    ExcMessage("You have requested the particle property <initial "
                               "composition>, but the number of compositional fields is 0. "
                               "Please add compositional fields to your model, or remove "
                               "this particle property."));

        std::vector<std::pair<std::string,unsigned int>> property_information;

        for (unsigned int i = 0; i < this->n_compositional_fields(); ++i)
          {
            std::ostringstream field_name;
            field_name << "initial " << this->introspection().name_for_compositional_index(i);
            property_information.emplace_back(field_name.str(),1);
          }

        return property_information;
      }

      template <int dim>
      void
      InitialComposition<dim>::update_particle_properties(const ParticleUpdateInputs<dim> &/*inputs*/,
                                                          typename ParticleHandler<dim>::particle_iterator_range &particles) const
      {
        const aspect::PrescribedSolution::Manager<dim> &prescribed_solution_manager = this->get_prescribed_solution();
        const auto &plugin_objects = prescribed_solution_manager.get_active_plugins();


        std::vector<Point<dim>> evaluation_points(1);
        std::vector<unsigned int> component_indices(1);
        std::vector<bool> component_is_constrained(1);
        std::vector<double> constrained_component_value(1);

        for (auto &particle : particles)
          {
            evaluation_points[0] = particle.get_location();

            for (unsigned int j = 0; j < this->n_compositional_fields(); ++j)
              {
                component_indices[0] = this->introspection().component_indices.compositional_fields[j];
                component_is_constrained[0] = false;
                constrained_component_value[0] = 0.0;

                for (const auto &plugin : plugin_objects)
                  {
                    typename DoFHandler<dim>::active_cell_iterator cell;
                    plugin->constrain_solution(cell, evaluation_points, component_indices, component_is_constrained, constrained_component_value);
                  }

                if (component_is_constrained[0])
                  particle.get_properties()[this->data_position+j] = constrained_component_value[0];
              }
          }
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
