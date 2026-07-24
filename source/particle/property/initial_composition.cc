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
        for (const unsigned int composition_index : selected_compositional_field_indices)
          data.push_back(this->get_initial_composition_manager().initial_composition(position,composition_index));
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
        Assert (property_component < selected_compositional_field_indices.size(),
                ExcInternalError());
        return AdvectionField::composition(selected_compositional_field_indices[property_component]);
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

        for (const unsigned int composition_index : selected_compositional_field_indices)
          {
            std::ostringstream field_name;
            field_name << "initial " << this->introspection().name_for_compositional_index(composition_index);
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


            unsigned int local_index=0;
            for (const unsigned int composition_index : selected_compositional_field_indices)
              {
                component_indices[0] = this->introspection().component_indices.compositional_fields[composition_index];
                component_is_constrained[0] = false;
                constrained_component_value[0] = 0.0;

                for (const auto &plugin : plugin_objects)
                  {
                    typename DoFHandler<dim>::active_cell_iterator cell;
                    plugin->constrain_solution(cell, evaluation_points, component_indices, component_is_constrained, constrained_component_value);
                  }

                if (component_is_constrained[0])
                  particle.get_properties()[this->data_position+local_index] = constrained_component_value[0];
                local_index++;
              }
          }
      }

      template<int dim>
      void InitialComposition<dim>::declare_parameters(ParameterHandler &prm)
      {
        prm.enter_subsection("Initial composition");
        {
          prm.declare_entry("Selected compositional fields",
                            "all",
                            Patterns::List(Patterns::Anything()),
                            "A list that determines which compositional fields to "
                            "track on this particle manager. The default value 'all' means "
                            "every composition field is tracked on every particle manager." );
        }
        prm.leave_subsection();
      }



      template<int dim>
      void InitialComposition<dim>::parse_parameters(ParameterHandler &prm)
      {
        prm.enter_subsection("Initial composition");
        {
          if (prm.get("Selected compositional fields") == "all")
            {
              for (unsigned int i=0; i<this->n_compositional_fields(); ++i)
                {
                  selected_compositional_field_indices.push_back(i);
                }
            }
          else
            {
              const std::vector<std::string> selected_fields = Utilities::split_string_list(prm.get("Selected compositional fields"));
              for (const std::string &selected_field : selected_fields)
                {
                  selected_compositional_field_indices.push_back(this->introspection().compositional_index_for_name(selected_field));
                }
            }
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
      ASPECT_REGISTER_PARTICLE_PROPERTY(InitialComposition,
                                        "initial composition",
                                        "Implementation of a plugin in which the particle "
                                        "property is given as the initial composition "
                                        "at the particle's initial position. The "
                                        "'Selected compositional fields' chooses which "
                                        "compositional fields to track on this particle "
                                        "manager. If no 'Selected compositional fields' "
                                        "are chosen, the particle manager gets as many "
                                        "properties as there are compositional fields.")
    }
  }
}
