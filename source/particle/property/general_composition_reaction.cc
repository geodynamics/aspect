/*
  Copyright (C) 2022 - 2024 by the authors of the ASPECT code.

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

#include <aspect/particle/property/general_composition_reaction.h>
#include <aspect/initial_composition/interface.h>
#include <aspect/material_model/interface.h>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      template <int dim>
      GeneralCompositionReaction<dim>::GeneralCompositionReaction ()
        :
        material_inputs(0,0),
        material_outputs(0,0)
      {}



      template <int dim>
      void
      GeneralCompositionReaction<dim>::initialize_one_particle_property(const Point<dim> &position,
                                                                        std::vector<double> &data) const
      {
        // Set the initial composition to selected compositional fields
        for (unsigned int composition_index : selected_compositional_field_indices)
          {
            data.push_back(this->get_initial_composition_manager().initial_composition(position, composition_index));
          }
      }



      template <int dim>
      void
      GeneralCompositionReaction<dim>::update_particle_properties(const ParticleUpdateInputs<dim> &inputs,
                                                                  typename ParticleHandler<dim>::particle_iterator_range &particles) const
      {
        material_inputs.resize(inputs.solution.size(), this->n_compositional_fields());
        material_outputs.resize(inputs.solution.size(), this->n_compositional_fields());
        material_inputs.requested_properties = MaterialModel::MaterialProperties::reaction_terms;
        material_inputs.current_cell = inputs.current_cell;

        unsigned int p = 0;
        for (auto particle: particles)
          {
            // Make sure all particles are in the same cell
            Assert(particle.get_surrounding_cell() == inputs.current_cell,
                   ExcMessage("All particles must be in the same cell."));

            material_inputs.position[p] = particle.get_location();
            material_inputs.temperature[p] = inputs.solution[p][this->introspection().component_indices.temperature];
            material_inputs.pressure[p] = inputs.solution[p][this->introspection().component_indices.pressure];

            for (unsigned int d = 0; d < dim; ++d)
              material_inputs.velocity[p][d] = inputs.solution[p][this->introspection().component_indices.velocities[d]];

            for (unsigned int i=0; i<this->n_compositional_fields(); ++i)
              material_inputs.composition[p][i] = inputs.solution[p][this->introspection().component_indices.compositional_fields[i]];

            // Loop through all selected compositional fields in this particle manager and overwrite them with the particle properties
            unsigned int local_index=0;
            for (unsigned int composition_index : selected_compositional_field_indices)
              {
                material_inputs.composition[p][composition_index] = particle.get_properties()[this->data_position + local_index];
                ++local_index;
              }

            Tensor<2,dim> grad_u;
            for (unsigned int d=0; d<dim; ++d)
              grad_u[d] = inputs.gradients[p][this->introspection().component_indices.velocities[d]];

            material_inputs.strain_rate[p] = symmetrize (grad_u);

            ++p;
          }

        // call the material model to evaluate the reaction terms
        this->get_material_model().evaluate(material_inputs,
                                            material_outputs);

        // update the particle properties with the reaction terms
        p = 0;
        for (auto &particle: particles)
          {
            unsigned int local_index = 0;
            // update the particle properties with the reaction terms for the selected compositional fields in this particle manager
            for (unsigned int composition_index : selected_compositional_field_indices)
              {
                particle.get_properties()[this->data_position + local_index] += material_outputs.reaction_terms[p][composition_index];
                ++local_index;
              }
            ++p;
          }
      }



      template <int dim>
      InitializationModeForLateParticles
      GeneralCompositionReaction<dim>::late_initialization_mode () const
      {
        return interpolate_respect_boundary;
      }



      template <int dim>
      AdvectionField
      GeneralCompositionReaction<dim>::advection_field_for_boundary_initialization (const unsigned int property_component) const
      {
        (void) property_component;
        Assert (property_component < selected_compositional_field_indices.size(),
                ExcInternalError());

        return AdvectionField::composition(selected_compositional_field_indices[property_component]);
      }



      template <int dim>
      UpdateTimeFlags
      GeneralCompositionReaction<dim>::need_update() const
      {
        return update_time_step;
      }



      template <int dim>
      UpdateFlags
      GeneralCompositionReaction<dim>::get_update_flags (const unsigned int component) const
      {
        if (this->introspection().component_masks.velocities[component] == true)
          return update_values | update_gradients;

        return update_values;
      }



      template <int dim>
      std::vector<std::pair<std::string, unsigned int>>
      GeneralCompositionReaction<dim>::get_property_information() const
      {
        std::vector<std::pair<std::string, unsigned int>>
        property_information;

        for (const unsigned int composition_index :
             selected_compositional_field_indices)
          {
            const std::string field_name = this->introspection().name_for_compositional_index(composition_index);
            property_information.emplace_back(field_name,1);
          }

        return property_information;
      }

      template <int dim>
      void
      GeneralCompositionReaction<dim>::declare_parameters(ParameterHandler &prm)
      {
        prm.enter_subsection("General composition reaction");
        {
          prm.declare_entry(
            "Selected compositional fields",
            "all",
            Patterns::List(Patterns::Anything()),
            "A list that deterines which Compositional fields are stored for particles of this kind, "
            "and on this particle manager. The value 'all' selects every compositional field.");
        }
        prm.leave_subsection();
      }

      template<int dim>
      void GeneralCompositionReaction<dim>::parse_parameters(
        ParameterHandler &prm)
      {
        prm.enter_subsection("General composition reaction");
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
      ASPECT_REGISTER_PARTICLE_PROPERTY(GeneralCompositionReaction,
                                        "general composition reaction",
                                        "A particle property that stores compositional fields "
                                        "and updates them using the 'reaction terms' returned "
                                        "by the material model. While this particle property "
                                        "allows fields to be tracked on different particle "
                                        "managers, for now it is recommended to have only one "
                                        "particle manager with this property in your simulation.");
    }
  }
}
