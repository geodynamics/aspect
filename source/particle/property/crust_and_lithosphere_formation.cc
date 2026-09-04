/*
  Copyright (C) 2025 by the authors of the ASPECT code.

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

#include <aspect/particle/property/crust_and_lithosphere_formation.h>
#include <aspect/material_model/reaction_model/crust_and_lithosphere_formation.h>
#include <aspect/initial_composition/interface.h>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      template <int dim>
      CrustLithosphereFormation<dim>::CrustLithosphereFormation ()
        :
        material_inputs(0,0),
        material_outputs(0,0)
      {}



      template <int dim>
      void
      CrustLithosphereFormation<dim>::initialize ()
      {
        material_inputs.resize (1, this->n_compositional_fields());
        material_outputs.resize (1, this->n_compositional_fields(), true);

        AssertThrow(this->introspection().compositional_name_exists("basalt") &&
                    this->introspection().compositional_name_exists("harzburgite"),
                    ExcMessage("The particle property <crust and lithosphere formation> "
                               "can only be used if there are compositional fields named "
                               "'basalt' and 'harzburgite'."));

        basalt_index = this->introspection().compositional_index_for_name("basalt");
        harzburgite_index = this->introspection().compositional_index_for_name("harzburgite");
      }



      template <int dim>
      void
      CrustLithosphereFormation<dim>::initialize_one_particle_property(const Point<dim> &position,
                                                                       std::vector<double> &data) const
      {
        // Set the initial composition to the initial basalt and harzburgite fractions.
        data.push_back(this->get_initial_composition_manager().initial_composition(position,basalt_index));
        data.push_back(this->get_initial_composition_manager().initial_composition(position,harzburgite_index));

        if (if_track_basalt_formation_time)
          {
            // First basalt formation time. It is initialized to a value smaller than
            // the model start time so that particles that have not yet been converted
            // to basalt can be distinguished.
            data.push_back(this->get_parameters().start_time - 1e10);
            // Last basalt formation time. It is initialized to a value smaller than
            // the model start time so that particles that have not yet been converted
            // to basalt a second time can be distinguished.
            data.push_back(this->get_parameters().start_time - 1e10);
          }
      }



      template <int dim>
      void
      CrustLithosphereFormation<dim>::update_particle_properties(const ParticleUpdateInputs<dim> &inputs,
                                                                 typename ParticleHandler<dim>::particle_iterator_range &particles) const
      {
        material_inputs.resize (inputs.solution.size(), this->n_compositional_fields());
        material_outputs.resize (inputs.solution.size(), this->n_compositional_fields(), true);
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

            for (unsigned int n = 0; n < this->n_compositional_fields(); ++n)
              material_inputs.composition[p][n] = inputs.solution[p][this->introspection().component_indices.compositional_fields[n]];

            material_inputs.composition[p][basalt_index] = particle.get_properties()[this->data_position];
            material_inputs.composition[p][harzburgite_index] = particle.get_properties()[this->data_position+1];

            Tensor<2,dim> grad_u;
            for (unsigned int d=0; d<dim; ++d)
              grad_u[d] = inputs.gradients[p][this->introspection().component_indices.velocities[d]];

            material_inputs.strain_rate[p] = symmetrize (grad_u);

            // We assume that the only reactions come from the
            // crust and lithosphere formation reaction model.
            material_outputs.reaction_terms[p][basalt_index] = 0.0;
            material_outputs.reaction_terms[p][harzburgite_index] = 0.0;

            ++p;
          }

        crust_lithosphere_formation->calculate_reaction_terms(material_inputs,
                                                              material_outputs);

        p = 0;
        for (auto &particle: particles)
          {
            particle.get_properties()[this->data_position] += material_outputs.reaction_terms[p][basalt_index];
            particle.get_properties()[this->data_position+1] += material_outputs.reaction_terms[p][harzburgite_index];

            if (if_track_basalt_formation_time)
              {
                const double current_time = this->get_time() / (this->convert_output_to_years() ? year_in_seconds : 1.0);

                const Particle::Property::Manager<dim> &manager = this->get_particle_manager(this->get_particle_manager_index()).get_property_manager();
                const unsigned int generation_time_index = manager.get_data_info().get_field_index_by_name("particle generation time");

                // For newly generated particles, their initial values are interpolated from
                // the surrounding particles. We reset their basalt formation times to the lowest value
                // so that they can be identified as not yet converted to basalt.
                if (particle.get_properties()[generation_time_index] == current_time)
                  {
                    particle.get_properties()[this->data_position + 2] = this->get_parameters().start_time - 1e10;
                    particle.get_properties()[this->data_position + 3] = this->get_parameters().start_time - 1e10;
                  }

                // If the particle starts as basalt or has already been converted to basalt,
                // the basalt reaction term remains zero, even if it satisfies the basalt formation condition.
                // Therefore, we use a fake "zero basalt" value to check whether the particle
                // would satisfy the condition for new basalt formation.
                material_inputs.composition[p][basalt_index] = 0;
                crust_lithosphere_formation->calculate_reaction_terms(material_inputs,material_outputs);

                if (material_outputs.reaction_terms[p][basalt_index] > 0)
                  {
                    // First generation:
                    // The first basalt formation time is initialized to the lowest value.
                    // A value smaller than the model start time indicates that the particle
                    // has not yet been converted to basalt. Newly generated particles are not converted
                    // to basalt in the same time step, because the basalt formation time is intended to
                    // record the first time the particle reaches the conversion condition.
                    if ((particle.get_properties()[this->data_position+2] < this->get_parameters().start_time) && (particle.get_properties()[generation_time_index] < current_time))
                      particle.get_properties()[this->data_position+2] = current_time;

                    // Last generation:
                    // We record the last generation time if the particle has been converted to basalt at least once before.
                    else if (particle.get_properties()[this->data_position+2] > this->get_parameters().start_time)
                      {
                        particle.get_properties()[this->data_position+3] = current_time;
                      }
                  }
              }
            ++p;
          }
      }



      template <int dim>
      InitializationModeForLateParticles
      CrustLithosphereFormation<dim>::late_initialization_mode () const
      {
        return interpolate_respect_boundary;
      }



      template <int dim>
      AdvectionField
      CrustLithosphereFormation<dim>::advection_field_for_boundary_initialization (const unsigned int property_component) const
      {
        if (property_component == 0)
          return AdvectionField::composition(basalt_index);
        else if (property_component == 1)
          return AdvectionField::composition(harzburgite_index);
        else
          {
            Assert(false,
                   ExcMessage("The crust and lithosphere formation particle property "
                              "only has two components: basalt and harzburgite."));
            // This line will never be reached but is needed to avoid compiler warnings.
            return AdvectionField::composition(0);
          }
      }



      template <int dim>
      UpdateTimeFlags
      CrustLithosphereFormation<dim>::need_update() const
      {
        return update_time_step;
      }



      template <int dim>
      UpdateFlags
      CrustLithosphereFormation<dim>::get_needed_update_flags () const
      {
        return update_values | update_gradients;
      }



      template <int dim>
      std::vector<std::pair<std::string, unsigned int>>
      CrustLithosphereFormation<dim>::get_property_information() const
      {
        if (if_track_basalt_formation_time)
          return {{"basalt",1},
            {"harzburgite",1},
            {"first_basalt_formation_time",1},
            {"last_basalt_formation_time",1}
          };

        else
          return {{"basalt",1}, {"harzburgite",1}};
      }



      template <int dim>
      void
      CrustLithosphereFormation<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Crust and lithosphere formation");
        {
          // Add parameter for 'track basalt formation' on or off
          MaterialModel::ReactionModel::CrustLithosphereFormation<dim>::declare_parameters(prm);

          prm.declare_entry ("Track basalt formation time", "false",
                             Patterns::Bool (),
                             "If true, the model time when a particle is first and last converted to basalt "
                             "will be recorded as particle properties. A value smaller than the model start time"
                             "indicates that the particle has not yet been converted to basalt.");

        }
        prm.leave_subsection();
      }


      template <int dim>
      void
      CrustLithosphereFormation<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Crust and lithosphere formation");
        {
          // Parse reaction model parameters
          crust_lithosphere_formation = std::make_unique<MaterialModel::ReactionModel::CrustLithosphereFormation<dim>>();
          crust_lithosphere_formation->initialize_simulator(this->get_simulator());
          crust_lithosphere_formation->parse_parameters(prm);

          if_track_basalt_formation_time = prm.get_bool("Track basalt formation time");

          // When tracking basalt formation times, we need to make sure that
          // the particle get a correct initialization time when they are created.
          // The particle generation time property is required to identify newly generated particles
          const auto &manager = this->get_particle_manager(this->get_particle_manager_index()).get_property_manager();
          AssertThrow(if_track_basalt_formation_time == false || manager.plugin_name_exists("particle generation time"),
                      ExcMessage("The 'particle generation time' plugin must be enabled to track "
                                 "the basalt formation time."));
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
      ASPECT_REGISTER_PARTICLE_PROPERTY(CrustLithosphereFormation,
                                        "crust and lithosphere formation",
                                        "A plugin in which the particle property is "
                                        "defined as the evolving chemical composition "
                                        "that results from the formation of oceanic crust "
                                        "and lithosphere as mantle material approaches the "
                                        "surface and melts. Note that this does not necessarily "
                                        "conserves the bulk chemical composition of the mantle, "
                                        "since the conversion only depends on the mantle flow. "
                                        "See the crust and lithosphere formation reaction model "
                                        "documentation for more detailed information.")

    }
  }
}
