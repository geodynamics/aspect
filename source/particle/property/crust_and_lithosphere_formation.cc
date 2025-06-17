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
        material_inputs(1,0),
        material_outputs(1,0)
      {}



      template <int dim>
      void
      CrustLithosphereFormation<dim>::initialize ()
      {
        material_inputs  = MaterialModel::MaterialModelInputs<dim>(1, this->n_compositional_fields());
        material_outputs = MaterialModel::MaterialModelOutputs<dim>(1, this->n_compositional_fields());

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
      }



      template <int dim>
      void
      CrustLithosphereFormation<dim>::update_particle_properties(const ParticleUpdateInputs<dim> &inputs,
                                                                 typename ParticleHandler<dim>::particle_iterator_range &particles) const
      {
        material_inputs  = MaterialModel::MaterialModelInputs<dim>(inputs.solution.size(), this->n_compositional_fields());
        material_outputs = MaterialModel::MaterialModelOutputs<dim>(inputs.solution.size(), this->n_compositional_fields());
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
        return {{"basalt",1}, {"harzburgite",1}};
      }



      template <int dim>
      void
      CrustLithosphereFormation<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Crust and lithosphere formation");
        {
          MaterialModel::ReactionModel::CrustLithosphereFormation<dim>::declare_parameters(prm);
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
