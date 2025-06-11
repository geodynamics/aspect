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

#include <aspect/particle/property/elastic_stress.h>
#include <aspect/material_model/visco_plastic.h>
#include <aspect/material_model/viscoelastic.h>
#include <aspect/initial_composition/interface.h>
#include <aspect/particle/world.h>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      template <int dim>
      ElasticStress<dim>::ElasticStress ()
        :
        material_inputs(1,0),
        material_outputs(1,0),
        material_inputs_cell(1,0),
        material_outputs_cell(1,0)
      {}



      template <int dim>
      void
      ElasticStress<dim>::initialize ()
      {
        AssertThrow((Plugins::plugin_type_matches<const MaterialModel::ViscoPlastic<dim>>(this->get_material_model())
                     ||
                     Plugins::plugin_type_matches<const MaterialModel::Viscoelastic<dim>>(this->get_material_model())),
                    ExcMessage("This particle property only makes sense in combination with the viscoelastic or visco_plastic material model."));

        AssertThrow(this->get_parameters().enable_elasticity == true,
                    ExcMessage ("This particle property should only be used if 'Enable elasticity' is set to true"));

        const auto &manager = this->get_particle_manager(this->get_particle_manager_index()).get_property_manager();
        AssertThrow(!manager.plugin_name_exists("composition"),
                    ExcMessage("The 'elastic stress' plugin cannot be used in combination with the 'composition' plugin."));

        material_inputs = MaterialModel::MaterialModelInputs<dim>(1, this->n_compositional_fields());

        material_outputs = MaterialModel::MaterialModelOutputs<dim>(1, this->n_compositional_fields());

        material_inputs.requested_properties = MaterialModel::MaterialProperties::reaction_terms;

        // Get the indices of those compositions that correspond to stress tensor elements.
        stress_field_indices = this->introspection().get_indices_for_fields_of_type(CompositionalFieldDescription::stress);
        AssertThrow((stress_field_indices.size() == 2*SymmetricTensor<2,dim>::n_independent_components),
                    ExcMessage("The number of stress tensor element fields in the 'elastic stress' plugin does not equal twice the number of independent components."));

        // Get the indices of all compositions that do not correspond to stress tensor elements.
        std::vector<unsigned int> all_field_indices(this->n_compositional_fields());
        std::iota (std::begin(all_field_indices), std::end(all_field_indices), 0);
        std::set_difference(all_field_indices.begin(), all_field_indices.end(),
                            stress_field_indices.begin(), stress_field_indices.end(),
                            std::inserter(non_stress_field_indices, non_stress_field_indices.begin()));

        // Connect to the signal after particles are restored at the beginning of
        // a nonlinear iteration of iterative advection schemes.
        this->get_signals().post_restore_particles.connect(
          [&](typename Particle::Manager<dim> &particle_manager)
        {
          this->update_particles(particle_manager);
        }
        );
      }



      template <int dim>
      void
      ElasticStress<dim>::update_particles(typename Particle::Manager<dim> &particle_manager) const
      {
        // There is no update of the stress to apply during the first (0th) timestep
        if (this->simulator_is_past_initialization() == false || this->get_timestep_number() == 0)
          return;

        // Determine the data position of the first stress tensor component
        const unsigned int data_position = particle_manager.get_property_manager().get_data_info().get_position_by_field_name("ve_stress_xx");

        // Get handler
        Particle::ParticleHandler<dim> &particle_handler = particle_manager.get_particle_handler();

        // TODO instead of calling the manager, and looping over all properties,
        // can we use this property's get_update_flags() function only?
        const std::vector<UpdateFlags> update_flags = particle_manager.get_property_manager().get_update_flags();

        // combine all update flags to a single flag, which is the required information
        // for the mapping inside the solution evaluator
        UpdateFlags mapping_flags = update_flags[0];
        for (unsigned int i=1; i<update_flags.size(); ++i)
          mapping_flags |= update_flags[i];

        // Create evaluators that get the old solution (= solution of previous timestep)
        // at the locations of the particles within one cell.
        // The particles have not been advected in the current timestep yet, or have been restored
        // to their pre-advection locations, so they are in their old locations corresponding to the old
        // solution.
        std::unique_ptr<SolutionEvaluator<dim>> evaluators = construct_solution_evaluator(*this,
                                                              mapping_flags);

        // FEPointEvaluation uses different evaluation flags than the common UpdateFlags.
        // Translate between the two.
        std::vector<EvaluationFlags::EvaluationFlags> evaluation_flags (update_flags.size(), EvaluationFlags::nothing);

        for (unsigned int i=0; i<update_flags.size(); ++i)
          {
            if (update_flags[i] & update_values)
              evaluation_flags[i] |= EvaluationFlags::values;

            if (update_flags[i] & update_gradients)
              evaluation_flags[i] |= EvaluationFlags::gradients;
          }

        // Vector to store the positions of all the particles in one cell
        small_vector<Point<dim>> positions;

        // Loop over all active and owned cells
        for (const auto &cell : this->get_dof_handler().active_cell_iterators())
          if (cell->is_locally_owned())
            {
              // Find which particles are in the current cell
              typename ParticleHandler<dim>::particle_iterator_range
              particles_in_cell = particle_handler.particles_in_cell(cell);

              // Only update particles, if there are any in this cell
              if (particles_in_cell.begin() != particles_in_cell.end())
                {
                  const unsigned int n_particles_in_cell = particle_handler.n_particles_in_cell(cell);
                  // Get the locations of the particles within the reference cell
                  positions.resize(n_particles_in_cell);
                  unsigned int p = 0;
                  for (const auto &particle : particles_in_cell)
                    {
                      positions[p] = particle.get_reference_location();
                      ++p;
                    }

                  // Resize the material model inputs to the number of particles in the current cell
                  material_inputs_cell  = MaterialModel::MaterialModelInputs<dim>(n_particles_in_cell, this->n_compositional_fields());
                  material_inputs_cell.current_cell = cell;
                  material_inputs_cell.requested_properties = MaterialModel::MaterialProperties::reaction_rates;
                  material_outputs_cell = MaterialModel::MaterialModelOutputs<dim>(n_particles_in_cell, this->n_compositional_fields());
                  // The reaction rates are stored in additional outputs
                  this->get_material_model().create_additional_named_outputs(material_outputs_cell);

                  const std::shared_ptr<MaterialModel::ReactionRateOutputs<dim>> reaction_rate_outputs
                    = material_outputs_cell.template get_additional_output_object<MaterialModel::ReactionRateOutputs<dim>>();

                  // Collect the values of the old solution restricted to the current cell's DOFs
                  small_vector<double> old_solution_values(this->get_fe().dofs_per_cell);
                  cell->get_dof_values(this->get_old_solution(),
                                       old_solution_values.begin(),
                                       old_solution_values.end());

                  EvaluationFlags::EvaluationFlags evaluation_flags_union = EvaluationFlags::nothing;
                  for (unsigned int i=0; i<evaluation_flags.size(); ++i)
                    evaluation_flags_union |= evaluation_flags[i];

                  // Update evaluators to the current cell
                  if (evaluation_flags_union & (EvaluationFlags::values | EvaluationFlags::gradients))
                    {
                      // Reinitialize and evaluate the requested solution values and gradients
                      evaluators->reinit(cell, {positions.data(), positions.size()});

                      evaluators->evaluate({old_solution_values.data(),old_solution_values.size()},
                                           evaluation_flags);
                    }

                  // To store the old solutions
                  std::vector<small_vector<double,50>> old_solution(n_particles_in_cell,small_vector<double,50>(evaluators->n_components(), numbers::signaling_nan<double>()));

                  // To store the old gradients
                  std::vector<small_vector<Tensor<1,dim>,50>> old_gradients(n_particles_in_cell,small_vector<Tensor<1,dim>,50>(evaluators->n_components(), numbers::signaling_nan<Tensor<1,dim>>()));

                  // Loop over all particles in the cell
                  auto particle = particles_in_cell.begin();
                  for (unsigned int i = 0; particle!=particles_in_cell.end(); ++particle,++i)
                    {
                      // Evaluate the old solution, but only if it is requested in the update_flags
                      if (evaluation_flags_union & EvaluationFlags::values)
                        evaluators->get_solution(i, {&old_solution[i][0],old_solution[i].size()}, evaluation_flags);

                      // Evaluate the old gradients, but only if they are requested in the update_flags
                      if (evaluation_flags_union & EvaluationFlags::gradients)
                        evaluators->get_gradients(i, {&old_gradients[i][0],old_gradients[i].size()}, evaluation_flags);

                      // Fill material model input
                      // Get the real location of the particle
                      material_inputs_cell.position[i] = particle->get_location();

                      material_inputs_cell.temperature[i] = old_solution[i][this->introspection().component_indices.temperature];

                      material_inputs_cell.pressure[i] = old_solution[i][this->introspection().component_indices.pressure];

                      for (unsigned int d = 0; d < dim; ++d)
                        material_inputs_cell.velocity[i][d] = old_solution[i][this->introspection().component_indices.velocities[d]];

                      // Fill the non-stress composition inputs with the old_solution.
                      for (const unsigned int &n : non_stress_field_indices)
                        material_inputs_cell.composition[i][n] = old_solution[i][this->introspection().component_indices.compositional_fields[n]];

                      // Retrieve the ve_stress_* values from the particles and fields and
                      // fill the material model inputs with a weighted average of the two.
                      // In some cases, using the field value leads to more stable results.
                      // After the particles have been restored, their properties have the values of the previous timestep.
                      for (unsigned int n = 0; n < 2*SymmetricTensor<2,dim>::n_independent_components; ++n)
                        {
                          const double particle_stress_value = particle->get_properties()[data_position + n];
                          const double field_stress_value = old_solution[i][this->introspection().component_indices.compositional_fields[stress_field_indices[n]]];
                          const double stress_value = particle_weight * particle_stress_value + (1-particle_weight) * field_stress_value;
                          material_inputs_cell.composition[i][stress_field_indices[n]] = stress_value;
                        }

                      Tensor<2,dim> grad_u;
                      for (unsigned int d=0; d<dim; ++d)
                        grad_u[d] = old_gradients[i][d];
                      material_inputs_cell.strain_rate[i] = symmetrize (grad_u);
                    }

                  // Evaluate the material model to get the reaction rates
                  // for all the particles in the current cell.
                  this->get_material_model().evaluate (material_inputs_cell,material_outputs_cell);

                  // Update all particles in the current cell.
                  particle = particles_in_cell.begin();
                  for (unsigned int i = 0; particle!=particles_in_cell.end(); ++particle,++i)
                    {
                      // Add the reaction_rates * timestep = update to the corresponding stress
                      // tensor components
                      for (unsigned int p = 0; p < 2*SymmetricTensor<2,dim>::n_independent_components ; ++p)
                        particle->get_properties()[data_position + p] = material_inputs_cell.composition[i][stress_field_indices[p]] + reaction_rate_outputs->reaction_rates[i][stress_field_indices[p]] * this->get_timestep();
                    }
                }
            }
      }



      template <int dim>
      void
      ElasticStress<dim>::initialize_one_particle_property(const Point<dim> &position,
                                                           std::vector<double> &data) const
      {
        // Give each elastic stress field its initial composition if one is prescribed.
        data.push_back(this->get_initial_composition_manager().initial_composition(position,this->introspection().compositional_index_for_name("ve_stress_xx")));

        data.push_back(this->get_initial_composition_manager().initial_composition(position,this->introspection().compositional_index_for_name("ve_stress_yy")));

        if (dim == 2)
          {
            data.push_back(this->get_initial_composition_manager().initial_composition(position,this->introspection().compositional_index_for_name("ve_stress_xy")));

          }
        else if (dim == 3)
          {
            data.push_back(this->get_initial_composition_manager().initial_composition(position,this->introspection().compositional_index_for_name("ve_stress_zz")));

            data.push_back(this->get_initial_composition_manager().initial_composition(position,this->introspection().compositional_index_for_name("ve_stress_xy")));

            data.push_back(this->get_initial_composition_manager().initial_composition(position,this->introspection().compositional_index_for_name("ve_stress_xz")));

            data.push_back(this->get_initial_composition_manager().initial_composition(position,this->introspection().compositional_index_for_name("ve_stress_yz")));
          }

        data.push_back(this->get_initial_composition_manager().initial_composition(position,this->introspection().compositional_index_for_name("ve_stress_xx_old")));

        data.push_back(this->get_initial_composition_manager().initial_composition(position,this->introspection().compositional_index_for_name("ve_stress_yy_old")));

        if (dim == 2)
          {
            data.push_back(this->get_initial_composition_manager().initial_composition(position,this->introspection().compositional_index_for_name("ve_stress_xy_old")));
          }
        else if (dim == 3)
          {
            data.push_back(this->get_initial_composition_manager().initial_composition(position,this->introspection().compositional_index_for_name("ve_stress_zz_old")));

            data.push_back(this->get_initial_composition_manager().initial_composition(position,this->introspection().compositional_index_for_name("ve_stress_xy_old")));

            data.push_back(this->get_initial_composition_manager().initial_composition(position,this->introspection().compositional_index_for_name("ve_stress_xz_old")));

            data.push_back(this->get_initial_composition_manager().initial_composition(position,this->introspection().compositional_index_for_name("ve_stress_yz_old")));
          }
      }



      template <int dim>
      void
      ElasticStress<dim>::update_particle_properties(const ParticleUpdateInputs<dim> &inputs,
                                                     typename ParticleHandler<dim>::particle_iterator_range &particles) const
      {
        unsigned int p = 0;
        for (auto &particle: particles)
          {
            material_inputs.position[0] = particle.get_location();


            material_inputs.current_cell = inputs.current_cell;

            material_inputs.temperature[0] = inputs.solution[p][this->introspection().component_indices.temperature];

            material_inputs.pressure[0] = inputs.solution[p][this->introspection().component_indices.pressure];

            for (unsigned int d = 0; d < dim; ++d)
              material_inputs.velocity[0][d] = inputs.solution[p][this->introspection().component_indices.velocities[d]];

            // Fill the non-stress composition inputs with the solution.
            for (const unsigned int &n : non_stress_field_indices)
              material_inputs.composition[0][n] = inputs.solution[p][this->introspection().component_indices.compositional_fields[n]];
            // For the stress composition we use the ve_stress_* stored on the particles.
            for (unsigned int n = 0; n < 2*SymmetricTensor<2,dim>::n_independent_components; ++n)
              material_inputs.composition[0][stress_field_indices[n]] = particle.get_properties()[this->data_position + n];

            Tensor<2,dim> grad_u;
            for (unsigned int d=0; d<dim; ++d)
              grad_u[d] = inputs.gradients[p][d];
            material_inputs.strain_rate[0] = symmetrize (grad_u);

            this->get_material_model().evaluate (material_inputs,material_outputs);

            // Apply the stress rotation to the ve_stress_* fields, not the ve_stress_*_old fields.
            for (unsigned int i = 0; i < SymmetricTensor<2,dim>::n_independent_components ; ++i)
              particle.get_properties()[this->data_position + i] += material_outputs.reaction_terms[0][stress_field_indices[i]];

            ++p;
          }
      }



      template <int dim>
      UpdateTimeFlags
      ElasticStress<dim>::need_update() const
      {
        return update_time_step;
      }



      template <int dim>
      UpdateFlags
      ElasticStress<dim>::get_update_flags (const unsigned int component) const
      {
        if (this->introspection().component_masks.velocities[component] == true)
          return update_values | update_gradients;

        return update_values;
      }



      template <int dim>
      std::vector<std::pair<std::string, unsigned int>>
      ElasticStress<dim>::get_property_information() const
      {
        std::vector<std::pair<std::string,unsigned int>> property_information;

        property_information.emplace_back("ve_stress_xx",1);
        property_information.emplace_back("ve_stress_yy",1);

        if (dim == 2)
          {
            property_information.emplace_back("ve_stress_xy",1);
          }
        else if (dim == 3)
          {
            property_information.emplace_back("ve_stress_zz",1);
            property_information.emplace_back("ve_stress_xy",1);
            property_information.emplace_back("ve_stress_xz",1);
            property_information.emplace_back("ve_stress_yz",1);
          }

        property_information.emplace_back("ve_stress_xx_old",1);
        property_information.emplace_back("ve_stress_yy_old",1);

        if (dim == 2)
          {
            property_information.emplace_back("ve_stress_xy_old",1);
          }
        else if (dim == 3)
          {
            property_information.emplace_back("ve_stress_zz_old",1);
            property_information.emplace_back("ve_stress_xy_old",1);
            property_information.emplace_back("ve_stress_xz_old",1);
            property_information.emplace_back("ve_stress_yz_old",1);
          }

        return property_information;
      }



      template <int dim>
      void
      ElasticStress<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Elastic stress");
        {
          prm.declare_entry ("Particle stress value weight", "1.0",
                             Patterns::Double(0.),
                             "The weight given to the value of the stress tensor components "
                             "stored on the particles in the weighted average of those "
                             "values and the values of the compositional fields evaluated on the "
                             "particle location. The average is used in the Material Model inputs "
                             "used to compute the reaction rates and to update the particle property "
                             "with the reaction rates. In some cases, using the field values "
                             "leads to more stable results.");

        }
        prm.leave_subsection();
      }



      template <int dim>
      void
      ElasticStress<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Elastic stress");
        {
          particle_weight = prm.get_double("Particle stress value weight");
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
      ASPECT_REGISTER_PARTICLE_PROPERTY(ElasticStress,
                                        "elastic stress",
                                        "A plugin in which the particle property tensor is "
                                        "defined as the total elastic stress a particle has "
                                        "accumulated. This plugin modifies the properties "
                                        "with the name ve_stress_*. It first applies the stress "
                                        "change resulting from system evolution during the previous "
                                        "computational timestep, and then the rotation of those "
                                        "stresses into the current timestep. "
                                        "See the viscoelastic or visco_plastic material model "
                                        "documentation for more detailed information.")
    }
  }
}
