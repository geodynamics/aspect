/*
  Copyright (C) 2015 - 2022 by the authors of the ASPECT code.

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
        material_outputs(1,0)
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

        const auto &manager = this->get_particle_world().get_property_manager();
        AssertThrow(!manager.plugin_name_exists("composition"),
                    ExcMessage("The 'elastic stress' plugin cannot be used in combination with the 'composition' plugin."));

        material_inputs = MaterialModel::MaterialModelInputs<dim>(1, this->n_compositional_fields());

        material_outputs = MaterialModel::MaterialModelOutputs<dim>(1, this->n_compositional_fields());

        // The reaction rates are stored in additional outputs
        this->get_material_model().create_additional_named_outputs(material_outputs);

        // Connect to the signal after particles are restored at the beginning of
        // a nonlinear iteration of iterative advection schemes.
        this->get_signals().post_restore_particles.connect(
          [&](typename Particle::World<dim> &particle_world)
        {
          this->update_particles(particle_world);
        }
        );
      }



      template <int dim>
      void
      ElasticStress<dim>::update_particles(typename Particle::World<dim> &particle_world) const
      {
        // There is no update of the stress to apply during the first (0th) timestep
        if (this->simulator_is_past_initialization() == false || this->get_timestep_number() == 0)
          return;

        // Determine the data position of the first stress tensor component
        const unsigned int data_position = particle_world.get_property_manager().get_data_info().get_position_by_field_name("ve_stress_xx");

        // Get handler
        Particle::ParticleHandler<dim> &particle_handler = particle_world.get_particle_handler();

        // Structure for the reaction rates
        MaterialModel::ReactionRateOutputs<dim> *reaction_rate_outputs
          = material_outputs.template get_additional_output<MaterialModel::ReactionRateOutputs<dim>>();

        const UpdateFlags update_flags = get_needed_update_flags();

        // Create evaluators that get the old solution (= solution of previous timestep)
        // at the locations of the particles within one cell.
        // The particles have not been advected in the current timestep yet, or have been restored
        // to their pre-advection locations, so they are in their old locations corresponding to the old
        // solution.

        // Only use deal.II FEPointEvaluation if its fast path is used
        bool use_fast_path = false;
        if (dynamic_cast<const MappingQGeneric<dim> *>(&this->get_mapping()) != nullptr ||
            dynamic_cast<const MappingCartesian<dim> *>(&this->get_mapping()) != nullptr)
          use_fast_path = true;

        // For fast path only
        std::unique_ptr<internal::SolutionEvaluators<dim>> evaluators;

        if (use_fast_path == true)
          evaluators = internal::construct_solution_evaluators(*this,
                                                               update_flags);

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

                  if (use_fast_path == true)
                    {
                      // Get the locations of the particles within the reference cell
                      std::vector<Point<dim>> positions;
                      positions.reserve(n_particles_in_cell);
                      for (auto particle = particles_in_cell.begin(); particle!=particles_in_cell.end(); ++particle)
                        positions.push_back(particle->get_reference_location());


                      // Collect the values of the old solution restricted to the current cell's DOFs
                      boost::container::small_vector<double, 100> old_solution_values(this->get_fe().dofs_per_cell);
                      cell->get_dof_values(this->get_old_solution(),
                                           old_solution_values.begin(),
                                           old_solution_values.end());

                      // Update evaluators to the current cell
                      if (update_flags & (update_values | update_gradients))
                        evaluators->reinit(cell, positions, {old_solution_values.data(), old_solution_values.size()}, update_flags);

                      // To store the old solutions
                      Vector<double> old_solution;
                      if (update_flags & update_values)
                        old_solution.reinit(this->introspection().n_components);

                      // To store the old gradients
                      std::vector<Tensor<1,dim>> old_gradients;
                      if (update_flags & update_gradients)
                        old_gradients.resize(this->introspection().n_components);

                      // Loop over all particles in the cell
                      auto particle = particles_in_cell.begin();
                      for (unsigned int i = 0; particle!=particles_in_cell.end(); ++particle,++i)
                        {
                          // Evaluate the old solution, but only if it is requested in the update_flags
                          if (update_flags & update_values)
                            evaluators->get_solution(i, old_solution);

                          // Evaluate the old gradients, but only if they are requested in the update_flags
                          if (update_flags & update_gradients)
                            evaluators->get_gradients(i, old_gradients);

                          // Fill material model input
                          // TODO this is very similar to the slow path -> refactor?

                          // Get the real location of the particle
                          material_inputs.position[0] = particle->get_location();

                          material_inputs.current_cell = typename DoFHandler<dim>::active_cell_iterator(*particle->get_surrounding_cell(),
                                                                                                        &(this->get_dof_handler()));

                          material_inputs.temperature[0] = old_solution[this->introspection().component_indices.temperature];

                          material_inputs.pressure[0] = old_solution[this->introspection().component_indices.pressure];

                          for (unsigned int d = 0; d < dim; ++d)
                            material_inputs.velocity[0][d] = old_solution[this->introspection().component_indices.velocities[d]];

                          // For the ve_stress_* fields, we use the values on the particles. After they have been restored,
                          // their properties have the values of the previous timestep.
                          for (unsigned int n = 0; n < this->n_compositional_fields(); ++n)
                            material_inputs.composition[0][n] = old_solution[this->introspection().component_indices.compositional_fields[n]];
                          for (unsigned int n = 0; n < 2*SymmetricTensor<2,dim>::n_independent_components; ++n)
                            material_inputs.composition[0][n] = particle->get_properties()[data_position + n];

                          Tensor<2,dim> grad_u;
                          for (unsigned int d=0; d<dim; ++d)
                            grad_u[d] = old_gradients[d];
                          material_inputs.strain_rate[0] = symmetrize (grad_u);

                          // Evaluate the material model to get the reaction rates
                          this->get_material_model().evaluate (material_inputs,material_outputs);

                          // Add the reaction_rates * timestep = update to the corresponding stress
                          // tensor components
                          for (unsigned int p = 0; p < 2*SymmetricTensor<2,dim>::n_independent_components ; ++p)
                            particle->get_properties()[data_position + p] += reaction_rate_outputs->reaction_rates[0][p] * this->get_timestep();
                        }
                    }
                  // Slow path
                  else
                    {
                      const unsigned int solution_components = this->introspection().n_components;

                      // The old solutions and old gradients of one particle
                      Vector<double>             old_value (solution_components);
                      std::vector<Tensor<1,dim>> old_gradient (solution_components,Tensor<1,dim>());

                      // Vector of old solutions and gradients for all particles in one cell
                      std::vector<Vector<double>>             old_solution(n_particles_in_cell,old_value);
                      std::vector<std::vector<Tensor<1,dim>>> old_gradients(n_particles_in_cell,old_gradient);
                      std::vector<Point<dim>>                 positions(n_particles_in_cell);

                      // Get the locations of the particles within the reference cell
                      auto it = particles_in_cell.begin();
                      for (unsigned int i = 0; it!=particles_in_cell.end(); ++it,++i)
                        {
                          positions[i] = it->get_reference_location();
                        }

                      // Set up a quadrature rule to evaluate the finite element
                      // solutions at the location of the particles
                      const Quadrature<dim> quadrature_formula(positions);
                      FEValues<dim> fe_value (this->get_mapping(),
                                              this->get_fe(),
                                              quadrature_formula,
                                              update_flags);

                      // Update to the current cell
                      fe_value.reinit (cell);

                      // Get the old solution values
                      if (update_flags & update_values)
                        fe_value.get_function_values (this->get_old_solution(),
                                                      old_solution);
                      // Get the old gradient values
                      if (update_flags & update_gradients)
                        fe_value.get_function_gradients (this->get_old_solution(),
                                                         old_gradients);

                      // Loop over all particles in the cell
                      auto particle = particles_in_cell.begin();
                      for (unsigned int i = 0; particle!=particles_in_cell.end(); ++particle,++i)
                        {
                          // Fill material model input
                          material_inputs.position[0] = particle->get_location();

                          material_inputs.current_cell = typename DoFHandler<dim>::active_cell_iterator(*particle->get_surrounding_cell(),
                                                                                                        &(this->get_dof_handler()));
                          material_inputs.temperature[0] = old_solution[i][this->introspection().component_indices.temperature];

                          material_inputs.pressure[0] = old_solution[i][this->introspection().component_indices.pressure];

                          for (unsigned int d = 0; d < dim; ++d)
                            material_inputs.velocity[0][d] = old_solution[i][this->introspection().component_indices.velocities[d]];

                          // For the ve_stress_* fields, we use the values on the particles
                          for (unsigned int n = 0; n < this->n_compositional_fields(); ++n)
                            material_inputs.composition[0][n] = old_solution[i][this->introspection().component_indices.compositional_fields[n]];
                          for (unsigned int n = 0; n < 2*SymmetricTensor<2,dim>::n_independent_components; ++n)
                            material_inputs.composition[0][n] = particle->get_properties()[data_position + n];

                          Tensor<2,dim> grad_u;
                          for (unsigned int d=0; d<dim; ++d)
                            grad_u[d] = old_gradients[i][d];
                          material_inputs.strain_rate[0] = symmetrize (grad_u);

                          // Evaluate the material model to get the reaction rates
                          this->get_material_model().evaluate (material_inputs,material_outputs);

                          // Add the reaction_rates * timestep = update to the corresponding stress
                          // tensor components
                          for (unsigned int p = 0; p < 2*SymmetricTensor<2,dim>::n_independent_components ; ++p)
                            particle->get_properties()[data_position + p] += reaction_rate_outputs->reaction_rates[0][p] * this->get_timestep();
                        }
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
      ElasticStress<dim>::update_particle_property(const unsigned int data_position,
                                                   const Vector<double> &solution,
                                                   const std::vector<Tensor<1,dim>> &gradients,
                                                   typename ParticleHandler<dim>::particle_iterator &particle) const
      {
        material_inputs.position[0] = particle->get_location();


        material_inputs.current_cell = typename DoFHandler<dim>::active_cell_iterator(*particle->get_surrounding_cell(),
                                                                                      &(this->get_dof_handler()));

        material_inputs.temperature[0] = solution[this->introspection().component_indices.temperature];

        material_inputs.pressure[0] = solution[this->introspection().component_indices.pressure];

        for (unsigned int d = 0; d < dim; ++d)
          material_inputs.velocity[0][d] = solution[this->introspection().component_indices.velocities[d]];

        // Instead of using the material model inputs, we use the ve_stress_* stored on particles. Other fields are copied from
        // the solution.
        for (unsigned int n = 0; n < this->n_compositional_fields(); ++n)
          material_inputs.composition[0][n] = solution[this->introspection().component_indices.compositional_fields[n]];
        for (unsigned int n = 0; n < 2*SymmetricTensor<2,dim>::n_independent_components; ++n)
          material_inputs.composition[0][n] = particle->get_properties()[data_position + n];

        Tensor<2,dim> grad_u;
        for (unsigned int d=0; d<dim; ++d)
          grad_u[d] = gradients[d];
        material_inputs.strain_rate[0] = symmetrize (grad_u);

        this->get_material_model().evaluate (material_inputs,material_outputs);

        // For the second set of stresses, the ve_stress_*_old fields, the update will be zero
        for (unsigned int i = 0; i < 2*SymmetricTensor<2,dim>::n_independent_components ; ++i)
          particle->get_properties()[data_position + i] += material_outputs.reaction_terms[0][i];
      }



      template <int dim>
      UpdateTimeFlags
      ElasticStress<dim>::need_update() const
      {
        return update_time_step;
      }



      template <int dim>
      UpdateFlags
      ElasticStress<dim>::get_needed_update_flags () const
      {
        return update_values | update_gradients;
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
