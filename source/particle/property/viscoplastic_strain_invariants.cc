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

#include <aspect/particle/property/viscoplastic_strain_invariants.h>
#include <aspect/material_model/visco_plastic.h>
#include <aspect/initial_composition/interface.h>


namespace aspect
{
  namespace Particle
  {
    namespace Property
    {

      template <int dim>
      ViscoPlasticStrainInvariant<dim>::ViscoPlasticStrainInvariant ()
        :
        n_components(0),
        material_inputs(1,0)
      {}



      template <int dim>
      void
      ViscoPlasticStrainInvariant<dim>::initialize ()
      {
        AssertThrow(Plugins::plugin_type_matches<const MaterialModel::ViscoPlastic<dim>>
                    (this->get_material_model()),
                    ExcMessage("This initial condition only makes sense in combination "
                               "with the visco_plastic material model."));

        n_components = 0;
        material_inputs = MaterialModel::MaterialModelInputs<dim>(1,this->n_compositional_fields());

        // Find out which fields are used.
        if (this->introspection().compositional_name_exists("plastic_strain"))
          n_components += 1;

        if (this->introspection().compositional_name_exists("viscous_strain"))
          n_components += 1;

        if (this->introspection().compositional_name_exists("total_strain"))
          n_components += 1;

        if (this->introspection().compositional_name_exists("noninitial_plastic_strain"))
          n_components += 1;

        if (n_components == 0)
          AssertThrow(false,
                      ExcMessage("This particle property requires a compositional "
                                 "strain field (plastic_strain, viscous_strain, "
                                 "or total_strain)."));
      }



      template <int dim>
      void
      ViscoPlasticStrainInvariant<dim>::initialize_one_particle_property(const Point<dim> &position,
                                                                         std::vector<double> &data) const
      {
        // Give each strain field its initial composition if one is prescribed.
        if (this->introspection().compositional_name_exists("plastic_strain"))
          data.push_back(this->get_initial_composition_manager().initial_composition(position,this->introspection().compositional_index_for_name("plastic_strain")));

        if (this->introspection().compositional_name_exists("viscous_strain"))
          data.push_back(this->get_initial_composition_manager().initial_composition(position,this->introspection().compositional_index_for_name("viscous_strain")));

        if (this->introspection().compositional_name_exists("total_strain"))
          data.push_back(this->get_initial_composition_manager().initial_composition(position,this->introspection().compositional_index_for_name("total_strain")));

        if (this->introspection().compositional_name_exists("noninitial_plastic_strain"))
          data.push_back(this->get_initial_composition_manager().initial_composition(position,this->introspection().compositional_index_for_name("noninitial_plastic_strain")));

      }



      template <int dim>
      void
      ViscoPlasticStrainInvariant<dim>::update_particle_properties(const ParticleUpdateInputs<dim> &inputs,
                                                                   typename ParticleHandler<dim>::particle_iterator_range &particles) const
      {
        // Find out plastic yielding by calling function in material model.
        const MaterialModel::ViscoPlastic<dim> &viscoplastic
          = Plugins::get_plugin_as_type<const MaterialModel::ViscoPlastic<dim>>(this->get_material_model());

        // Current timestep
        const double dt = this->get_timestep();
        const unsigned int data_position = this->data_position;

        unsigned int p = 0;
        for (auto &particle: particles)
          {

            // Velocity gradients
            Tensor<2,dim> grad_u;
            for (unsigned int d=0; d<dim; ++d)
              grad_u[d] = inputs.gradients[p][d];

            material_inputs.pressure[0] = inputs.solution[p][this->introspection().component_indices.pressure];
            material_inputs.temperature[0] = inputs.solution[p][this->introspection().component_indices.temperature];
            material_inputs.position[0] = particle.get_location();

            // Calculate strain rate from velocity gradients
            material_inputs.strain_rate[0] = symmetrize (grad_u);

            // Put compositional fields into single variable
            for (unsigned int i = 0; i < this->n_compositional_fields(); ++i)
              {
                material_inputs.composition[0][i] = inputs.solution[p][this->introspection().component_indices.compositional_fields[i]];
              }


            const bool plastic_yielding = viscoplastic.is_yielding(material_inputs);

            // Next take the integrated strain invariant from the prior time step.
            const auto data = particle.get_properties();

            // Calculate strain rate second invariant
            const double edot_ii = std::sqrt(std::max(-second_invariant(deviator(material_inputs.strain_rate[0])), 0.));

            // Calculate strain invariant magnitude over the last time step
            const double strain_update = dt*edot_ii;

            /* Update the strain values that are used in the simulation, which use the following assumptions
             * to identify the correct position in the data vector for each value:
             * (1) Total strain cannot be used in combination with any other strain field
             * (2) If plastic strain is tracked, it will always be in the first data position
             * (3) If noninitial plastic strain is tracked, it will always be in the last data position
             * (4) If noninitial plastic strain is tracked, plastic strain is also being tracked
             * (5) If only viscous strain is tracked, it will be in the first data position.
             * (6) If both viscous and plastic strain are tracked, viscous strain will be in the second data position
             * If these assumptions change in the future, they will need to be updated.
             * */

            if (this->introspection().compositional_name_exists("plastic_strain") && plastic_yielding == true)
              data[data_position] += strain_update;

            if (this->introspection().compositional_name_exists("viscous_strain") && plastic_yielding == false)
              {
                // Not yielding and only one field, which tracks the viscous strain.
                if (n_components == 1)
                  data[data_position] += strain_update;

                // Not yielding and either two or three fields are tracked. If two fields are tracked,
                // they represent plastic strain (first data position) and viscous strain (second data
                // data position, updated below). If three fields are tracked, they represent plastic
                // strain (first data position), viscous strain (second data position, updated below),
                // and noninitial plastic strain (third data position). In either case, the viscous
                // strain is in the second data position, allowing us to use a single expression.
                if (n_components > 1)
                  data[data_position+1] += strain_update;
              }

            // Only one field, which tracks total strain and is updated regardless of whether the
            // material is yielding or not.
            if (this->introspection().compositional_name_exists("total_strain"))
              data[data_position] += strain_update;

            // Yielding, and noninitial plastic strain (last data position, updated below) is tracked.
            if (this->introspection().compositional_name_exists("noninitial_plastic_strain") && plastic_yielding == true)
              data[data_position+(n_components-1)] += strain_update;

            ++p;
          }
      }



      template <int dim>
      UpdateTimeFlags
      ViscoPlasticStrainInvariant<dim>::need_update() const
      {
        return update_time_step;
      }



      template <int dim>
      UpdateFlags
      ViscoPlasticStrainInvariant<dim>::get_update_flags (const unsigned int component) const
      {
        if (this->introspection().component_masks.velocities[component] == true)
          return update_values | update_gradients;

        return update_values;
      }

      template <int dim>
      std::vector<std::pair<std::string, unsigned int>>
      ViscoPlasticStrainInvariant<dim>::get_property_information() const
      {
        std::vector<std::pair<std::string,unsigned int>> property_information;

        //Check which fields are used in model and make an output for each.
        if (this->introspection().compositional_name_exists("plastic_strain"))
          property_information.emplace_back("plastic_strain", 1);

        if (this->introspection().compositional_name_exists("viscous_strain"))
          property_information.emplace_back("viscous_strain", 1);

        if (this->introspection().compositional_name_exists("total_strain"))
          property_information.emplace_back("total_strain", 1);

        if (this->introspection().compositional_name_exists("noninitial_plastic_strain"))
          property_information.emplace_back("noninitial_plastic_strain", 1);

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
      ASPECT_REGISTER_PARTICLE_PROPERTY(ViscoPlasticStrainInvariant,
                                        "viscoplastic strain invariants",
                                        "A plugin that calculates the finite strain invariant a particle has "
                                        "experienced and assigns it to either the plastic and/or viscous strain field based "
                                        "on whether the material is plastically yielding, or the total strain field "
                                        "used in the visco plastic material model. The implementation of this property "
                                        "is equivalent to the implementation for compositional fields that is located in "
                                        "the plugin in \\texttt{benchmarks/buiter\\_et\\_al\\_2008\\_jgr/plugin/},"
                                        "and is effectively the same as what the visco plastic material model uses for compositional fields.")
    }
  }
}
