/*
  Copyright (C) 2015 - 2021 by the authors of the ASPECT code.

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
        material_inputs(1,0),
        material_outputs(1,0)
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
        material_outputs = MaterialModel::MaterialModelOutputs<dim>(1, this->n_compositional_fields());


        // Find out which fields are used.
        if (this->introspection().compositional_name_exists("plastic_strain"))
          n_components += 1;

        if (this->introspection().compositional_name_exists("viscous_strain"))
          n_components += 1;

        if (this->introspection().compositional_name_exists("total_strain"))
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

      }



      template <int dim>
      void
      ViscoPlasticStrainInvariant<dim>::update_particle_property(const unsigned int data_position,
                                                                 const Vector<double> &solution,
                                                                 const std::vector<Tensor<1,dim>> &gradients,
                                                                 typename ParticleHandler<dim>::particle_iterator &particle) const
      {
        // Velocity gradients
        Tensor<2,dim> grad_u;
        for (unsigned int d=0; d<dim; ++d)
          grad_u[d] = gradients[d];

        material_inputs.pressure[0] = solution[this->introspection().component_indices.pressure];
        material_inputs.temperature[0] = solution[this->introspection().component_indices.temperature];
        material_inputs.position[0] = particle->get_location();

        // Calculate strain rate from velocity gradients
        material_inputs.strain_rate[0] = symmetrize (grad_u);

        // Put compositional fields into single variable
        for (unsigned int i = 0; i < this->n_compositional_fields(); i++)
          {
            material_inputs.composition[0][i] = solution[this->introspection().component_indices.compositional_fields[i]];
          }

        // Evaluate directly in the viscoplastic material model and modify the reaction outputs
        this->get_material_model().evaluate (material_inputs,material_outputs);

        const int plastic_strain_index = this->introspection().compositional_index_for_name("plastic_strain"); 
        const int viscous_strain_index = this->introspection().compositional_index_for_name("viscous_strain");
        const int total_strain_index = this->introspection().compositional_index_for_name("total_strain");
        
        if (this->introspection().compositional_name_exists("plastic_strain"))
          particle->get_properties()[data_position] += material_outputs.reaction_terms[0][plastic_strain_index];
        
        if (this->introspection().compositional_name_exists("viscous_strain"))
          particle->get_properties()[data_position] += material_outputs.reaction_terms[0][viscous_strain_index];

        if (this->introspection().compositional_name_exists("total_strain"))
          particle->get_properties()[data_position] += material_outputs.reaction_terms[0][total_strain_index];  

        //for (unsigned int i = 0; i < n_components-1 ; ++i)
        //  particle->get_properties()[data_position + i] += material_outputs.reaction_terms[0][i];
      }



      template <int dim>
      UpdateTimeFlags
      ViscoPlasticStrainInvariant<dim>::need_update() const
      {
        return update_time_step;
      }



      template <int dim>
      UpdateFlags
      ViscoPlasticStrainInvariant<dim>::get_needed_update_flags () const
      {
        // Need to update both of these to send into material model.
        return update_values | update_gradients;
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
