/*
  Copyright (C) 2015 - 2017 by the authors of the ASPECT code.

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
      {}

      template <int dim>
      void
      ViscoPlasticStrainInvariant<dim>::initialize ()
      {
        AssertThrow(dynamic_cast<const MaterialModel::ViscoPlastic<dim> *>(&this->get_material_model()) != NULL,
                    ExcMessage("This initial condition only makes sense in combination with the visco_plastic material model."));

        n_components = 0;

        // Find out which fields are used.
        if (this->introspection().compositional_name_exists("plastic_strain"))
          n_components += 1;

        if (this->introspection().compositional_name_exists("viscous_strain"))
          n_components += 1;
        else if (this->introspection().compositional_name_exists("total_strain") && !this->introspection().compositional_name_exists("plastic_strain"))
          n_components = 1;

        if (n_components == 0)
          AssertThrow(false, ExcMessage("This particle property requires a compositional strain field (plastic_strain, viscous_strain, or total_strain)."));
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
        else if (this->introspection().compositional_name_exists("total_strain") && !this->introspection().compositional_name_exists("plastic_strain"))
          data.push_back(this->get_initial_composition_manager().initial_composition(position,this->introspection().compositional_index_for_name("total_strain")));

      }

      template <int dim>
      void
      ViscoPlasticStrainInvariant<dim>::update_one_particle_property(const unsigned int data_position,
                                                                     const Point<dim> &,
                                                                     const Vector<double> &solution,
                                                                     const std::vector<Tensor<1,dim> > &gradients,
                                                                     const ArrayView<double> &data) const
      {
        // Current timestep
        const double dt = this->get_timestep();

        // Velocity gradients
        Tensor<2,dim> grad_u;
        for (unsigned int d=0; d<dim; ++d)
          grad_u[d] = gradients[d];

        // Calculate strain rate from velocity gradients
        const SymmetricTensor<2,dim> strain_rate = symmetrize (grad_u);


        //Put compositional fields into single variable
        std::vector<double> composition(this->n_compositional_fields());
        for (unsigned int i = 0; i < this->n_compositional_fields(); i++)
          {
            composition[i] = solution[this->introspection().component_indices.compositional_fields[i]];
          }

        // Find out plastic yielding by calling function in material model.
        const MaterialModel::ViscoPlastic<dim> *viscoplastic
          = dynamic_cast<const MaterialModel::ViscoPlastic<dim> *>(&this->get_material_model());

        bool plastic_yielding = false;
        plastic_yielding = viscoplastic->get_plastic_yielding(solution[this->introspection().component_indices.pressure],
                                                              solution[this->introspection().component_indices.temperature],
                                                              composition,
                                                              strain_rate);


        /* Integrated strain invariant from prior time step, grab second data component
         * for viscous strain if both viscous and plastic are used.
         */
        double old_strain = data[data_position];
        if (n_components == 2 && plastic_yielding == false)
          old_strain = data[data_position+(n_components-1)];


        // Calculate strain rate second invariant
        const double edot_ii = std::sqrt(std::fabs(second_invariant(deviator(strain_rate))));

        // New strain is the old strain plus dt*edot_ii
        const double new_strain = old_strain + dt*edot_ii;


        // Put strain into correct data position.
        if (this->introspection().compositional_name_exists("plastic_strain") && plastic_yielding == true)
          data[data_position] = new_strain;

        if (this->introspection().compositional_name_exists("viscous_strain") && plastic_yielding == false)
          data[data_position+(n_components-1)] = new_strain;
        else if (this->introspection().compositional_name_exists("total_strain") && !this->introspection().compositional_name_exists("plastic_strain"))
          data[data_position] = new_strain;

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
      std::vector<std::pair<std::string, unsigned int> >
      ViscoPlasticStrainInvariant<dim>::get_property_information() const
      {

        std::vector<std::pair<std::string,unsigned int> > property_information;

        //Check which fields are used in model and make an output for each.
        if (this->introspection().compositional_name_exists("plastic_strain"))
          property_information.emplace_back("plastic_strain",1);

        if (this->introspection().compositional_name_exists("viscous_strain"))
          property_information.emplace_back("viscous_strain",1);
        else if (this->introspection().compositional_name_exists("total_strain") && !this->introspection().compositional_name_exists("plastic_strain"))
          property_information.emplace_back("total_strain",1);

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
                                        "A plugin which calculates the the integrated strain invariant "
                                        "(as described in the integrated_strain_invariant particle property) "
                                        "and assigns the strain to the plastic, viscous, or total strain fields "
                                        "used in the visco_plastic material model.")
    }
  }
}
