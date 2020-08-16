/*
  Copyright (C) 2015 - 2020 by the authors of the ASPECT code.

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
        material_inputs = MaterialModel::MaterialModelInputs<dim>(1, this->n_compositional_fields());

        material_outputs = MaterialModel::MaterialModelOutputs<dim>(1, this->n_compositional_fields());

        AssertThrow((Plugins::plugin_type_matches<const MaterialModel::ViscoPlastic<dim>>(this->get_material_model())
                     ||
                     Plugins::plugin_type_matches<const MaterialModel::Viscoelastic<dim>>(this->get_material_model())),
                    ExcMessage("This particle property only makes sense in combination with the viscoelastic or visco_plastic material model."));

        AssertThrow(this->get_parameters().enable_elasticity == true,
                    ExcMessage ("This particle property should only be used if 'Enable elasticity' is set to true"));

      }



      template <int dim>
      void
      ElasticStress<dim>::initialize_one_particle_property(const Point<dim> &position,
                                                           std::vector<double> &data) const
      {
        // Give each elastic stress field its initial composition if one is prescribed.
        data.push_back(this->get_initial_composition_manager().initial_composition(position,this->introspection().compositional_index_for_name("stress_xx")));

        data.push_back(this->get_initial_composition_manager().initial_composition(position,this->introspection().compositional_index_for_name("stress_yy")));

        if (dim == 2)
          {
            data.push_back(this->get_initial_composition_manager().initial_composition(position,this->introspection().compositional_index_for_name("stress_xy")));
          }
        else if (dim == 3)
          {
            data.push_back(this->get_initial_composition_manager().initial_composition(position,this->introspection().compositional_index_for_name("stress_zz")));

            data.push_back(this->get_initial_composition_manager().initial_composition(position,this->introspection().compositional_index_for_name("stress_xy")));

            data.push_back(this->get_initial_composition_manager().initial_composition(position,this->introspection().compositional_index_for_name("stress_xz")));

            data.push_back(this->get_initial_composition_manager().initial_composition(position,this->introspection().compositional_index_for_name("stress_yz")));
          }

      }



      template <int dim>
      void
      ElasticStress<dim>::update_particle_property(const unsigned int data_position,
                                                   const Vector<double> &solution,
                                                   const std::vector<Tensor<1,dim> > &gradients,
                                                   typename ParticleHandler<dim>::particle_iterator &particle) const
      {
        material_inputs.position[0] = particle->get_location();

        material_inputs.current_cell = typename DoFHandler<dim>::active_cell_iterator(*particle->get_surrounding_cell(this->get_triangulation()),
                                                                                      &(this->get_dof_handler()));

        material_inputs.temperature[0] = solution[this->introspection().component_indices.temperature];

        material_inputs.pressure[0] = solution[this->introspection().component_indices.pressure];

        for (unsigned int d = 0; d < dim; ++d)
          material_inputs.velocity[0][d] = solution[this->introspection().component_indices.velocities[d]];

        for (unsigned int n = 0; n < this->n_compositional_fields(); ++n)
          material_inputs.composition[0][n] = solution[this->introspection().component_indices.compositional_fields[n]];

        Tensor<2,dim> grad_u;
        for (unsigned int d=0; d<dim; ++d)
          grad_u[d] = gradients[d];
        material_inputs.strain_rate[0] = symmetrize (grad_u);

        this->get_material_model().evaluate (material_inputs,material_outputs);

        for (unsigned int i = 0; i < SymmetricTensor<2,dim>::n_independent_components ; ++i)
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
      std::vector<std::pair<std::string, unsigned int> >
      ElasticStress<dim>::get_property_information() const
      {
        std::vector<std::pair<std::string,unsigned int> > property_information;

        //Check which fields are used in model and make an output for each.
        if (this->introspection().compositional_name_exists("stress_xx"))
          property_information.emplace_back("stress_xx",1);

        if (this->introspection().compositional_name_exists("stress_yy"))
          property_information.emplace_back("stress_yy",1);

        if (dim == 2)
          {
            if (this->introspection().compositional_name_exists("stress_xy"))
              property_information.emplace_back("stress_xy",1);
          }
        else if (dim == 3)
          {
            if (this->introspection().compositional_name_exists("stress_zz"))
              property_information.emplace_back("stress_zz",1);

            if (this->introspection().compositional_name_exists("stress_xy"))
              property_information.emplace_back("stress_xy",1);

            if (this->introspection().compositional_name_exists("stress_xz"))
              property_information.emplace_back("stress_xz",1);

            if (this->introspection().compositional_name_exists("stress_yz"))
              property_information.emplace_back("stress_yz",1);
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
                                        "accumulated. See the viscoelastic material model "
                                        "documentation for more detailed information.")

    }
  }
}

