/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

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


#include <aspect/global.h>
#include <aspect/utilities.h>
#include <aspect/heating_model/interface.h>
#include <aspect/heating_model/adiabatic_heating.h>
#include <aspect/heating_model/shear_heating.h>

#include <deal.II/base/signaling_nan.h>

#include <tuple>
#include <list>


namespace aspect
{
  namespace HeatingModel
  {
    template <int dim>
    void
    Interface<dim>::
    create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> & /*outputs*/) const
    {}



    template <int dim>
    void
    Interface<dim>::
    create_additional_material_model_inputs(MaterialModel::MaterialModelInputs<dim> & /*inputs*/) const
    {}


    template <int dim>
    MaterialModel::MaterialProperties::Property
    Interface<dim>::
    get_required_properties() const
    {
      return MaterialModel::MaterialProperties::all_properties;
    }



    // ------------------------------ Manager -----------------------------



    template <int dim>
    bool
    Manager<dim>::adiabatic_heating_enabled() const
    {
      return this->template has_matching_active_plugin<HeatingModel::AdiabaticHeating<dim>>() ;
    }



    template <int dim>
    bool
    Manager<dim>::shear_heating_enabled() const
    {
      return this->template has_matching_active_plugin<HeatingModel::ShearHeating<dim>>() ;
    }



    namespace
    {
      std::tuple
      <aspect::internal::Plugins::UnusablePluginList,
      aspect::internal::Plugins::UnusablePluginList,
      aspect::internal::Plugins::PluginList<Interface<2>>,
      aspect::internal::Plugins::PluginList<Interface<3>>> registered_plugins;
    }


    template <int dim>
    void
    Manager<dim>::register_heating_model (const std::string &name,
                                          const std::string &description,
                                          void (*declare_parameters_function) (ParameterHandler &),
                                          std::unique_ptr<Interface<dim>> (*factory_function) ())
    {
      std::get<dim>(registered_plugins).register_plugin (name,
                                                         description,
                                                         declare_parameters_function,
                                                         factory_function);
    }


    template <int dim>
    void
    Manager<dim>::parse_parameters (ParameterHandler &prm)
    {
      // find out which plugins are requested and the various other
      // parameters we declare here
      prm.enter_subsection ("Heating model");
      {
        this->plugin_names
          = Utilities::split_string_list(prm.get("List of model names"));

        AssertThrow(Utilities::has_unique_entries(this->plugin_names),
                    ExcMessage("The list of strings for the parameter "
                               "'Heating model/List of model names' contains entries more than once. "
                               "This is not allowed. Please check your parameter file."));
      }
      prm.leave_subsection ();

      // go through the list, create objects and let them parse
      // their own parameters
      for (auto &model_name : this->plugin_names)
        {
          this->plugin_objects.emplace_back (std::get<dim>(registered_plugins)
                                             .create_plugin (model_name,
                                                             "Heating model::Model names"));

          if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(&*this->plugin_objects.back()))
            sim->initialize_simulator (this->get_simulator());

          this->plugin_objects.back()->parse_parameters (prm);
          this->plugin_objects.back()->initialize ();
        }
    }


    template <int dim>
    void
    Manager<dim>::evaluate (const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
                            const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
                            HeatingModel::HeatingModelOutputs &heating_model_outputs) const
    {
      // Initialize all outputs to zero, because heating_model_outputs is
      // often reused in loops over all cells
      for (unsigned int q=0; q<heating_model_outputs.heating_source_terms.size(); ++q)
        {
          heating_model_outputs.heating_source_terms[q] = 0.0;
          heating_model_outputs.lhs_latent_heat_terms[q] = 0.0;
          heating_model_outputs.rates_of_temperature_change[q] = 0.0;
        }

      HeatingModel::HeatingModelOutputs individual_heating_outputs(material_model_inputs.n_evaluation_points(),
                                                                   this->n_compositional_fields());

      for (const auto &heating_model : this->plugin_objects)
        {
          heating_model->evaluate(material_model_inputs, material_model_outputs, individual_heating_outputs);
          for (unsigned int q=0; q<heating_model_outputs.heating_source_terms.size(); ++q)
            {
              heating_model_outputs.heating_source_terms[q] += individual_heating_outputs.heating_source_terms[q];
              heating_model_outputs.lhs_latent_heat_terms[q] += individual_heating_outputs.lhs_latent_heat_terms[q];

              if (!this->get_parameters().use_operator_splitting)
                Assert(individual_heating_outputs.rates_of_temperature_change[q] == 0.0,
                       ExcMessage("Rates of temperature change heating model outputs have to be zero "
                                  "if the model does not use operator splitting."));
              heating_model_outputs.rates_of_temperature_change[q] += individual_heating_outputs.rates_of_temperature_change[q];
            }
          individual_heating_outputs.reset();
        }

      // If the heating model does not get the reaction rate outputs, it can not correctly compute
      // the rates of temperature change. To make sure these (incorrect) values are never used anywhere,
      // overwrite them with signaling_NaNs.
      const std::shared_ptr<const MaterialModel::ReactionRateOutputs<dim>> reaction_rate_outputs
        = material_model_outputs.template get_additional_output_object<MaterialModel::ReactionRateOutputs<dim>>();

      if (reaction_rate_outputs == nullptr)
        for (double &q : heating_model_outputs.rates_of_temperature_change)
          q = numbers::signaling_nan<double>();
    }



    template <int dim>
    void
    Manager<dim>::
    create_additional_material_model_inputs_and_outputs(MaterialModel::MaterialModelInputs<dim>  &material_model_inputs,
                                                        MaterialModel::MaterialModelOutputs<dim> &material_model_outputs) const
    {
      for (const auto &heating_model : this->plugin_objects)
        {
          heating_model->create_additional_material_model_inputs(material_model_inputs);
        }

      for (const auto &heating_model : this->plugin_objects)
        {
          if ((heating_model->get_required_properties() & MaterialModel::MaterialProperties::additional_outputs) != 0)
            heating_model->create_additional_material_model_outputs(material_model_outputs);
        }
    }



    template <int dim>
    const std::vector<std::string> &
    Manager<dim>::get_active_heating_model_names () const
    {
      return this->plugin_names;
    }


    template <int dim>
    const std::list<std::unique_ptr<Interface<dim>>> &
    Manager<dim>::get_active_heating_models () const
    {
      return this->plugin_objects;
    }


    template <int dim>
    void
    Manager<dim>::declare_parameters (ParameterHandler &prm)
    {
      // declare the actual entry in the parameter file
      prm.enter_subsection ("Heating model");
      {
        const std::string pattern_of_names
          = std::get<dim>(registered_plugins).get_pattern_of_names ();

        prm.declare_entry("List of model names",
                          "",
                          Patterns::MultipleSelection(pattern_of_names),
                          "A comma separated list of heating models that "
                          "will be used to calculate the heating terms in the energy "
                          "equation. The results of each of these criteria, i.e., "
                          "the heating source terms and the latent heat terms for the "
                          "left hand side will be added.\n\n"
                          "The following heating models are available:\n\n"
                          +
                          std::get<dim>(registered_plugins).get_description_string());
      }
      prm.leave_subsection ();

      std::get<dim>(registered_plugins).declare_parameters (prm);
    }



    template <int dim>
    void
    Manager<dim>::write_plugin_graph (std::ostream &out)
    {
      std::get<dim>(registered_plugins).write_plugin_graph ("Heating model interface",
                                                            out);
    }



    HeatingModelOutputs::HeatingModelOutputs(const unsigned int n_points,
                                             const unsigned int)
      :
      heating_source_terms(n_points,numbers::signaling_nan<double>()),
      // initialize the reaction terms with zeroes because they are not filled
      // in all heating models
      rates_of_temperature_change(n_points,0.0),
      lhs_latent_heat_terms(n_points,numbers::signaling_nan<double>())
    {
    }



    void
    HeatingModelOutputs::reset()
    {
      for (unsigned int q=0; q<heating_source_terms.size(); ++q)
        {
          heating_source_terms[q] = numbers::signaling_nan<double>();
          lhs_latent_heat_terms[q] = numbers::signaling_nan<double>();
          rates_of_temperature_change[q] = 0.0;
        }
    }



    template <int dim>
    std::string
    get_valid_model_names_pattern ()
    {
      return std::get<dim>(registered_plugins).get_pattern_of_names ();
    }

  }
}

// explicit instantiations
namespace aspect
{
  namespace internal
  {
    namespace Plugins
    {
      template <>
      std::list<internal::Plugins::PluginList<HeatingModel::Interface<2>>::PluginInfo> *
      internal::Plugins::PluginList<HeatingModel::Interface<2>>::plugins = nullptr;
      template <>
      std::list<internal::Plugins::PluginList<HeatingModel::Interface<3>>::PluginInfo> *
      internal::Plugins::PluginList<HeatingModel::Interface<3>>::plugins = nullptr;
    }
  }

  namespace HeatingModel
  {
#define INSTANTIATE(dim) \
  template class Interface<dim>; \
  template class Manager<dim>; \
  \
  template \
  std::string \
  get_valid_model_names_pattern<dim> ();

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
