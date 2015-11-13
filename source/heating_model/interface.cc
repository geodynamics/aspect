/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include <aspect/global.h>
#include <aspect/heating_model/interface.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/std_cxx11/tuple.h>

#include <list>


namespace aspect
{
  namespace HeatingModel
  {
    template <int dim>
    Interface<dim>::~Interface ()
    {}


    template <int dim>
    void
    Interface<dim>::initialize ()
    {}



    template <int dim>
    void
    Interface<dim>::update ()
    {}



#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
    template <int dim>
    void
    Interface<dim>::evaluate (const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
                              const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
                              HeatingModel::HeatingModelOutputs &heating_model_outputs) const
    {
      Assert(heating_model_outputs.heating_source_terms.size() == material_model_inputs.position.size(),
             ExcMessage ("Heating outputs need to have the same number of entries as the material model inputs."));
      for (unsigned int q=0; q<heating_model_outputs.heating_source_terms.size(); ++q)
        {
          heating_model_outputs.heating_source_terms[q] = specific_heating_rate(material_model_inputs.temperature[q],
                                                                                material_model_inputs.pressure[q],
                                                                                material_model_inputs.composition[q],
                                                                                material_model_inputs.position[q])
                                                          * material_model_outputs.densities[q];
          heating_model_outputs.lhs_latent_heat_terms[q] = 0.0;
        }
    }
#pragma GCC diagnostic pop


    template <int dim>
    double
    Interface<dim>::specific_heating_rate (const double,
                                           const double,
                                           const std::vector<double> &,
                                           const Point<dim> &) const
    {
      Assert(false,
             ExcMessage ("There is no 'evaluate()' or 'specific_heating_rate()' function implemented in the heating model!"));
      return 0.0;
    }


    template <int dim>
    void
    Interface<dim>::
    declare_parameters (dealii::ParameterHandler &)
    {}


    template <int dim>
    void
    Interface<dim>::parse_parameters (dealii::ParameterHandler &)
    {}


    // ------------------------------ Manager -----------------------------

    template <int dim>
    Manager<dim>::~Manager()
    {}


    namespace
    {
      std_cxx11::tuple
      <void *,
      void *,
      internal::Plugins::PluginList<Interface<2> >,
      internal::Plugins::PluginList<Interface<3> > > registered_plugins;
    }


    template <int dim>
    void
    Manager<dim>::register_heating_model (const std::string &name,
                                          const std::string &description,
                                          void (*declare_parameters_function) (ParameterHandler &),
                                          Interface<dim> *(*factory_function) ())
    {
      std_cxx11::get<dim>(registered_plugins).register_plugin (name,
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
        model_names
          = Utilities::split_string_list(prm.get("List of model names"));

        const std::string model_name = prm.get ("Model name");

        AssertThrow (model_name == "unspecified" || model_names.size() == 0,
                     ExcMessage ("The parameter 'Model name' is only used for reasons"
                                 "of backwards compatibility and can not be used together with "
                                 "the new functionality 'List of model names'. Please add your "
                                 "heating model to the list instead."));

        if (!(model_name == "unspecified"))
          model_names.push_back(model_name);

      }
      prm.leave_subsection ();

      prm.enter_subsection ("Model settings");
      {
        const bool include_shear_heating = prm.get_bool ("Include shear heating");
        Assert(!(include_shear_heating && std::find(model_names.begin(), model_names.end(), "shear heating") != model_names.end()),
               ExcMessage ("Deprecated: The old functionality 'Include shear heating'"
                           "is only allowed for reasons of backwards compatibility and "
                           "can not be used together with the new functionality 'List "
                           "of model names'. Please remove the 'Include shear heating'"
                           "setting."));
        if (include_shear_heating)
          model_names.push_back("shear heating");

        const bool include_adiabatic_heating = prm.get_bool ("Include adiabatic heating");
        Assert(!(include_adiabatic_heating && std::find(model_names.begin(), model_names.end(), "adiabatic heating") != model_names.end()),
               ExcMessage ("Deprecated: The old functionality 'Include adiabatic heating'"
                           "is only allowed for reasons of backwards compatibility and "
                           "can not be used together with the new functionality 'List "
                           "of model names'. Please remove the 'Include adiabatic heating'"
                           "setting."));
        if (include_adiabatic_heating)
          model_names.push_back("adiabatic heating");

        const bool include_latent_heat = prm.get_bool ("Include latent heat");
        Assert(!(include_latent_heat && std::find(model_names.begin(), model_names.end(), "latent heat") != model_names.end()),
               ExcMessage ("Deprecated: The old functionality 'Include latent heat'"
                           "is only allowed for reasons of backwards compatibility and "
                           "can not be used together with the new functionality 'List "
                           "of model names'. Please remove the 'Include latent heat'"
                           "setting."));
        if (include_latent_heat)
          model_names.push_back("latent heat");
      }
      prm.leave_subsection ();



      // go through the list, create objects and let them parse
      // their own parameters
      for (unsigned int name=0; name<model_names.size(); ++name)
        {
          heating_model_objects.push_back (std_cxx11::shared_ptr<Interface<dim> >
                                           (std_cxx11::get<dim>(registered_plugins)
                                            .create_plugin (model_names[name],
                                                            "Heating model::Model names")));

          if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(&*heating_model_objects.back()))
            sim->initialize_simulator (this->get_simulator());

          heating_model_objects.back()->parse_parameters (prm);
          heating_model_objects.back()->initialize ();
        }
    }


    template <int dim>
    void
    Manager<dim>::update ()
    {
      for (typename std::list<std_cxx11::shared_ptr<HeatingModel::Interface<dim> > >::const_iterator
           heating_model = heating_model_objects.begin();
           heating_model != heating_model_objects.end(); ++heating_model)
        {
          (*heating_model)->update();
        }
    }

    template <int dim>
    void
    Manager<dim>::evaluate (const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
                            const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
                            HeatingModel::HeatingModelOutputs &heating_model_outputs) const
    {
      // the heating outputs are initialized with zeros, so there is no heating if they are not set
      // in the individual plugins
      HeatingModel::HeatingModelOutputs individual_heating_outputs(material_model_inputs.position.size(),
                                                                   this->n_compositional_fields());

      for (typename std::list<std_cxx11::shared_ptr<HeatingModel::Interface<dim> > >::const_iterator
           heating_model = heating_model_objects.begin();
           heating_model != heating_model_objects.end(); ++heating_model)
        {
          (*heating_model)->evaluate(material_model_inputs, material_model_outputs, individual_heating_outputs);
          for (unsigned int q=0; q<heating_model_outputs.heating_source_terms.size(); ++q)
            {
              heating_model_outputs.heating_source_terms[q] += individual_heating_outputs.heating_source_terms[q];
              heating_model_outputs.lhs_latent_heat_terms[q] += individual_heating_outputs.lhs_latent_heat_terms[q];
            }
        }
    }


    template <int dim>
    const std::vector<std::string> &
    Manager<dim>::get_active_heating_model_names () const
    {
      return model_names;
    }


    template <int dim>
    const std::list<std_cxx11::shared_ptr<Interface<dim> > > &
    Manager<dim>::get_active_heating_models () const
    {
      return heating_model_objects;
    }


    template <int dim>
    void
    Manager<dim>::declare_parameters (ParameterHandler &prm)
    {
      // declare the actual entry in the parameter file
      prm.enter_subsection ("Heating model");
      {
        const std::string pattern_of_names
          = std_cxx11::get<dim>(registered_plugins).get_pattern_of_names ();

        prm.declare_entry("List of model names",
                          "",
                          Patterns::MultipleSelection(pattern_of_names),
                          "A comma separated list of heating models that "
                          "will be used to calculate the heating terms in the energy"
                          "equation. The results of each of these criteria , i.e., "
                          "the heating source terms and the latent heat terms for the"
                          "left hand side will be added.\n\n"
                          "The following heating models are available:\n\n"
                          +
                          std_cxx11::get<dim>(registered_plugins).get_description_string());

        prm.declare_entry ("Model name", "unspecified",
                           Patterns::Selection (pattern_of_names+"|unspecified"),
                           "Select one of the following models:\n\n"
                           "Warning: This is the old formulation of specifying "
                           "heating models and shouldn't be used. Please use 'List of"
                           "model names' instead."
                           +
                           std_cxx11::get<dim>(registered_plugins).get_description_string());
      }
      prm.leave_subsection ();

      prm.enter_subsection ("Model settings");
      {
        prm.declare_entry ("Include shear heating", "false",
                           Patterns::Bool (),
                           "Whether to include shear heating into the model or not. From a "
                           "physical viewpoint, shear heating should always be used but may "
                           "be undesirable when comparing results with known benchmarks that "
                           "do not include this term in the temperature equation."
                           "Warning: deprecated! Add 'shear heating' to the 'List of model "
                           "names' instead.");
        prm.declare_entry ("Include adiabatic heating", "false",
                           Patterns::Bool (),
                           "Whether to include adiabatic heating into the model or not. From a "
                           "physical viewpoint, adiabatic heating should always be used but may "
                           "be undesirable when comparing results with known benchmarks that "
                           "do not include this term in the temperature equation."
                           "Warning: deprecated! Add 'adiabatic heating' to the 'List of model "
                           "names' instead.");
        prm.declare_entry ("Include latent heat", "false",
                           Patterns::Bool (),
                           "Whether to include the generation of latent heat at phase transitions "
                           "into the model or not. From a physical viewpoint, latent heat should "
                           "always be used but may be undesirable when comparing results with known "
                           "benchmarks that do not include this term in the temperature equation "
                           "or when dealing with a model without phase transitions."
                           "Warning: deprecated! Add 'latent heat' to the 'List of model "
                           "names' instead.");
      }
      prm.leave_subsection ();

      std_cxx11::get<dim>(registered_plugins).declare_parameters (prm);
    }


    HeatingModelOutputs::HeatingModelOutputs(const unsigned int n_points,
                                             const unsigned int)
      :
      heating_source_terms(n_points),
      lhs_latent_heat_terms(n_points)
    {
    }


    template <int dim>
    std::string
    get_valid_model_names_pattern ()
    {
      return std_cxx11::get<dim>(registered_plugins).get_pattern_of_names ();
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
      std::list<internal::Plugins::PluginList<HeatingModel::Interface<2> >::PluginInfo> *
      internal::Plugins::PluginList<HeatingModel::Interface<2> >::plugins = 0;
      template <>
      std::list<internal::Plugins::PluginList<HeatingModel::Interface<3> >::PluginInfo> *
      internal::Plugins::PluginList<HeatingModel::Interface<3> >::plugins = 0;
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
  }
}
