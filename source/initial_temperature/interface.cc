/*
  Copyright (C) 2011 - 2019 by the authors of the ASPECT code.

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
#include <aspect/initial_temperature/interface.h>

#include <deal.II/base/exceptions.h>
#include <tuple>

#include <list>


namespace aspect
{
  namespace InitialTemperature
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
    Interface<dim>::
    declare_parameters (dealii::ParameterHandler &)
    {}


    template <int dim>
    void
    Interface<dim>::parse_parameters (dealii::ParameterHandler &)
    {}



    // ------------------------------ Manager -----------------------------
    // -------------------------------- Deal with registering initial_temperature models and automating
    // -------------------------------- their setup and selection at run time

    template <int dim>
    Manager<dim>::~Manager()
    {}



    namespace
    {
      std::tuple
      <void *,
      void *,
      aspect::internal::Plugins::PluginList<Interface<2> >,
      aspect::internal::Plugins::PluginList<Interface<3> > > registered_plugins;
    }


    template <int dim>
    void
    Manager<dim>::register_initial_temperature (const std::string &name,
                                                const std::string &description,
                                                void (*declare_parameters_function) (ParameterHandler &),
                                                Interface<dim> *(*factory_function) ())
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
      prm.enter_subsection ("Initial temperature model");
      {
        model_names
          = Utilities::split_string_list(prm.get("List of model names"));

        AssertThrow(Utilities::has_unique_entries(model_names),
                    ExcMessage("The list of strings for the parameter "
                               "'Initial temperature model/List of model names' contains entries more than once. "
                               "This is not allowed. Please check your parameter file."));

        const std::string model_name = prm.get ("Model name");

        AssertThrow (model_name == "unspecified" || model_names.size() == 0,
                     ExcMessage ("The parameter 'Model name' is only used for reasons"
                                 "of backwards compatibility and can not be used together with "
                                 "the new functionality 'List of model names'. Please add your "
                                 "initial temperature model to the list instead."));

        if (!(model_name == "unspecified"))
          model_names.push_back(model_name);

        // create operator list
        const std::vector<std::string> model_operator_names =
          Utilities::possibly_extend_from_1_to_N (Utilities::split_string_list(prm.get("List of model operators")),
                                                  model_names.size(),
                                                  "List of model operators");
        model_operators = Utilities::create_model_operator_list(model_operator_names);
      }
      prm.leave_subsection ();

      // go through the list, create objects and let them parse
      // their own parameters
      for (auto &model_name : model_names)
        {
          // create initial temperature objects
          initial_temperature_objects.push_back (std::unique_ptr<Interface<dim> >
                                                 (std::get<dim>(registered_plugins)
                                                  .create_plugin (model_name,
                                                                  "Initial temperature model::Model names")));

          if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(&*initial_temperature_objects.back()))
            sim->initialize_simulator (this->get_simulator());

          initial_temperature_objects.back()->parse_parameters (prm);
          initial_temperature_objects.back()->initialize ();
        }
    }


    template <int dim>
    double
    Manager<dim>::initial_temperature (const Point<dim> &position) const
    {
      double temperature = 0.0;
      int i = 0;

      for (const auto &initial_temperature_object : initial_temperature_objects)
        {
          temperature = model_operators[i](temperature,
                                           initial_temperature_object->initial_temperature(position));
          i++;
        }
      return temperature;
    }


    template <int dim>
    const std::vector<std::string> &
    Manager<dim>::get_active_initial_temperature_names () const
    {
      return model_names;
    }


    template <int dim>
    const std::list<std::unique_ptr<Interface<dim> > > &
    Manager<dim>::get_active_initial_temperature_conditions () const
    {
      return initial_temperature_objects;
    }


    template <int dim>
    void
    Manager<dim>::declare_parameters (ParameterHandler &prm)
    {
      // declare the entry in the parameter file
      prm.enter_subsection ("Initial temperature model");
      {
        const std::string pattern_of_names
          = std::get<dim>(registered_plugins).get_pattern_of_names ();

        prm.declare_entry("List of model names",
                          "",
                          Patterns::MultipleSelection(pattern_of_names),
                          "A comma-separated list of initial temperature models that "
                          "will be used to initialize the temperature. "
                          "These plugins are loaded in the order given, and modify the "
                          "existing temperature field via the operators listed "
                          "in 'List of model operators'.\n\n"
                          "The following initial temperature models are available:\n\n"
                          +
                          std::get<dim>(registered_plugins).get_description_string());

        prm.declare_entry("List of model operators", "add",
                          Patterns::MultipleSelection(Utilities::get_model_operator_options()),
                          "A comma-separated list of operators that "
                          "will be used to append the listed temperature models onto "
                          "the previous models. If only one operator is given, "
                          "the same operator is applied to all models.");

        prm.declare_entry ("Model name", "unspecified",
                           Patterns::Selection (pattern_of_names+"|unspecified"),
                           "Select one of the following models:\n\n"
                           +
                           std::get<dim>(registered_plugins).get_description_string()
                           + "\n\n" +
                           "\\textbf{Warning}: This parameter provides an old and "
                           "deprecated way of specifying "
                           "initial temperature models and shouldn't be used. "
                           "Please use 'List of model names' instead.");
      }
      prm.leave_subsection ();

      std::get<dim>(registered_plugins).declare_parameters (prm);
    }



    template <int dim>
    std::string
    get_valid_model_names_pattern ()
    {
      return std::get<dim>(registered_plugins).get_pattern_of_names ();
    }



    template <int dim>
    void
    Manager<dim>::write_plugin_graph (std::ostream &out)
    {
      std::get<dim>(registered_plugins).write_plugin_graph ("Initial temperature interface",
                                                            out);
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
      std::list<internal::Plugins::PluginList<InitialTemperature::Interface<2> >::PluginInfo> *
      internal::Plugins::PluginList<InitialTemperature::Interface<2> >::plugins = nullptr;
      template <>
      std::list<internal::Plugins::PluginList<InitialTemperature::Interface<3> >::PluginInfo> *
      internal::Plugins::PluginList<InitialTemperature::Interface<3> >::plugins = nullptr;
    }
  }

  namespace InitialTemperature
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
