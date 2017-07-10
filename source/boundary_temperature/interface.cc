/*
  Copyright (C) 2011 - 2016 by the authors of the ASPECT code.

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
#include <aspect/boundary_temperature/interface.h>
#include <aspect/simulator_access.h>

#include <aspect/utilities.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/std_cxx11/tuple.h>

#include <list>


namespace aspect
{
  namespace BoundaryTemperature
  {
    template <int dim>
    Interface<dim>::~Interface ()
    {}


    template <int dim>
    void
    Interface<dim>::update ()
    {}

    template <int dim>
    void
    Interface<dim>::initialize ()
    {}

    template <int dim>
    double
    Interface<dim>::boundary_temperature (const types::boundary_id /*boundary_indicator*/,
                                          const Point<dim>        &/*position*/) const
    {
      AssertThrow(false,
                  ExcMessage("The boundary temperature plugin has to implement a function called `temperature' "
                             "with three arguments or a function `boundary_temperature' with two arguments. "
                             "The function with three arguments is deprecated and will "
                             "be removed in a later version of ASPECT."));
      return numbers::signaling_nan<double>();
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
    // -------------------------------- Deal with registering boundary_temperature models and automating
    // -------------------------------- their setup and selection at run time

    template <int dim>
    Manager<dim>::~Manager()
    {}



    template <int dim>
    void
    Manager<dim>::update ()
    {
      for (unsigned int i=0; i<boundary_temperature_objects.size(); ++i)
        {
          boundary_temperature_objects[i]->update();
        }
      return;
    }



    namespace
    {
      std_cxx11::tuple
      <void *,
      void *,
      aspect::internal::Plugins::PluginList<Interface<2> >,
      aspect::internal::Plugins::PluginList<Interface<3> > > registered_plugins;
    }


    template <int dim>
    void
    Manager<dim>::register_boundary_temperature (const std::string &name,
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
      prm.enter_subsection ("Boundary temperature model");
      {
        model_names
          = Utilities::split_string_list(prm.get("List of model names"));

        AssertThrow(Utilities::has_unique_entries(model_names),
                    ExcMessage("The list of strings for the parameter "
                               "'Boundary temperature model/List of model names' contains entries more than once. "
                               "This is not allowed. Please check your parameter file."));

        const std::string model_name = prm.get ("Model name");

        AssertThrow (model_name == "unspecified" || model_names.size() == 0,
                     ExcMessage ("The parameter 'Model name' is only used for reasons"
                                 "of backwards compatibility and can not be used together with "
                                 "the new functionality 'List of model names'. Please add your "
                                 "boundary temperature model to the list instead."));

        if (!(model_name == "unspecified"))
          model_names.push_back(model_name);

        // create operator list
        std::vector<std::string> model_operator_names =
          Utilities::possibly_extend_from_1_to_N (Utilities::split_string_list(prm.get("List of model operators")),
                                                  model_names.size(),
                                                  "List of model operators");
        model_operators = Utilities::create_model_operator_list(model_operator_names);
      }
      prm.leave_subsection ();

      // go through the list, create objects and let them parse
      // their own parameters
      for (unsigned int i=0; i<model_names.size(); ++i)
        {
          // create boundary temperature objects
          boundary_temperature_objects.push_back (std_cxx11::shared_ptr<Interface<dim> >
                                                  (std_cxx11::get<dim>(registered_plugins)
                                                   .create_plugin (model_names[i],
                                                                   "Boundary temperature::Model names")));

          if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(boundary_temperature_objects.back().get()))
            sim->initialize_simulator (this->get_simulator());

          boundary_temperature_objects.back()->parse_parameters (prm);
          boundary_temperature_objects.back()->initialize ();
        }
    }



    template <int dim>
    double
    Manager<dim>::boundary_temperature (const types::boundary_id boundary_indicator,
                                        const Point<dim> &position) const
    {
      double temperature = 0.0;

      for (unsigned int i=0; i<boundary_temperature_objects.size(); ++i)
        temperature = model_operators[i](temperature,
                                         boundary_temperature_objects[i]->boundary_temperature(boundary_indicator,
                                             position));

      return temperature;
    }



    template <int dim>
    double
    Manager<dim>::minimal_temperature (const std::set<types::boundary_id> &fixed_boundary_ids) const
    {
      double temperature = std::numeric_limits<double>::max();

      for (unsigned int i=0; i<boundary_temperature_objects.size(); ++i)
        temperature = std::min(temperature,
                               boundary_temperature_objects[i]->minimal_temperature(fixed_boundary_ids));

      return temperature;
    }



    template <int dim>
    double
    Manager<dim>::maximal_temperature (const std::set<types::boundary_id> &fixed_boundary_ids) const
    {
      double temperature = 0.0;

      for (unsigned int i=0; i<boundary_temperature_objects.size(); ++i)
        temperature = std::max(temperature,
                               boundary_temperature_objects[i]->maximal_temperature(fixed_boundary_ids));

      return temperature;
    }



    template <int dim>
    const std::vector<std::string> &
    Manager<dim>::get_active_boundary_temperature_names () const
    {
      return model_names;
    }


    template <int dim>
    const std::vector<std_cxx11::shared_ptr<Interface<dim> > > &
    Manager<dim>::get_active_boundary_temperature_conditions () const
    {
      return boundary_temperature_objects;
    }


    template <int dim>
    void
    Manager<dim>::declare_parameters (ParameterHandler &prm)
    {
      // declare the entry in the parameter file
      prm.enter_subsection ("Boundary temperature model");
      {
        const std::string pattern_of_names
          = std_cxx11::get<dim>(registered_plugins).get_pattern_of_names ();

        prm.declare_entry("List of model names",
                          "",
                          Patterns::MultipleSelection(pattern_of_names),
                          "A comma-separated list of boundary temperature models that "
                          "will be used to initialize the temperature. "
                          "These plugins are loaded in the order given, and modify the "
                          "existing temperature field via the operators listed "
                          "in 'List of model operators'.\n\n"
                          "The following boundary temperature models are available:\n\n"
                          +
                          std_cxx11::get<dim>(registered_plugins).get_description_string());

        prm.declare_entry("List of model operators", "add",
                          Patterns::MultipleSelection("add|subtract|minimum|maximum"),
                          "A comma-separated list of operators that "
                          "will be used to append the listed temperature models onto "
                          "the previous models. If only one operator is given, "
                          "the same operator is applied to all models.");

        prm.declare_entry ("Model name", "unspecified",
                           Patterns::Selection (pattern_of_names+"|unspecified"),
                           "Select one of the following models:\n\n"
                           +
                           std_cxx11::get<dim>(registered_plugins).get_description_string()
                           + "\n\n" +
                           "\\textbf{Warning}: This parameter provides an old and "
                           "deprecated way of specifying "
                           "boundary temperature models and shouldn't be used. "
                           "Please use 'List of model names' instead.");
      }
      prm.leave_subsection ();

      std_cxx11::get<dim>(registered_plugins).declare_parameters (prm);
    }



    template <int dim>
    void
    Manager<dim>::write_plugin_graph (std::ostream &out)
    {
      std_cxx11::get<dim>(registered_plugins).write_plugin_graph ("Boundary temperature interface",
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
      std::list<internal::Plugins::PluginList<BoundaryTemperature::Interface<2> >::PluginInfo> *
      internal::Plugins::PluginList<BoundaryTemperature::Interface<2> >::plugins = 0;
      template <>
      std::list<internal::Plugins::PluginList<BoundaryTemperature::Interface<3> >::PluginInfo> *
      internal::Plugins::PluginList<BoundaryTemperature::Interface<3> >::plugins = 0;
    }
  }

  namespace BoundaryTemperature
  {
#define INSTANTIATE(dim) \
  template class Interface<dim>; \
  template class Manager<dim>;

    ASPECT_INSTANTIATE(INSTANTIATE)
  }
}
