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
#include <aspect/utilities.h>
#include <aspect/initial_composition/interface.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/std_cxx11/tuple.h>

#include <list>


namespace aspect
{
  namespace InitialComposition
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
    // ------------------------------ Deal with registering initial composition models and automating
    // ------------------------------ their setup and selection at run time

    template <int dim>
    Manager<dim>::~Manager()
    {}



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
    Manager<dim>::register_initial_composition (const std::string &name,
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
      prm.enter_subsection ("Initial composition model");
      {
        model_names
          = Utilities::split_string_list(prm.get("List of model names"));

        AssertThrow(Utilities::has_unique_entries(model_names),
                    ExcMessage("The list of strings for the parameter "
                               "'Compositional initial conditions/List of model names' contains entries more than once. "
                               "This is not allowed. Please check your parameter file."));

        const std::string model_name = prm.get ("Model name");

        AssertThrow (model_name == "unspecified" || model_names.size() == 0,
                     ExcMessage ("The parameter 'Model name' is only used for reasons"
                                 "of backwards compatibility and can not be used together with "
                                 "the new functionality 'List of model names'. Please add your "
                                 "initial composition model to the list instead."));

        if (!(model_name == "unspecified"))
          model_names.push_back(model_name);

      }
      prm.leave_subsection ();

      if (model_names.size() > 0)
        AssertThrow(this->n_compositional_fields() > 0,
                    ExcMessage("A plugin for the initial composition condition was specified, but there "
                               "is no compositional field. This can lead to errors within the initialization of "
                               "the initial composition plugin and is therefore not supported. Please remove "
                               "the initial composition plugin or add a compositional field."));

      // go through the list, create objects and let them parse
      // their own parameters
      for (unsigned int name=0; name<model_names.size(); ++name)
        {
          initial_composition_objects.push_back (std_cxx11::shared_ptr<Interface<dim> >
                                                 (std_cxx11::get<dim>(registered_plugins)
                                                  .create_plugin (model_names[name],
                                                                  "Compositional initial conditions::Model names")));

          if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(&*initial_composition_objects.back()))
            sim->initialize_simulator (this->get_simulator());

          initial_composition_objects.back()->parse_parameters (prm);
          initial_composition_objects.back()->initialize ();
        }
    }



    template <int dim>
    double
    Manager<dim>::initial_composition (const Point<dim> &position,
                                       const unsigned int n_comp) const
    {
      double composition = 0.0;
      for (typename std::list<std_cxx11::shared_ptr<InitialComposition::Interface<dim> > >::const_iterator
           initial_composition_object = initial_composition_objects.begin();
           initial_composition_object != initial_composition_objects.end();
           ++initial_composition_object)
        {
          composition += (*initial_composition_object)->initial_composition(position,n_comp);
        }
      return composition;
    }


    template <int dim>
    const std::vector<std::string> &
    Manager<dim>::get_active_initial_composition_names () const
    {
      return model_names;
    }


    template <int dim>
    const std::list<std_cxx11::shared_ptr<Interface<dim> > > &
    Manager<dim>::get_active_initial_composition_conditions () const
    {
      return initial_composition_objects;
    }


    template <int dim>
    void
    Manager<dim>::declare_parameters (ParameterHandler &prm)
    {
      // declare the entry in the parameter file
      prm.enter_subsection ("Initial composition model");
      {
        const std::string pattern_of_names
          = std_cxx11::get<dim>(registered_plugins).get_pattern_of_names ();

        prm.declare_entry("List of model names",
                          "",
                          Patterns::MultipleSelection(pattern_of_names),
                          "A comma separated list of initial composition models that "
                          "describe the initial composition field. The results "
                          "of each of these criteria will be added.\n\n"
                          "The following heating models are available:\n\n"
                          +
                          std_cxx11::get<dim>(registered_plugins).get_description_string());

        prm.declare_entry ("Model name", "unspecified",
                           Patterns::Selection (pattern_of_names+"|unspecified"),
                           "Select one of the following models:\n\n"
                           "Warning: This is the old formulation of specifying "
                           "initial composition models and shouldn't be used. "
                           "Please use 'List of model names' instead."
                           +
                           std_cxx11::get<dim>(registered_plugins).get_description_string());
      }
      prm.leave_subsection ();

      std_cxx11::get<dim>(registered_plugins).declare_parameters (prm);
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
      std::list<internal::Plugins::PluginList<InitialComposition::Interface<2> >::PluginInfo> *
      internal::Plugins::PluginList<InitialComposition::Interface<2> >::plugins = 0;
      template <>
      std::list<internal::Plugins::PluginList<InitialComposition::Interface<3> >::PluginInfo> *
      internal::Plugins::PluginList<InitialComposition::Interface<3> >::plugins = 0;
    }
  }

  namespace InitialComposition
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
