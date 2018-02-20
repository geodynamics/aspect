/*
  Copyright (C) 2011 - 2018 by the authors of the ASPECT code.

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
    void
    Interface<dim>::
    declare_parameters (ParameterHandler &)
    {}



    template <int dim>
    void
    Interface<dim>::parse_parameters (ParameterHandler &)
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
    }



    template <int dim>
    const std::set<types::boundary_id> &
    Manager<dim>::get_fixed_temperature_boundary_indicators() const
    {
      return fixed_temperature_boundary_indicators;
    }



    template <int dim>
    const std::vector<std::string> &
    Manager<dim>::get_active_boundary_temperature_names () const
    {
      return boundary_temperature_names;
    }



    template <int dim>
    const std::vector<std_cxx11::shared_ptr<Interface<dim> > > &
    Manager<dim>::get_active_boundary_temperature_conditions () const
    {
      return boundary_temperature_objects;
    }



    template <int dim>
    double
    Manager<dim>::minimal_temperature (const std::set<types::boundary_id> &fixed_boundary_ids) const
    {
      const std::set<types::boundary_id> &boundary_ids = fixed_boundary_ids.empty()
                                                         ?
                                                         fixed_temperature_boundary_indicators
                                                         :
                                                         fixed_boundary_ids;

      double temperature = std::numeric_limits<double>::max();
      for (typename std::vector<std_cxx11::shared_ptr<Interface<dim> > >::const_iterator
           boundary = boundary_temperature_objects.begin();
           boundary != boundary_temperature_objects.end(); ++boundary)
        temperature = std::min(temperature,
                               (*boundary)->minimal_temperature(boundary_ids));

      return temperature;
    }



    template <int dim>
    double
    Manager<dim>::maximal_temperature (const std::set<types::boundary_id> &fixed_boundary_ids) const
    {
      const std::set<types::boundary_id> &boundary_ids = fixed_boundary_ids.empty()
                                                         ?
                                                         fixed_temperature_boundary_indicators
                                                         :
                                                         fixed_boundary_ids;

      double temperature = 0.0;
      for (typename std::vector<std_cxx11::shared_ptr<Interface<dim> > >::const_iterator
           boundary = boundary_temperature_objects.begin();
           boundary != boundary_temperature_objects.end(); ++boundary)
        temperature = std::max(temperature,
                               (*boundary)->maximal_temperature(boundary_ids));

      return temperature;
    }



    template <int dim>
    double
    Manager<dim>::boundary_temperature (const types::boundary_id boundary_indicator,
                                        const Point<dim> &position) const
    {
      double temperature = 0.0;

      std::pair<std::multimap<types::boundary_id,unsigned int>::const_iterator,
          std::multimap<types::boundary_id,unsigned int>::const_iterator> relevant_plugins =
            boundary_temperature_map.equal_range(boundary_indicator);

      for (; relevant_plugins.first != relevant_plugins.second; ++relevant_plugins.first)
        {
          const unsigned int plugin_index = relevant_plugins.first->second;
          temperature = boundary_temperature_operators[plugin_index](temperature,
                                                                     boundary_temperature_objects[plugin_index]->boundary_temperature(boundary_indicator,
                                                                         position));
        }

      return temperature;
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
                          "will be used to set the boundary temperature. "
                          "These plugins are loaded in the order given, and modify the "
                          "boundary temperature via the operators listed "
                          "in 'List of model operators'. Which boundary model "
                          "is used for which boundary is controlled by the "
                          "'Boundary temperature model/Fixed temperature boundary indicators' parameter.\n\n"
                          "The following boundary temperature models are available:\n\n"
                          +
                          std_cxx11::get<dim>(registered_plugins).get_description_string());

        prm.declare_entry("List of model operators", "add",
                          Patterns::MultipleSelection("add|subtract|minimum|maximum"),
                          "A comma-separated list of operators that "
                          "will be used to append the listed temperature models onto "
                          "the previous models. If only one operator is given, "
                          "the same operator is applied to all models.");

        prm.declare_entry ("Fixed temperature boundary indicators", "",
                           Patterns::List (Patterns::Anything()),
                           "A comma separated list of names denoting those boundaries "
                           "on which the temperature is fixed and described by the "
                           "boundary temperature objects selected in the "
                           "<Boundary temperature model/List of model names> parameter "
                           "of this input file. All boundary indicators used by the geometry "
                           "but not explicitly listed here will end up with no-flux "
                           "(insulating) boundary conditions."
                           "\n\n"
                           "The names of the boundaries listed here can either by "
                           "numbers (in which case they correspond to the numerical "
                           "boundary indicators assigned by the geometry object), or they "
                           "can correspond to any of the symbolic names the geometry object "
                           "may have provided for each part of the boundary. You may want "
                           "to compare this with the documentation of the geometry model you "
                           "use in your model."
                           "\n\n"
                           "Each entry in this list has to be either a single boundary "
                           "identifier as described above, or a pattern of the type "
                           "'boundary identifier : plugin name', where plugin name is "
                           "one of the names given in the "
                           "<Boundary temperature model/List of model names> parameter. "
                           "If only a boundary identifier is given, the boundary temperature "
                           "of all active boundary temperature plugins will be combined to "
                           "compute the boundary temperature. If only one plugin is given, "
                           "only this plugin is used. The same boundary identifier can "
                           "appear multiple times with different plugin names, in which "
                           "case all given plugins will be combined to compute the boundary "
                           "temperature.");
      }
      prm.leave_subsection ();

      std_cxx11::get<dim>(registered_plugins).declare_parameters (prm);
    }



    template <int dim>
    void
    Manager<dim>::parse_parameters (ParameterHandler &prm)
    {
      // find out which plugins are requested and the various other
      // parameters we declare here
      prm.enter_subsection ("Boundary temperature model");
      {
        boundary_temperature_names
          = Utilities::split_string_list(prm.get("List of model names"));

        AssertThrow(Utilities::has_unique_entries(boundary_temperature_names),
                    ExcMessage("The list of strings for the parameter "
                               "'Boundary temperature model/List of model names' contains entries more than once. "
                               "This is not allowed. Please check your parameter file."));

        // create operator list
        boundary_temperature_operators = Utilities::create_model_operator_list(
                                           Utilities::possibly_extend_from_1_to_N (
                                             Utilities::split_string_list(prm.get("List of model operators")),
                                             boundary_temperature_names.size(),
                                             "List of model operators"),
                                           "Boundary temperature model/List of model operators");

        // Finally, figure out which model is used for which boundary
        const std::vector<std::string> x_fixed_temperature_boundary_indicators
          = Utilities::split_string_list (prm.get ("Fixed temperature boundary indicators"));

        for (unsigned int i = 0; i < x_fixed_temperature_boundary_indicators.size(); ++i)
          {
            // each entry has the format (white space is optional):
            // <id> : <value (might have spaces)>
            //
            // first tease apart the two halves
            const std::vector<std::string> split_parts = Utilities::split_string_list (x_fixed_temperature_boundary_indicators[i], ':');

            // try to translate the key into a boundary_id
            types::boundary_id boundary_id;
            try
              {
                boundary_id = this->get_geometry_model().translate_symbolic_boundary_name_to_id(split_parts[0]);
                fixed_temperature_boundary_indicators.insert(boundary_id);
              }
            catch (const std::string &error)
              {
                AssertThrow (false, ExcMessage ("While parsing the entry <Boundary temperature model/Prescribed "
                                                "temperature indicators>, there was an error. Specifically, "
                                                "the conversion function complained as follows: "
                                                + error));
              }

            // If only a boundary_id is given, use all models for this boundary
            if (split_parts.size() == 1)
              {
                AssertThrow(boundary_temperature_map.find(boundary_id) == boundary_temperature_map.end(),
                            ExcMessage("You can not specify a single boundary indicator in the "
                                       "<Boundary temperature model/Prescribed temperature indicators> parameter "
                                       "if you already specified something for the same boundary. "
                                       "For each boundary, either use all "
                                       "plugins by only specifying the boundary indicator, or select one or "
                                       "several plugins by name."));

                for (unsigned int j = 0; j < boundary_temperature_names.size(); ++j)
                  boundary_temperature_map.insert(std::make_pair(boundary_id,j));
              }
            // Add the selected plugin otherwise.
            else if (split_parts.size() == 2)
              {
                // Make sure the two parts are something like "boundary_id : plugin_name"
                // We checked the first part above already, now check the second.
                std::vector<std::string>::iterator plugin_name = std::find(boundary_temperature_names.begin(),
                                                                           boundary_temperature_names.end(),
                                                                           split_parts[1]);

                AssertThrow(plugin_name != boundary_temperature_names.end(),
                            ExcMessage("The plugin name " + split_parts[1] + " in the parameter "
                                       "<Boundary temperature model/Fixed temperature boundary indicators> has to be one of "
                                       "the selected plugin names in the parameter "
                                       "<Boundary temperature model/List of model names>. This seems to be not the "
                                       "case. Please check your input file."));

                const unsigned int plugin_position = std::distance(boundary_temperature_names.begin(),
                                                                   plugin_name);

                // Make sure this plugin was not already added for this boundary
                {
                  std::pair<std::multimap<types::boundary_id,unsigned int>::const_iterator,
                      std::multimap<types::boundary_id,unsigned int>::const_iterator> relevant_plugins =
                        boundary_temperature_map.equal_range(boundary_id);

                  for (; relevant_plugins.first != relevant_plugins.second; ++relevant_plugins.first)
                    {
                      const unsigned int existing_plugin_position = relevant_plugins.first->second;
                      AssertThrow(plugin_position != existing_plugin_position,
                                  ExcMessage("You specified the same plugin multiple times for "
                                             "the same boundary in the <Boundary temperature model/Fixed "
                                             "temperature boundary indicators>. This is not supported."));
                    }
                }

                // Add the plugin
                boundary_temperature_map.insert(std::make_pair(boundary_id,plugin_position));
              }
            else
              AssertThrow(false,
                          ExcMessage("Each entry in <Boundary temperature model/Fixed "
                                     "temperature boundary indicators> has to be either "
                                     "a boundary name or a pair of boundary name : plugin name"));
          }
      }
      prm.leave_subsection();

      // Now create and initialize all of the requested models. We can only do this here at the end,
      // because some models might ask, which boundary they are responsible for.
      for (unsigned int i = 0; i<boundary_temperature_names.size(); ++i)
        {
          boundary_temperature_objects.push_back(std_cxx11::shared_ptr<Interface<dim> > (std_cxx11::get<dim>(registered_plugins)
                                                                                         .create_plugin (boundary_temperature_names[i],
                                                                                             "Boundary temperature::Model names")));
          if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(boundary_temperature_objects.back().get()))
            sim->initialize_simulator (this->get_simulator());
          boundary_temperature_objects.back()->parse_parameters (prm);
          boundary_temperature_objects.back()->initialize ();
        }
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
