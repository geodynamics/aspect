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
#include <aspect/boundary_composition/interface.h>

#include <aspect/utilities.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/std_cxx11/tuple.h>

#include <list>


namespace aspect
{
  namespace BoundaryComposition
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
    // -------------------------------- Deal with registering boundary_composition models and automating
    // -------------------------------- their setup and selection at run time

    template <int dim>
    Manager<dim>::~Manager()
    {}



    template <int dim>
    void
    Manager<dim>::update ()
    {
      for (unsigned int i=0; i<boundary_composition_objects.size(); ++i)
        {
          boundary_composition_objects[i]->update();
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
    Manager<dim>::register_boundary_composition (const std::string &name,
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
      prm.enter_subsection ("Boundary composition model");
      {
        boundary_composition_names
          = Utilities::split_string_list(prm.get("List of model names"));

        AssertThrow(Utilities::has_unique_entries(boundary_composition_names),
                    ExcMessage("The list of strings for the parameter "
                               "'Boundary composition model/List of model names' contains entries more than once. "
                               "This is not allowed. Please check your parameter file."));

        // create operator list
        std::vector<std::string> model_operator_names =
          Utilities::possibly_extend_from_1_to_N (Utilities::split_string_list(prm.get("List of model operators")),
                                                  boundary_composition_names.size(),
                                                  "List of model operators");
        boundary_composition_operators = Utilities::create_model_operator_list(model_operator_names,
                                                                               "Boundary composition model/List of model operators");
      }
      prm.leave_subsection ();

      // Now, figure out which model is used for which boundary and which composition
      prm.enter_subsection("Model settings");
      {
        boundary_composition_maps.resize(this->n_compositional_fields());

        const std::vector<std::string> x_fixed_composition_boundary_indicators
          = Utilities::split_string_list (prm.get ("Fixed composition boundary indicators"));

        for (unsigned int i = 0; i < x_fixed_composition_boundary_indicators.size(); ++i)
          {
            // each entry has the format (white space is optional):
            // <boundary_id>
            // or
            // <boundary_id> : <plugin_name>
            // or
            // <boundary_id> : <composition_id>
            // or
            // <boundary_id> : <composition_id> : <plugin_name>.
            // Parse the entry. We first tease apart the list.
            const std::vector<std::string> split_parts = Utilities::split_string_list (x_fixed_composition_boundary_indicators[i], ':');

            // try to translate the first part into a boundary_id
            types::boundary_id boundary_id;
            try
              {
                boundary_id = this->get_geometry_model().translate_symbolic_boundary_name_to_id(split_parts[0]);
                fixed_composition_boundary_indicators.insert(boundary_id);
              }
            catch (const std::string &error)
              {
                AssertThrow (false, ExcMessage ("While parsing the entry <Model settings/Prescribed "
                                                "composition indicators>, there was an error. Specifically, "
                                                "the conversion function complained as follows: "
                                                + error));
              }

            // If only a boundary_id is given, use all models for all compositions for this boundary
            if (split_parts.size() == 1)
              {
                for (unsigned int composition = 0; composition < this->n_compositional_fields(); ++composition)
                  {
                    AssertThrow(boundary_composition_maps[composition].find(boundary_id) == boundary_composition_maps[composition].end(),
                                ExcMessage("You can not specify a single boundary indicator in the "
                                           "<Model settings/Prescribed composition indicators> parameter "
                                           "if you already specified something for the same boundary. "
                                           "For each boundary, either use all "
                                           "plugins by only specifying the boundary indicator, or select one or "
                                           "several plugins by name."));

                    for (unsigned int j = 0; j < boundary_composition_names.size(); ++j)
                      boundary_composition_maps[composition].insert(std::make_pair(boundary_id,j));
                  }
              }
            else if (split_parts.size() >= 2)
              {
                // Check if the two parts are something like "boundary_id : plugin_name"
                std::vector<std::string>::iterator plugin_name = std::find(boundary_composition_names.begin(),
                                                                           boundary_composition_names.end(),
                                                                           split_parts[1]);

                if (plugin_name != boundary_composition_names.end())
                  {
                    // The second part is a plugin name. Add this plugin for all compositions.
                    const unsigned int plugin_position = std::distance(boundary_composition_names.begin(),
                                                                       plugin_name);

                    // Make sure this plugin was not already added for this boundary
                    for (unsigned int composition = 0; composition < this->n_compositional_fields(); ++composition)
                      {
                        std::pair<std::multimap<types::boundary_id,unsigned int>::const_iterator,
                            std::multimap<types::boundary_id,unsigned int>::const_iterator> relevant_plugins =
                              boundary_composition_maps[composition].equal_range(boundary_id);

                        for (; relevant_plugins.first != relevant_plugins.second; ++relevant_plugins.first)
                          {
                            const unsigned int existing_plugin_position = relevant_plugins.first->second;
                            AssertThrow(plugin_position != existing_plugin_position,
                                        ExcMessage("You specified the same plugin multiple times for "
                                                   "the same boundary in the <Model settings/Fixed "
                                                   "composition boundary indicators>. This is not supported."));
                          }

                        boundary_composition_maps[composition].insert(std::make_pair(boundary_id,plugin_position));
                      }
                  }
                else
                  {
                    // The second part is a compositional field. Find the compositional index
                    // for the field. This call will fail if the string
                    // does not identify a compositional field.
                    const unsigned int composition = this->introspection().compositional_index_for_string(split_parts[1]);

                    if (split_parts.size() == 2)
                      {
                        // No specific plugins assigned, so add all of them, make sure nothing was assigned before.
                        AssertThrow(boundary_composition_maps[composition].count(boundary_id) == 0,
                                    ExcMessage("You specified all plugins for composition " + split_parts[1] +
                                               " for boundary " + split_parts[0] +
                                               " in the <Model settings/Fixed "
                                               "composition boundary indicators>, but this composition already has "
                                               "boundary conditions assigned for this boundary. Check your "
                                               "input file for double assignments."));

                        // Add the plugins
                        for (unsigned int plugin_position = 0; plugin_position != boundary_composition_names.size(); ++plugin_position)
                          boundary_composition_maps[composition].insert(std::make_pair(boundary_id,plugin_position));
                      }
                    else if (split_parts.size() == 3)
                      {
                        // The third part specified a plugin. Look for it and add it.
                        std::vector<std::string>::iterator plugin_name = std::find(boundary_composition_names.begin(),
                                                                                   boundary_composition_names.end(),
                                                                                   split_parts[2]);

                        AssertThrow(plugin_name != boundary_composition_names.end(),
                                    ExcMessage("The plugin name " + split_parts[2] + " in the parameter "
                                               "<Model settings/Fixed composition boundary indicators> has to be one of "
                                               "the selected plugin names in the parameter "
                                               "<Boundary composition model/List of model names>. This seems to be not the "
                                               "case. Please check your input file."));

                        const unsigned int plugin_position = std::distance(boundary_composition_names.begin(),
                                                                           plugin_name);

                        std::pair<std::multimap<types::boundary_id,unsigned int>::const_iterator,
                            std::multimap<types::boundary_id,unsigned int>::const_iterator> relevant_plugins =
                              boundary_composition_maps[composition].equal_range(boundary_id);

                        for (; relevant_plugins.first != relevant_plugins.second; ++relevant_plugins.first)
                          {
                            const unsigned int existing_plugin_position = relevant_plugins.first->second;
                            AssertThrow(plugin_position != existing_plugin_position,
                                        ExcMessage("You specified the same plugin multiple times for "
                                                   "the same boundary in the <Model settings/Fixed "
                                                   "composition boundary indicators>. This is not supported."));
                          }

                        boundary_composition_maps[composition].insert(std::make_pair(boundary_id,plugin_position));
                      }
                  }
              }
          }
      }
      prm.leave_subsection();

      // For now we can not support boundary conditions for one field, but not for others at the same
      // boundary, because we use a single matrix block for all of them, and different Dirichlet
      // boundary conditions require a different sparsity pattern. Thus for now all boundary_composition_maps
      // need to have the same set of keys.
      //
      // TODO: Lift this restriction by keeping constrained degrees of freedom in the sparsity pattern
      // of the composition matrix block (might require changes to deal.II).
      for (std::set<types::boundary_id>::const_iterator boundary = fixed_composition_boundary_indicators.begin();
           boundary != fixed_composition_boundary_indicators.end(); ++boundary)
        {
          for (unsigned int i=0; i<boundary_composition_maps.size(); ++i)
            AssertThrow(boundary_composition_maps[i].find(*boundary) != boundary_composition_maps[i].end(),
                        ExcMessage("You have specified Dirichlet boundary conditions for any compositional field "
                                   "at boundary id <" + this->get_geometry_model().translate_id_to_symbol_name(*boundary)
                                   + "> but not for compositional field number "
                                   + Utilities::to_string(i) + ". If you specify Dirichlet boundary conditions for a "
                                   "compositional field at a boundary you have to do so for every field."));
        }

      // Now that everything is parsed, go through the list, create objects
      // and let them parse their own parameters.
      for (unsigned int i=0; i<boundary_composition_names.size(); ++i)
        {
          boundary_composition_objects.push_back (std_cxx11::shared_ptr<Interface<dim> >
                                                  (std_cxx11::get<dim>(registered_plugins)
                                                   .create_plugin (boundary_composition_names[i],
                                                                   "Boundary composition::Model names")));

          if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(boundary_composition_objects.back().get()))
            sim->initialize_simulator (this->get_simulator());

          boundary_composition_objects.back()->parse_parameters (prm);
          boundary_composition_objects.back()->initialize ();
        }
    }



    template <int dim>
    double
    Manager<dim>::boundary_composition (const types::boundary_id boundary_indicator,
                                        const Point<dim> &position,
                                        const unsigned int compositional_field) const
    {
      double composition = 0.0;

      std::pair<std::multimap<types::boundary_id,unsigned int>::const_iterator,
          std::multimap<types::boundary_id,unsigned int>::const_iterator> relevant_plugins =
            boundary_composition_maps[compositional_field].equal_range(boundary_indicator);

      for (; relevant_plugins.first != relevant_plugins.second; ++relevant_plugins.first)
        {
          const unsigned int plugin_index = relevant_plugins.first->second;
          const double plugin_composition = boundary_composition_objects[plugin_index]->boundary_composition(boundary_indicator,
                                            position,
                                            compositional_field);

          composition =
            boundary_composition_operators[plugin_index](composition,
                                                         plugin_composition);
        }

      return composition;
    }



    template <int dim>
    bool
    Manager<dim>::has_boundary_composition(const types::boundary_id boundary_id,
                                           const unsigned int       compositional_index) const
    {
      Assert(compositional_index < this->n_compositional_fields(),
             ExcMessage("There is no compositional field with index " + Utilities::to_string(compositional_index)));

      return (boundary_composition_maps[compositional_index].find(boundary_id) !=
              boundary_composition_maps[compositional_index].end());
    }



    template <int dim>
    const std::set<types::boundary_id> &
    Manager<dim>::get_fixed_composition_boundary_indicators() const
    {
      return fixed_composition_boundary_indicators;
    }



    template <int dim>
    const std::vector<std::string> &
    Manager<dim>::get_active_boundary_composition_names () const
    {
      return boundary_composition_names;
    }


    template <int dim>
    const std::vector<std_cxx11::shared_ptr<Interface<dim> > > &
    Manager<dim>::get_active_boundary_composition_conditions () const
    {
      return boundary_composition_objects;
    }


    template <int dim>
    void
    Manager<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Model settings");
      {
        prm.declare_entry ("Fixed composition boundary indicators", "",
                           Patterns::List (Patterns::Anything()),
                           "A comma separated list of names denoting those boundaries "
                           "on which the composition is fixed and described by the "
                           "boundary composition object selected in its own section "
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
                           "This parameter only describes which boundaries have a fixed "
                           "composition, but not what composition should hold on these "
                           "boundaries. The latter piece of information needs to be "
                           "implemented in a plugin in the BoundaryComposition "
                           "group, unless an existing implementation in this group "
                           "already provides what you want.");
      }
      prm.leave_subsection();

      // declare the entry in the parameter file
      prm.enter_subsection ("Boundary composition model");
      {
        const std::string pattern_of_names
          = std_cxx11::get<dim>(registered_plugins).get_pattern_of_names ();

        prm.declare_entry("List of model names",
                          "",
                          Patterns::MultipleSelection(pattern_of_names),
                          "A comma-separated list of boundary composition models that "
                          "will be used to initialize the composition. "
                          "These plugins are loaded in the order given, and modify the "
                          "existing composition field via the operators listed "
                          "in 'List of model operators'.\n\n"
                          "The following boundary composition models are available:\n\n"
                          +
                          std_cxx11::get<dim>(registered_plugins).get_description_string());

        prm.declare_entry("List of model operators", "add",
                          Patterns::MultipleSelection("add|subtract|minimum|maximum"),
                          "A comma-separated list of operators that "
                          "will be used to append the listed composition models onto "
                          "the previous models. If only one operator is given, "
                          "the same operator is applied to all models.");
      }
      prm.leave_subsection ();

      std_cxx11::get<dim>(registered_plugins).declare_parameters (prm);
    }



    template <int dim>
    void
    Manager<dim>::write_plugin_graph (std::ostream &out)
    {
      std_cxx11::get<dim>(registered_plugins).write_plugin_graph ("Boundary composition interface",
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
      std::list<internal::Plugins::PluginList<BoundaryComposition::Interface<2> >::PluginInfo> *
      internal::Plugins::PluginList<BoundaryComposition::Interface<2> >::plugins = 0;
      template <>
      std::list<internal::Plugins::PluginList<BoundaryComposition::Interface<3> >::PluginInfo> *
      internal::Plugins::PluginList<BoundaryComposition::Interface<3> >::plugins = 0;
    }
  }

  namespace BoundaryComposition
  {
#define INSTANTIATE(dim) \
  template class Interface<dim>; \
  template class Manager<dim>;

    ASPECT_INSTANTIATE(INSTANTIATE)
  }
}
