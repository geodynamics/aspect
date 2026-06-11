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
#include <aspect/boundary_composition/interface.h>

#include <aspect/utilities.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/signaling_nan.h>

#include <tuple>
#include <list>


namespace aspect
{
  namespace BoundaryComposition
  {
    // ------------------------------ Manager -----------------------------
    // -------------------------------- Deal with registering boundary_composition models and automating
    // -------------------------------- their setup and selection at run time

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
    Manager<dim>::register_boundary_composition (const std::string &name,
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
      prm.enter_subsection ("Boundary composition model");
      {
        // First split the list into its separate entries that each
        // specify a boundary indicator, and, potentially, the plugins that are active
        // on that boundary for the listed compositional fields.
        const std::vector<std::string> x_fixed_composition_boundary_indicators
          = Utilities::split_string_list(prm.get("Fixed composition boundary indicators"), ',');

        // For each plugin object, for each boundary it is active on,
        // we want to get a component_mask that specifies for which
        // compositional fields it is active (`masks_fields`).
        // We also need a vector of boundary indicators for each plugin
        // that specifies on which boundaries it is active.
        for (const auto &p : x_fixed_composition_boundary_indicators)
          {
            // Each entry has the format:
            // <boundary_name>
            // or
            // <boundary_name> [field_name_1|field_name_2|field_name_5] : plugin_name_1+plugin_name_2
            // That means
            // <boundary_name> [field_name_1|field_name_2]
            // is not allowed; plugins need to be specified if field names are specified.
            //
            // First tease apart the two halves, if there are two.
            const std::vector<std::string> split_parts = Utilities::split_string_list (p, ':');

            // If there is only one part, it should be a boundary name.
            if (split_parts.size() == 1)
              {
                types::boundary_id boundary_indicator;
                try
                  {
                    boundary_indicator = this->get_geometry_model().translate_symbolic_boundary_name_to_id(split_parts[0]);
                  }
                catch (const std::string &error)
                  {
                    AssertThrow (false, ExcMessage ("While parsing the entry <Boundary composition model/Fixed "
                                                    "boundary composition indicators>, there was an error. Specifically, "
                                                    "the conversion function complained as follows:\n\n"
                                                    + error));
                  }

                // Since only the boundary indicator is specified,
                // all fields need to be fixed on this boundary.
                ComponentMask component_mask(this->n_compositional_fields(),
                                             true);

                for (unsigned int f = 0; f<this->n_compositional_fields(); ++f)
                  fixed_compositional_fields.insert(f);

                // The plugins that should be used to compute the
                // fixed composition are given as one list.
                // Any boundary for which no fields and plugins are
                // specified, will get the same list of plugins.
                // TODO instead of rereading, copy?
                std::vector<std::string> list_of_plugins = Utilities::split_string_list(prm.get("List of model names"));

                AssertThrow(Utilities::has_unique_entries(list_of_plugins),
                            ExcMessage("The list of strings for the parameter "
                                       "'Boundary composition model/List of model names' contains entries more than once. "
                                       "This is not allowed. Please check your parameter file."));

                // Create the operator list to combine the plugin results.
                std::vector<std::string> model_operator_names =
                  Utilities::possibly_extend_from_1_to_N (Utilities::split_string_list(prm.get("List of model operators")),
                                                          list_of_plugins.size(),
                                                          "List of model operators");

                std::vector<aspect::Utilities::Operator> list_of_model_operators = Utilities::create_model_operator_list(model_operator_names);

                // Loop over all the listed plugins and check whether
                // this plugin was already listed for another boundary/field.
                unsigned int mn = 0;
                for (auto &model_name : list_of_plugins)
                  {
                    // Multiple boundaries can list the same plugins.
                    // We do not need multiple entries in plugin_names, however.
                    const auto plugin_name_it = std::find(this->plugin_names.begin(), this->plugin_names.end(), model_name);
                    if (plugin_name_it == this->plugin_names.end())
                      {
                        // The plugin was not listed before, so add it.
                        this->plugin_names.push_back(model_name);

                        // Also add the compositional field masks for this plugin.
                        std::vector<ComponentMask> vec_component_mask;
                        vec_component_mask.push_back(component_mask);
                        masks_fields.push_back(vec_component_mask);

                        // Insert the boundary indicator that belongs to this plugin into the set.
                        fixed_composition_boundary_indicators.insert(boundary_indicator);

                        // Insert the boundary indicator that belongs to this plugin into the vector.
                        // It is possible that more boundary indicators are affected by this plugin,
                        // they would be added to the vector of indicators for this plugin
                        // in the other conditional branch.
                        std::vector<types::boundary_id> vec_boundary_indicator;
                        vec_boundary_indicator.push_back(boundary_indicator);
                        boundary_indicators.push_back(vec_boundary_indicator);

                        // Add the model operator for this plugin
                        std::vector<aspect::Utilities::Operator> vec_model_operator;
                        vec_model_operator.push_back(list_of_model_operators[mn]);
                        model_operators.push_back(vec_model_operator);
                      }
                    // The plugin was already listed for another boundary or field, so add the
                    // boundary indicator and corresponding masks at the existing plugin entry.
                    else
                      {
                        // Check that the plugin does not specify the same field at the same boundary.
                        for (unsigned int bi = 0; bi < boundary_indicators[plugin_name_it - this->plugin_names.begin()].size(); ++bi)
                          if (boundary_indicators[plugin_name_it - this->plugin_names.begin()][bi] == boundary_indicator)
                            for (unsigned int mf = 0; mf < masks_fields[plugin_name_it - this->plugin_names.begin()][bi].size(); ++mf)
                              if (masks_fields[plugin_name_it - this->plugin_names.begin()][bi][mf] == true)
                                AssertThrow (masks_fields[plugin_name_it - this->plugin_names.begin()][bi][mf] != component_mask[mf], ExcMessage ("The same plugin already specifies the same field on the same boundary."));

                        // Add the masks for this plugin and this boundary indicator.
                        masks_fields[plugin_name_it - this->plugin_names.begin()].push_back(component_mask);

                        // Add the boundary indicator for this plugin.
                        boundary_indicators[plugin_name_it - this->plugin_names.begin()].push_back(boundary_indicator);

                        // Also insert this boundary indicator for the existing plugin.
                        fixed_composition_boundary_indicators.insert(boundary_indicator);

                        // Add the operator for this boundary for the existing plugin.
                        model_operators[plugin_name_it - this->plugin_names.begin()].push_back(list_of_model_operators[mn]);
                      }
                    ++mn;
                  }

              }
            // If there are two parts, the second part holds plugin names.
            else if (split_parts.size() == 2)
              {
                // Get the plugin names from the second part.
                const std::vector<std::string> plugins = Utilities::split_string_list(split_parts[1], '+');

                // All the plugins get the operator 'add'.
                // TODO This reduces functionality of the operators.
                std::vector<std::string> model_operator_names =
                  Utilities::possibly_extend_from_1_to_N (Utilities::split_string_list("add"),
                                                          plugins.size(),
                                                          "add");

                std::vector<aspect::Utilities::Operator> list_of_model_operators = Utilities::create_model_operator_list(model_operator_names);

                // Split the first part into boundary indicator and compositional field names.
                const std::vector<std::string> first_split_parts = Utilities::split_string_list (split_parts[0], '[');

                AssertThrow (first_split_parts.size() == 2, ExcMessage ("While parsing the entry <Boundary composition model/Fixed "
                                                                        "boundary composition indicators>, there was an error. Specifically, "
                                                                        "the boundary indicator and corresponding plugins were incorrect in "
                                                                        + split_parts[0] + "."));

                types::boundary_id boundary_indicator;
                try
                  {
                    boundary_indicator = this->get_geometry_model().translate_symbolic_boundary_name_to_id(first_split_parts[0]);
                  }
                catch (const std::string &error)
                  {
                    AssertThrow (false, ExcMessage ("While parsing the entry <Boundary composition model/Fixed "
                                                    "boundary composition indicators>, there was an error. Specifically, "
                                                    "the conversion function complained as follows:\n\n"
                                                    + error));
                  }

                // Get the field names from the second part.
                // TODO what if indices are used?
                // Remove the closing square bracket first.
                std::string tmp_field_names = first_split_parts[1];
                AssertThrow(tmp_field_names[tmp_field_names.size()-1] == ']', ExcMessage ("While parsing the entry <Boundary composition model/Fixed "
                                                                                          "boundary composition indicators>, there was an error. The list of fields for a "
                                                                                          "specific boundary needs to be enclosed in square brackets. The closing bracket "
                                                                                          "is missing."));
                tmp_field_names.erase (--tmp_field_names.end());

                const std::vector<std::string> field_names = Utilities::split_string_list(tmp_field_names, '|');

                // Loop over the field names and set their mask to true if they exist.
                ComponentMask component_mask(this->n_compositional_fields(),
                                             false);

                for (auto &field_name : field_names)
                  {
                    AssertThrow(this->introspection().compositional_name_exists(field_name),
                                ExcMessage ("The compositional field name "
                                            + field_name
                                            + "listed in `Fixed composition boundary indicators' does not exist."));
                    component_mask.set(this->introspection().compositional_index_for_name(field_name), true);

                    fixed_compositional_fields.insert(this->introspection().compositional_index_for_name(field_name));

                  }

                unsigned int mn = 0;
                for (auto &model_name : plugins)
                  {
                    // Multiple boundaries can list the same plugins.
                    // We do not need multiple entries in plugin_names, however.
                    const auto plugin_name_it = std::find(this->plugin_names.begin(), this->plugin_names.end(), model_name);
                    if (plugin_name_it == this->plugin_names.end())
                      {
                        // The plugin was not listed before, so add it.
                        this->plugin_names.push_back(model_name);

                        // Also add the compositional field masks for this plugin.
                        std::vector<ComponentMask> vec_component_mask;
                        vec_component_mask.push_back(component_mask);
                        masks_fields.push_back(vec_component_mask);

                        // Insert the boundary indicator that belongs to this plugin into the set.
                        fixed_composition_boundary_indicators.insert(boundary_indicator);

                        // Insert the boundary indicator that belongs to this plugin into the vector.
                        // It is possible that more boundary indicators are affected by this plugin,
                        // they would be added to the vector of indicators for this plugin
                        // in the other conditional branch.
                        std::vector<types::boundary_id> vec_boundary_indicator;
                        vec_boundary_indicator.push_back(boundary_indicator);
                        boundary_indicators.push_back(vec_boundary_indicator);

                        // Add the model operator for this plugin.
                        std::vector<aspect::Utilities::Operator> vec_model_operator;
                        vec_model_operator.push_back(list_of_model_operators[mn]);
                        model_operators.push_back(vec_model_operator);
                      }
                    // The plugin was already listed for another boundary or field,
                    // so update the fields masks to all be true and add the boundary
                    // indicator at the existing entry.
                    else
                      {
                        // Make sure the existing plugin entry does not prescribe the same fields
                        // on the same boundary.
                        for (unsigned int bi = 0; bi < boundary_indicators[plugin_name_it - this->plugin_names.begin()].size(); ++ bi)
                          {
                            if (boundary_indicators[plugin_name_it - this->plugin_names.begin()][bi] == boundary_indicator)
                              {
                                for (unsigned int mf = 0; mf < masks_fields[plugin_name_it - this->plugin_names.begin()][bi].size(); ++mf)
                                  {
                                    if (masks_fields[plugin_name_it - this->plugin_names.begin()][bi][mf] == true)
                                      AssertThrow (masks_fields[plugin_name_it - this->plugin_names.begin()][bi][mf] != component_mask[mf],
                                                   ExcMessage ("The plugin " + model_name + " is listed multiple times for the same compositional field " + this->introspection().name_for_compositional_index(mf) + " on the same boundary " + dealii::Utilities::int_to_string(boundary_indicator) + "."));
                                  }
                              }
                          }

                        // Add the masks for this plugin and this boundary indicator.
                        masks_fields[plugin_name_it - this->plugin_names.begin()].push_back(component_mask);

                        // Add the boundary indicator for this plugin.
                        boundary_indicators[plugin_name_it - this->plugin_names.begin()].push_back(boundary_indicator);

                        // Also insert this boundary indicator for the existing plugin.
                        fixed_composition_boundary_indicators.insert(boundary_indicator);

                        // Add the operator for this boundary for the existing plugin.
                        model_operators[plugin_name_it - this->plugin_names.begin()].push_back(list_of_model_operators[mn]);
                      }
                    ++mn;
                  }

              }
            else
              AssertThrow (false, ExcMessage ("The format for fixed composition boundary indicators "
                                              "requires that each entry consists of either a boundary name "
                                              "or `<boundary_name> [field_name_1|field_name_2] : plugin_name_1+plugin_name_2'."
                                              "The entry "
                                              + p
                                              + "does not appear to follow this format."));


            // TODO deal with deprecated parameter
            const std::string model_name = prm.get ("Model name");

            AssertThrow (model_name == "unspecified" || this->plugin_names.size() == 0,
                         ExcMessage ("The parameter 'Model name' is only used for reasons"
                                     "of backwards compatibility and can not be used together with "
                                     "the new functionality 'List of model names'. Please add your "
                                     "boundary composition model to the list instead."));

            if (!(model_name == "unspecified"))
              this->plugin_names.push_back(model_name);


            if (prm.get ("Allow fixed composition on outflow boundaries") == "true")
              allow_fixed_composition_on_outflow_boundaries = true;
            else if (prm.get ("Allow fixed composition on outflow boundaries") == "false")
              allow_fixed_composition_on_outflow_boundaries = false;
            else if (prm.get ("Allow fixed composition on outflow boundaries") == "false for models without melt")
              allow_fixed_composition_on_outflow_boundaries = this->get_parameters().include_melt_transport;
            else
              AssertThrow(false, ExcMessage("'Allow fixed composition on outflow boundaries' "
                                            "must be set to 'true' or 'false', or to its default value."));
          }
      }
      prm.leave_subsection ();

      // Check whether some boundaries are fixed for only a subset of fields.
      // TODO no need for iterator i.
      auto fcbi = fixed_composition_boundary_indicators.begin();
      for (unsigned int i=0; i<fixed_composition_boundary_indicators.size(); ++i, ++fcbi)
        {
          if (get_fixed_fields_on_boundary(*fcbi).size() < this->n_compositional_fields())
            {
              boundaries_with_fixed_subset_of_fields = true;
              break;
            }
        }

      // go through the list, create objects and let them parse
      // their own parameters
      for (auto &model_name : this->plugin_names)
        {
          // create boundary composition objects
          this->plugin_objects.emplace_back (std::get<dim>(registered_plugins)
                                             .create_plugin (model_name,
                                                             "Boundary composition::Model names"));

          if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(this->plugin_objects.back().get()))
            sim->initialize_simulator (this->get_simulator());

          this->plugin_objects.back()->parse_parameters (prm);
          this->plugin_objects.back()->initialize ();
        }

      // Checks on sizes of the different vectors
      Assert (this->plugin_objects.size() == boundary_indicators.size(), ExcMessage ("The number of boundary composition objects does not agree with the number of boundary indicator entries for the plugins."));
      Assert (this->plugin_objects.size() == masks_fields.size(), ExcMessage ("The number of boundary composition objects does not agree with the number of field mask entries for the plugins."));
      auto p = this->plugin_objects.begin();
      for (unsigned int i=0; i<this->plugin_objects.size(); ++p, ++i)
        {
          Assert (boundary_indicators[i].size() == masks_fields[i].size(), ExcMessage ("The number of boundary indicators does not agree with the number of field mask entries for each plugin."));
        }
    }



    template <int dim>
    double
    Manager<dim>::boundary_composition (const types::boundary_id boundary_indicator,
                                        const Point<dim> &position,
                                        const unsigned int compositional_field) const
    {
      // Loop over the plugin objects and for each object,
      // check whether it applies to the given boundary and field.
      // TODO use field_is_fixed_on_boundary function
      double composition = 0.0;
      bool found_plugin = false;

      auto p = this->plugin_objects.begin();
      for (unsigned int i=0; i<this->plugin_objects.size(); ++p, ++i)
        {
          for (unsigned int bi=0; bi<boundary_indicators[i].size(); ++bi)
            {
              // If the plugin is active for the given boundary ...
              if (boundary_indicators[i][bi] == boundary_indicator)
                {
                  // check that the mask is true for the given field.
                  for (unsigned int c = 0; c<this->n_compositional_fields(); ++c)
                    {
                      if (c == compositional_field && masks_fields[i][bi][c] == true)
                        {
                          found_plugin = true;
                          composition = model_operators[i][bi](composition,
                                                               (*p)->boundary_composition(boundary_indicator,
                                                                                          position,
                                                                                          compositional_field));
                        }
                    }
                }
            }
        }

      (void) found_plugin;
      Assert(found_plugin == true,
             ExcMessage("The boundary composition manager class was asked for the "
                        "boundary composition at boundary " + dealii::Utilities::int_to_string(boundary_indicator) +
                        " which contains no active boundary composition plugin for field " +
                        dealii::Utilities::int_to_string(compositional_field) + "."));

      return composition;
    }



    template <int dim>
    const std::vector<std::string> &
    Manager<dim>::get_active_boundary_composition_names () const
    {
      return this->plugin_names;
    }



    template <int dim>
    const std::list<std::unique_ptr<Interface<dim>>> &
    Manager<dim>::get_active_boundary_composition_conditions () const
    {
      return this->plugin_objects;
    }



    template <int dim>
    const std::set<types::boundary_id> &
    Manager<dim>::get_fixed_composition_boundary_indicators() const
    {
      return fixed_composition_boundary_indicators;
    }



    template <int dim>
    bool
    Manager<dim>::allows_fixed_composition_on_outflow_boundaries() const
    {
      return allow_fixed_composition_on_outflow_boundaries;
    }



    template <int dim>
    bool
    Manager<dim>::field_is_fixed_on_boundary(const types::boundary_id boundary_id,
                                             const unsigned int compositional_field) const
    {
      Assert(fixed_composition_boundary_indicators.find(boundary_id) != fixed_composition_boundary_indicators.end(),
             ExcMessage("The boundary composition manager class was asked whether the "
                        "compositional field <"
                        +
                        this->introspection().name_for_compositional_index(compositional_field)
                        +
                        "> is fixed on boundary indicator <"
                        +
                        Utilities::int_to_string(boundary_id)
                        +
                        "> with symbolic name <"
                        +
                        this->get_geometry_model().translate_id_to_symbol_name(boundary_id)
                        +
                        ">, but this boundary is not part of the active boundary composition plugins."));

      // Since the component masks of plugins at the same boundary can vary,
      // loop over all the plugins of this boundary and see if any masks for
      // the given compositional field are true.
      bool field_set_on_boundary = false;
      auto p = this->plugin_objects.begin();
      for (unsigned int i=0; i<this->plugin_objects.size(); ++p, ++i)
        {
          for (unsigned int bi=0; bi<boundary_indicators[i].size(); ++bi)
            {
              // If the plugin is active for the given boundary ...
              if (boundary_indicators[i][bi] == boundary_id)
                {
                  // check that the mask is true for the given field.
                  for (unsigned int c = 0; c<this->n_compositional_fields(); ++c)
                    {
                      if (c == compositional_field && masks_fields[i][bi][c] == true)
                        {
                          field_set_on_boundary = true;
                        }
                    }
                }
            }
        }
      return field_set_on_boundary;
    }



    template <int dim>
    bool
    Manager<dim>::boundaries_with_fixed_subset_of_fields_exist() const
    {
      return boundaries_with_fixed_subset_of_fields;
    }



    template <int dim>
    std::vector<unsigned int>
    Manager<dim>::get_fixed_fields_on_boundary (const types::boundary_id boundary_id) const
    {
      Assert(fixed_composition_boundary_indicators.find(boundary_id) != fixed_composition_boundary_indicators.end(),
             ExcMessage("The boundary composition manager class was asked for the "
                        "compositional fields that are fixed on boundary indicator <"
                        +
                        Utilities::int_to_string(boundary_id)
                        +
                        "> with symbolic name <"
                        +
                        this->get_geometry_model().translate_id_to_symbol_name(boundary_id)
                        +
                        ">, but this boundary is not part of the active boundary composition plugins."));

      // Since the component masks of plugins at the same boundary can vary,
      // loop over all the plugins of this boundary and see if any masks for
      // the given compositional field are true.
      std::vector<unsigned int> fixed_fields;
      auto p = this->plugin_objects.begin();
      for (unsigned int i=0; i<this->plugin_objects.size(); ++p, ++i)
        {
          for (unsigned int bi=0; bi<boundary_indicators[i].size(); ++bi)
            {
              // If the plugin is active for the given boundary ...
              if (boundary_indicators[i][bi] == boundary_id)
                {
                  // check whether any masks are true.
                  for (unsigned int c = 0; c<this->n_compositional_fields(); ++c)
                    {
                      if (masks_fields[i][bi][c] == true)
                        {
                          // There can be multiple plugins active for a field on a certain boundary,
                          // so only add the field if it is not listed already.
                          if (std::find(fixed_fields.begin(), fixed_fields.end(), c) == fixed_fields.end())
                            {
                              fixed_fields.push_back(c);
                            }
                        }
                    }
                }
            }
        }

      return fixed_fields;
    }



    template <int dim>
    std::set<types::boundary_id>
    Manager<dim>::get_fixed_boundaries_for_field (const unsigned int compositional_field) const
    {
      Assert(this->introspection().compositional_name_exists(this->introspection().name_for_compositional_index(compositional_field)),
             ExcMessage("The boundary composition manager class was asked for the "
                        "boundaries on which compositional field <"
                        +
                        this->introspection().name_for_compositional_index(compositional_field)
                        +
                        "> is fixed, but this field does not exist."));

      std::set<types::boundary_id> fixed_boundaries;

      // Since the component masks of plugins at the same boundary can vary,
      // loop over all the plugins of this boundary and see if any masks for
      // the given compositional field are true.
      auto p = this->plugin_objects.begin();
      for (unsigned int i=0; i<this->plugin_objects.size(); ++p, ++i)
        {
          for (unsigned int bi=0; bi<boundary_indicators[i].size(); ++bi)
            {
              // check whether the mask is true for the given field.
              for (unsigned int c = 0; c<this->n_compositional_fields(); ++c)
                {
                  if (c == compositional_field && masks_fields[i][bi][c] == true)
                    {
                      fixed_boundaries.insert(boundary_indicators[i][bi]);
                    }
                }
            }
        }
      return fixed_boundaries;
    }


    template <int dim>
    std::set<unsigned int>
    Manager<dim>::get_fixed_compositional_fields () const
    {
      return fixed_compositional_fields;
    }



    template <int dim>
    std::set<unsigned int>
    Manager<dim>::get_fixed_compositional_fields_for_plugin (const std::string plugin_name) const
    {
// If there are no boundaries for which only a subset of fields
// is fixed, then all plugins prescribe all fields.
      if (!boundaries_with_fixed_subset_of_fields)
        return fixed_compositional_fields;
      // Loop over the plugin names and for the matching name,
      // store which fields it prescribes.
      std::set<unsigned int> fixed_compositional_fields_for_plugin;

      auto p = this->plugin_names.begin();
      for (unsigned int i=0; i<this->plugin_names.size(); ++p, ++i)
        {
          if (*p == plugin_name)
            {
              for (unsigned int bi=0; bi<boundary_indicators[i].size(); ++bi)
                {
                  // check that the mask is true for the given field.
                  for (unsigned int c = 0; c<this->n_compositional_fields(); ++c)
                    {
                      if (masks_fields[i][bi][c] == true)
                        {
                          fixed_compositional_fields_for_plugin.insert(c);
                        }
                    }
                }
            }
        }

      return fixed_compositional_fields_for_plugin;
    }


    template <int dim>
    std::set<unsigned int>
    Manager<dim>::get_fixed_compositional_fields_for_plugin_on_boundary (const std::string plugin_name, const types::boundary_id boundary_id) const
    {
// If there are no boundaries for which only a subset of fields
// is fixed, then all plugins prescribe all fields.
      if (!boundaries_with_fixed_subset_of_fields)
        return fixed_compositional_fields;

      // Loop over the plugin names and for the matching name and
      // boundary indicator, store which fields it prescribes.
      std::set<unsigned int> fixed_compositional_fields_for_plugin;

      auto p = this->plugin_names.begin();
      for (unsigned int i=0; i<this->plugin_names.size(); ++p, ++i)
        {
          if (*p == plugin_name)
            {
              for (unsigned int bi=0; bi<boundary_indicators[i].size(); ++bi)
                {
                  if (boundary_indicators[i][bi] == boundary_id)
                    {
                      // check that the mask is true for the given field.
                      for (unsigned int c = 0; c<this->n_compositional_fields(); ++c)
                        {
                          if (masks_fields[i][bi][c] == true)
                            {
                              fixed_compositional_fields_for_plugin.insert(c);
                            }
                        }
                    }
                }
            }
        }

      return fixed_compositional_fields_for_plugin;
    }



    template <int dim>
    void
    Manager<dim>::declare_parameters (ParameterHandler &prm)
    {
      // declare the entry in the parameter file
      prm.enter_subsection ("Boundary composition model");
      {
        const std::string pattern_of_names
          = std::get<dim>(registered_plugins).get_pattern_of_names ();

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
                          std::get<dim>(registered_plugins).get_description_string());

        prm.declare_entry("List of model operators", "add",
                          Patterns::MultipleSelection(Utilities::get_model_operator_options()),
                          "A comma-separated list of operators that "
                          "will be used to append the listed composition models onto "
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
                           "boundary composition models and shouldn't be used. "
                           "Please use 'List of model names' instead.");

        prm.declare_entry ("Fixed composition boundary indicators", "",
                           Patterns::List (Patterns::Anything()),
                           "A comma separated list of names denoting those boundaries "
                           "on which the composition is fixed and described by the "
                           "boundary composition object selected in its own section "
                           "of this input file or in this list. "
                           "All boundary indicators used by the geometry "
                           "but not explicitly listed here will end up with no-flux "
                           "(insulating) boundary conditions."
                           "\n\n"
                           "The names of the boundaries listed here can either be "
                           "numbers (in which case they correspond to the numerical "
                           "boundary indicators assigned by the geometry object), or they "
                           "can correspond to any of the symbolic names the geometry object "
                           "may have provided for each part of the boundary. You may want "
                           "to compare this with the documentation of the geometry model you "
                           "use in your model."
                           "\n\n"
                           "To specify different boundary conditions for different compositions "
                           "on different boundaries, "
                           "you can use the format <boundary_name> [field_name_1|field_name_2]:"
                           "plugin_name_1+plugin_name_2. This also means that not all fields need "
                           "to have fixed compositions on a specific boundary.");
        prm.declare_entry ("Allow fixed composition on outflow boundaries", "false for models without melt",
                           Patterns::Selection("true|false|false for models without melt"),
                           "When the composition is fixed on a given boundary as determined "
                           "by the list of 'Fixed composition boundary indicators', there "
                           "might be parts of the boundary where material flows out and "
                           "one may want to prescribe the composition only on those parts of "
                           "the boundary where there is inflow. This parameter determines "
                           "if compositions are only prescribed at these inflow parts of the "
                           "boundary (if false) or everywhere on a given boundary, independent "
                           "of the flow direction (if true). By default, this parameter is set "
                           "to false, except in models with melt transport (see below). "
                           "Note that in this context, `fixed' refers to the fact that these "
                           "are the boundary indicators where Dirichlet boundary conditions are "
                           "applied, and does not imply that the boundary composition is "
                           "time-independent. "
                           "\n\n"
                           "Mathematically speaking, the compositional fields satisfy an "
                           "advection equation that has no diffusion. For this equation, one "
                           "can only impose Dirichlet boundary conditions (i.e., prescribe a "
                           "fixed compositional field value at the boundary) at those boundaries "
                           "where material flows in. This would correspond to the ``false'' "
                           "setting of this parameter, which is correspondingly the default. "
                           "On the other hand, on a finite dimensional discretization such as "
                           "the one one obtains from the finite element method, it is possible "
                           "to also prescribe values on outflow boundaries, even though this may "
                           "make no physical sense. This would then correspond to the ``true'' "
                           "setting of this parameter. Note however that this parameter is only "
                           "taken into account for the continuous field method and is not "
                           "applied to the Discontinuous Galerkin (DG) field method. "
                           "\n\n"
                           "A warning for models with melt transport: In models with fluid flow, "
                           "some compositional fields (in particular the porosity) might be "
                           "transported with the fluid velocity, and would need to set the "
                           "constraints based on the fluid velocity. However, this is currently "
                           "not possible, because we reuse the same matrix for all compositional "
                           "fields, and therefore can not use different constraints for different "
                           "fields. Consequently, we set this parameter to true by default in "
                           "models where melt transport is enabled. Be aware that if you change "
                           "this default setting, you will not use the melt velocity, but the solid "
                           "velocity to determine on which parts of the boundaries there is outflow.");
      }
      prm.leave_subsection ();

      std::get<dim>(registered_plugins).declare_parameters (prm);
    }



    template <int dim>
    void
    Manager<dim>::write_plugin_graph (std::ostream &out)
    {
      std::get<dim>(registered_plugins).write_plugin_graph ("Boundary composition interface",
                                                            out);
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace BoundaryComposition
  {
#define INSTANTIATE(dim) \
  template class Interface<dim>; \
  template class Manager<dim>;

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
