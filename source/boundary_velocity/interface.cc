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
#include <aspect/boundary_velocity/interface.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/exceptions.h>
#include <tuple>

#include <list>


namespace aspect
{
  namespace BoundaryVelocity
  {
    // ------------------------------ Manager -----------------------------
    // -------------------------------- Deal with registering boundary_velocity models and automating
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
    Manager<dim>::register_boundary_velocity (const std::string &name,
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
    Tensor<1,dim>
    Manager<dim>::boundary_velocity (const types::boundary_id boundary_indicator,
                                     const Point<dim> &position) const
    {
      Tensor<1,dim> velocity;

      bool found_plugin = false;
      unsigned int i=0;
      for (const auto &plugin: this->plugin_objects)
        {
          if (boundary_indicators[i] == boundary_indicator)
            {
              found_plugin = true;
              const Tensor<1,dim> plugin_velocity = plugin->boundary_velocity(boundary_indicator,
                                                                              position);
              for (unsigned int d=0; d<dim; ++d)
                if (component_masks[i][d] == true)
                  velocity[d] += plugin_velocity[d];
            }

          ++i;
        }

      (void) found_plugin;
      Assert(found_plugin == true,
             ExcMessage("The boundary velocity manager class was asked for the "
                        "boundary velocity at a boundary that contains no active "
                        "boundary velocity plugin."));

      return velocity;
    }



    template <int dim>
    DEAL_II_DEPRECATED
    const std::map<types::boundary_id, std::pair<std::string,std::vector<std::string>>> &
    Manager<dim>::get_active_boundary_velocity_names () const
    {
      return boundary_velocity_indicators;
    }



    template <int dim>
    DEAL_II_DEPRECATED
    const std::map<types::boundary_id,std::vector<std::unique_ptr<BoundaryVelocity::Interface<dim>>>> &
    Manager<dim>::get_active_boundary_velocity_conditions () const
    {
      AssertThrow(false, ExcMessage("This function has been removed. Use the function "
                                    "get_active_plugins() of the base class ManagerBase "
                                    "instead."));
      return boundary_velocity_objects;
    }



    template <int dim>
    const std::vector<types::boundary_id> &
    Manager<dim>::get_active_plugin_boundary_indicators() const
    {
      return boundary_indicators;
    }



    template <int dim>
    ComponentMask
    Manager<dim>::get_component_mask(const types::boundary_id boundary_id) const
    {
      Assert(prescribed_velocity_boundary_indicators.find(boundary_id) != prescribed_velocity_boundary_indicators.end(),
             ExcMessage("The boundary velocity manager class was asked for the "
                        "component mask of boundary indicator <"
                        +
                        Utilities::int_to_string(boundary_id)
                        +
                        "> with symbolic name <"
                        +
                        this->get_geometry_model().translate_id_to_symbol_name(boundary_id)
                        +
                        ">, but this boundary is not part of the active boundary velocity plugins."));

      // Since all component masks of plugins at the same boundary are identical, we can use
      // the component mask of the first plugin we find that is responsible for this boundary.
      for (unsigned int i=0; i<boundary_indicators.size(); ++i)
        if (boundary_indicators[i] == boundary_id)
          return component_masks[i];

      // We should never get here if plugins and boundary indicators were set up correctly.
      AssertThrow(false, ExcInternalError());
      return ComponentMask();
    }



    template <int dim>
    const std::set<types::boundary_id> &
    Manager<dim>::get_prescribed_boundary_velocity_indicators () const
    {
      return prescribed_velocity_boundary_indicators;
    }



    template <int dim>
    const std::set<types::boundary_id> &
    Manager<dim>::get_zero_boundary_velocity_indicators () const
    {
      return zero_velocity_boundary_indicators;
    }



    template <int dim>
    const std::set<types::boundary_id> &
    Manager<dim>::get_tangential_boundary_velocity_indicators () const
    {
      return tangential_velocity_boundary_indicators;
    }



    template <int dim>
    void
    Manager<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Boundary velocity model");
      {
        prm.declare_entry ("Prescribed velocity boundary indicators", "",
                           Patterns::Map (Patterns::Anything(),
                                          Patterns::Selection(std::get<dim>(registered_plugins).get_pattern_of_names ())),
                           "A comma separated list denoting those boundaries "
                           "on which the velocity is prescribed, i.e., where unknown "
                           "external forces act to prescribe a particular velocity. This is "
                           "often used to prescribe a velocity that equals that of "
                           "overlying plates."
                           "\n\n"
                           "The format of valid entries for this parameter is that of a map "
                           "given as ``key1 [selector]: value1, key2 [selector]: value2, key3: value3, ...'' where "
                           "each key must be a valid boundary indicator (which is either an "
                           "integer or the symbolic name the geometry model in use may have "
                           "provided for this part of the boundary) "
                           "and each value must be one of the currently implemented boundary "
                           "velocity models. ``selector'' is an optional string given as a subset "
                           "of the letters `xyz' that allows you to apply the boundary conditions "
                           "only to the components listed. As an example, '1 y: function' applies "
                           "the type `function' to the y component on boundary 1. Without a selector "
                           "it will affect all components of the velocity."
                           "\n\n"
                           "Note that the no-slip boundary condition is "
                           "a special case of the current one where the prescribed velocity "
                           "happens to be zero. It can thus be implemented by indicating that "
                           "a particular boundary is part of the ones selected "
                           "using the current parameter and using ``zero velocity'' as "
                           "the boundary values. Alternatively, you can simply list the "
                           "part of the boundary on which the velocity is to be zero with "
                           "the parameter ``Zero velocity boundary indicator'' in the "
                           "current parameter section."
                           "\n\n"
                           "Note that when ``Use years in output instead of seconds'' is set "
                           "to true, velocity should be given in m/yr. "
                           "The following boundary velocity models are available:\n\n"
                           +
                           std::get<dim>(registered_plugins).get_description_string());
        prm.declare_entry ("Zero velocity boundary indicators", "",
                           Patterns::List (Patterns::Anything()),
                           "A comma separated list of names denoting those boundaries "
                           "on which the velocity is zero."
                           "\n\n"
                           "The names of the boundaries listed here can either by "
                           "numbers (in which case they correspond to the numerical "
                           "boundary indicators assigned by the geometry object), or they "
                           "can correspond to any of the symbolic names the geometry object "
                           "may have provided for each part of the boundary. You may want "
                           "to compare this with the documentation of the geometry model you "
                           "use in your model.");
        prm.declare_entry ("Tangential velocity boundary indicators", "",
                           Patterns::List (Patterns::Anything()),
                           "A comma separated list of names denoting those boundaries "
                           "on which the velocity is tangential and unrestrained, i.e., free-slip where "
                           "no external forces act to prescribe a particular tangential "
                           "velocity (although there is a force that requires the flow to "
                           "be tangential)."
                           "\n\n"
                           "The names of the boundaries listed here can either by "
                           "numbers (in which case they correspond to the numerical "
                           "boundary indicators assigned by the geometry object), or they "
                           "can correspond to any of the symbolic names the geometry object "
                           "may have provided for each part of the boundary. You may want "
                           "to compare this with the documentation of the geometry model you "
                           "use in your model.");
      }
      prm.leave_subsection ();

      std::get<dim>(registered_plugins).declare_parameters (prm);
    }



    template <int dim>
    void
    Manager<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Boundary velocity model");
      {
        try
          {
            const std::vector<types::boundary_id> x_zero_velocity_boundary_indicators
              = this->get_geometry_model().translate_symbolic_boundary_names_to_ids(Utilities::split_string_list
                                                                                    (prm.get ("Zero velocity boundary indicators")));
            zero_velocity_boundary_indicators
              = std::set<types::boundary_id> (x_zero_velocity_boundary_indicators.begin(),
                                              x_zero_velocity_boundary_indicators.end());
          }
        catch (const std::string &error)
          {
            AssertThrow (false, ExcMessage ("While parsing the entry <Boundary velocity model/Zero velocity "
                                            "boundary indicators>, there was an error. Specifically, "
                                            "the conversion function complained as follows:\n\n"
                                            + error));
          }

        try
          {
            const std::vector<types::boundary_id> x_tangential_velocity_boundary_indicators
              = this->get_geometry_model().translate_symbolic_boundary_names_to_ids(Utilities::split_string_list
                                                                                    (prm.get ("Tangential velocity boundary indicators")));
            tangential_velocity_boundary_indicators
              = std::set<types::boundary_id> (x_tangential_velocity_boundary_indicators.begin(),
                                              x_tangential_velocity_boundary_indicators.end());
          }
        catch (const std::string &error)
          {
            AssertThrow (false, ExcMessage ("While parsing the entry <Boundary velocity model/Tangential velocity "
                                            "boundary indicators>, there was an error. Specifically, "
                                            "the conversion function complained as follows:\n\n"
                                            + error));
          }

        // find out which plugins are requested and the various other
        // parameters we declare here
        const std::vector<std::string> x_boundary_velocity_indicators
          = Utilities::split_string_list(prm.get("Prescribed velocity boundary indicators"));

        for (const auto &p : x_boundary_velocity_indicators)
          {
            // each entry has the format (white space is optional):
            // <id> [x][y][z] : <value (might have spaces)>
            //
            // first tease apart the two halves
            const std::vector<std::string> split_parts = Utilities::split_string_list (p, ':');
            AssertThrow (split_parts.size() == 2,
                         ExcMessage ("The format for prescribed velocity boundary indicators "
                                     "requires that each entry has the form `"
                                     "<id> [x][y][z] : <value>', but there does not "
                                     "appear to be a colon in the entry <"
                                     + p
                                     + ">."));

            // the easy part: get the value
            const std::string &value = split_parts[1];

            // now for the rest. since we don't know whether there is a
            // component selector, start reading at the end and subtracting
            // letters x, y and z
            std::string key_and_comp = split_parts[0];
            std::string comp;
            while ((key_and_comp.size()>0) &&
                   ((key_and_comp[key_and_comp.size()-1] == 'x')
                    ||
                    (key_and_comp[key_and_comp.size()-1] == 'y')
                    ||
                    ((key_and_comp[key_and_comp.size()-1] == 'z') && (dim==3))))
              {
                comp += key_and_comp[key_and_comp.size()-1];
                key_and_comp.erase (--key_and_comp.end());
              }

            // we've stopped reading component selectors now. there are three
            // possibilities:
            // - no characters are left. this means that key_and_comp only
            //   consisted of a single word that only consisted of 'x', 'y'
            //   and 'z's. then this would have been a mistake to classify
            //   as a component selector, and we better undo it
            // - the last character of key_and_comp is not a whitespace. this
            //   means that the last word in key_and_comp ended in an 'x', 'y'
            //   or 'z', but this was not meant to be a component selector.
            //   in that case, put these characters back.
            // - otherwise, we split successfully. eat spaces that may be at
            //   the end of key_and_comp to get key
            if (key_and_comp.size() == 0)
              key_and_comp.swap (comp);
            else if (key_and_comp[key_and_comp.size()-1] != ' ')
              {
                key_and_comp += comp;
                comp = "";
              }
            else
              {
                while ((key_and_comp.size()>0) && (key_and_comp[key_and_comp.size()-1] == ' '))
                  key_and_comp.erase (--key_and_comp.end());
              }

            // finally, try to translate the key into a boundary_id. then
            // make sure we haven't seen it yet
            types::boundary_id boundary_id;
            try
              {
                boundary_id = this->get_geometry_model().translate_symbolic_boundary_name_to_id(key_and_comp);
              }
            catch (const std::string &error)
              {
                AssertThrow (false, ExcMessage ("While parsing the entry <Boundary velocity model/Prescribed "
                                                "velocity indicators>, there was an error. Specifically, "
                                                "the conversion function complained as follows:\n\n"
                                                + error));
              }

            if (boundary_velocity_indicators.find(boundary_id) != boundary_velocity_indicators.end())
              {
                Assert(boundary_velocity_indicators[boundary_id].first == comp,
                       ExcMessage("Different velocity plugins for the same boundary have to have the same component selector. "
                                  "This was not the case for boundary: " + key_and_comp +
                                  ", for plugin: " + value + ", with component selector: " + comp));

                // finally, put it into the list
                boundary_velocity_indicators[boundary_id].second.push_back(value);
              }
            else
              {
                boundary_velocity_indicators[boundary_id] = std::make_pair(comp,std::vector<std::string>(1,value));
              }

            this->plugin_names.push_back(value);
            boundary_indicators.push_back(boundary_id);
            prescribed_velocity_boundary_indicators.insert(boundary_id);

            ComponentMask component_mask(this->introspection().n_components,
                                         false);

            if (comp.empty() || comp.find('x') != std::string::npos)
              component_mask.set(this->introspection().component_indices.velocities[0],true);
            if (comp.empty() || comp.find('y') != std::string::npos)
              component_mask.set(this->introspection().component_indices.velocities[1],true);
            if (dim == 3 && (comp.empty() || comp.find('z') != std::string::npos))
              component_mask.set(this->introspection().component_indices.velocities[2],true);

            component_masks.push_back(component_mask);
          }
      }
      prm.leave_subsection();

      // go through the list, create objects and let them parse
      // their own parameters
      for (const auto &plugin_name: this->plugin_names)
        {
          // create boundary traction objects
          this->plugin_objects.push_back(std::get<dim>(registered_plugins)
                                         .create_plugin (plugin_name,
                                                         "Boundary velocity::Model names"));

          if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(this->plugin_objects.back().get()))
            sim->initialize_simulator (this->get_simulator());

          this->plugin_objects.back()->parse_parameters (prm);
          this->plugin_objects.back()->initialize ();
        }
    }



    template <int dim>
    void
    Manager<dim>::write_plugin_graph (std::ostream &out)
    {
      std::get<dim>(registered_plugins).write_plugin_graph ("Boundary velocity interface",
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
      std::list<internal::Plugins::PluginList<BoundaryVelocity::Interface<2>>::PluginInfo> *
      internal::Plugins::PluginList<BoundaryVelocity::Interface<2>>::plugins = nullptr;
      template <>
      std::list<internal::Plugins::PluginList<BoundaryVelocity::Interface<3>>::PluginInfo> *
      internal::Plugins::PluginList<BoundaryVelocity::Interface<3>>::plugins = nullptr;
    }
  }

  namespace BoundaryVelocity
  {
#define INSTANTIATE(dim) \
  template class Interface<dim>; \
  template class Manager<dim>;

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
