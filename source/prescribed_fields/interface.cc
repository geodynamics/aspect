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

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>
#include <deal.II/fe/fe_values.h>
#include <aspect/global.h>
#include <aspect/utilities.h>
#include <aspect/simulator_signals.h>
#include <aspect/parameters.h>
#include <aspect/prescribed_fields/interface.h>

namespace aspect
{
  namespace PrescribedFields
  {

    /**
    * A set of helper functions that either return the point passed to it (if
    * the current dimension is the same) or return a dummy value (otherwise).
    */
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
    register_prescribed_fields_model (const std::string &name,
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
    parse_parameters (const Parameters<dim> &, ParameterHandler &prm)
    {
      prescribe_internal_fields = prm.get_bool ("Prescribe fields");
      prm.enter_subsection ("Prescribed fields");
      {
        plugin_names
          = Utilities::split_string_list(prm.get("List of model names"));

        AssertThrow(Utilities::has_unique_entries(plugin_names),
                    ExcMessage("The list of strings for the parameter "
                               "'Initial temperature model/List of model names' contains entries more than once. "
                               "This is not allowed. Please check your parameter file."));
      }
      prm.leave_subsection ();

      // go through the list, create objects and let them parse
      // their own parameters
      if (dim == 2)
        for (auto &model_name : plugin_names)
          {
            // create initial temperature objects
            plugin_objects_2d.emplace_back (std::get<2>(registered_plugins)
                                            .create_plugin (model_name,
                                                            "Prescribed field model::Model names"));

            plugin_objects_2d.back()->parse_parameters (prm);
            plugin_objects_2d.back()->initialize ();
          }
      else
        for (auto &model_name : plugin_names)
          {
            // create initial temperature objects
            plugin_objects_3d.emplace_back (std::get<3>(registered_plugins)
                                            .create_plugin (model_name,
                                                            "Prescribed field model::Model names"));

            plugin_objects_3d.back()->parse_parameters (prm);
            plugin_objects_3d.back()->initialize ();
          }
    }

    template <int dim>
    void
    declare_parameters (ParameterHandler &prm)
    {
      prm.declare_entry ("Prescribe fields", "false",
                         Patterns::Bool (),
                         "Whether or not to use any Prescribed internal fields. "
                         "Locations in which to prescribe velocities are defined "
                         "in section ``Prescribed fields/Indicator function'' "
                         "and the velocities are defined in section ``Prescribed "
                         "velocities/Velocity function''. Indicators are evaluated "
                         "at the center of each cell, and all DOFs associated with "
                         "the specified velocity component at the indicated cells "
                         "are constrained."
                        );

      prm.enter_subsection ("Prescribed fields");
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
      }
      prm.leave_subsection ();

      std::get<dim>(registered_plugins).declare_parameters (prm);
    }

    void
    declare_parameters_signal (const unsigned int dim, ParameterHandler &prm)
    {
      if (dim == 2)
        declare_parameters<2>(prm);
      else
        declare_parameters<3>(prm);
    }

    template <int dim>
    std::string
    get_valid_model_names_pattern ()
    {
      return std::get<dim>(registered_plugins).get_pattern_of_names ();
    }

    template <int dim>
    void
    write_plugin_graph (std::ostream &out)
    {
      std::get<dim>(registered_plugins).write_plugin_graph ("Prescribed fields interface",
                                                            out);
    }

    template <>
    void constrain_all_prescribed_internal_fields (const SimulatorAccess<2> &simulator_access,
                                                   AffineConstraints<double> &current_constraints)
    {
      for (auto &p: plugin_objects_2d)
        {
          p->update(simulator_access);
          p->constrain_internal_fields(simulator_access, current_constraints);
        }
    }

    template <>
    void constrain_all_prescribed_internal_fields (const SimulatorAccess<3> &simulator_access,
                                                   AffineConstraints<double> &current_constraints)
    {
      for (auto &p: plugin_objects_3d)
        {
          p->update(simulator_access);
          p->constrain_internal_fields(simulator_access, current_constraints);
        }
    }

    // Connect declare_parameters and parse_parameters to appropriate signals.
    void parameter_connector ()
    {
      SimulatorSignals<2>::declare_additional_parameters.connect (&declare_parameters_signal);
      SimulatorSignals<3>::declare_additional_parameters.connect (&declare_parameters_signal);

      SimulatorSignals<2>::parse_additional_parameters.connect (&parse_parameters<2>);
      SimulatorSignals<3>::parse_additional_parameters.connect (&parse_parameters<3>);
    }

    // Connect constraints function to correct signal.
    template <int dim>
    void signal_connector (SimulatorSignals<dim> &signals)
    {
      signals.post_constraints_creation.connect (&constrain_all_prescribed_internal_fields<dim>);
    }

    // Tell ASPECT to send signals to the connector functions
    ASPECT_REGISTER_SIGNALS_PARAMETER_CONNECTOR(parameter_connector)
    ASPECT_REGISTER_SIGNALS_CONNECTOR(signal_connector<2>, signal_connector<3>)
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
      std::list<internal::Plugins::PluginList<PrescribedFields::Interface<2>>::PluginInfo> *
      internal::Plugins::PluginList<PrescribedFields::Interface<2>>::plugins = nullptr;
      template <>
      std::list<internal::Plugins::PluginList<PrescribedFields::Interface<3>>::PluginInfo> *
      internal::Plugins::PluginList<PrescribedFields::Interface<3>>::plugins = nullptr;
    }
  }

  namespace PrescribedFields
  {
#define INSTANTIATE(dim) \
  template class Interface<dim>;\
  \
  template \
  void \
  register_prescribed_fields_model<dim> (const std::string &, \
                                         const std::string &, \
                                         void ( *) (ParameterHandler &), \
                                         std::unique_ptr<Interface<dim>>( *) ()); \
  \
  template  \
  void \
  declare_parameters<dim> (ParameterHandler &); \
  \
  template \
  void \
  write_plugin_graph<dim> (std::ostream &);

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }

}
