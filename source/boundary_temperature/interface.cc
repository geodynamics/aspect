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
#include <aspect/boundary_temperature/interface.h>
#include <aspect/simulator_access.h>

#include <aspect/utilities.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/signaling_nan.h>

#include <list>
#include <tuple>


namespace aspect
{
  namespace BoundaryTemperature
  {
    // ------------------------------ Manager -----------------------------
    // -------------------------------- Deal with registering boundary_temperature models and automating
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
    Manager<dim>::register_boundary_temperature (const std::string &name,
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
      prm.enter_subsection ("Boundary temperature model");
      {
        this->plugin_names
          = Utilities::split_string_list(prm.get("List of model names"));

        AssertThrow(Utilities::has_unique_entries(this->plugin_names),
                    ExcMessage("The list of strings for the parameter "
                               "'Boundary temperature model/List of model names' contains entries more than once. "
                               "This is not allowed. Please check your parameter file."));

        const std::string model_name = prm.get ("Model name");

        AssertThrow (model_name == "unspecified" || this->plugin_names.size() == 0,
                     ExcMessage ("The parameter 'Model name' is only used for reasons"
                                 "of backwards compatibility and can not be used together with "
                                 "the new functionality 'List of model names'. Please add your "
                                 "boundary temperature model to the list instead."));

        if (!(model_name == "unspecified"))
          this->plugin_names.push_back(model_name);

        // create operator list
        std::vector<std::string> model_operator_names =
          Utilities::possibly_extend_from_1_to_N (Utilities::split_string_list(prm.get("List of model operators")),
                                                  this->plugin_names.size(),
                                                  "List of model operators");
        model_operators = Utilities::create_model_operator_list(model_operator_names);

        try
          {
            const std::vector<types::boundary_id> x_fixed_temperature_boundary_indicators
              = this->get_geometry_model().translate_symbolic_boundary_names_to_ids(Utilities::split_string_list
                                                                                    (prm.get ("Fixed temperature boundary indicators")));
            fixed_temperature_boundary_indicators
              = std::set<types::boundary_id> (x_fixed_temperature_boundary_indicators.begin(),
                                              x_fixed_temperature_boundary_indicators.end());

            // If no fixed temperature boundary indicators have been set, there should be no model_names chosen either.
            // If that is indeed the case, clear the model_operators vector. Otherwise, raise an exception.
            if (fixed_temperature_boundary_indicators.size() == 0)
              {
                AssertThrow(this->plugin_names.size() == 0,
                            ExcMessage ("You have indicated that you wish to apply a boundary temperature "
                                        "model, but the <Fixed temperature boundary indicators> parameter "
                                        "is empty. Please use this parameter to specify the boundaries "
                                        "on which the model(s) should be applied."));

                model_operators.clear();
              }
          }
        catch (const std::string &error)
          {
            AssertThrow (false, ExcMessage ("While parsing the entry <Model settings/Fixed temperature "
                                            "boundary indicators>, there was an error. Specifically, "
                                            "the conversion function complained as follows:\n\n"
                                            + error));
          }

        allow_fixed_temperature_on_outflow_boundaries = prm.get_bool ("Allow fixed temperature on outflow boundaries");
      }
      prm.leave_subsection ();

      // go through the list, create objects and let them parse
      // their own parameters
      for (auto &model_name : this->plugin_names)
        {
          // create boundary temperature objects
          this->plugin_objects.emplace_back (std::get<dim>(registered_plugins)
                                             .create_plugin (model_name,
                                                             "Boundary temperature::Model names"));

          if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(this->plugin_objects.back().get()))
            sim->initialize_simulator (this->get_simulator());

          this->plugin_objects.back()->parse_parameters (prm);
          this->plugin_objects.back()->initialize ();
        }
    }



    template <int dim>
    double
    Manager<dim>::boundary_temperature (const types::boundary_id boundary_indicator,
                                        const Point<dim> &position) const
    {
      double temperature = 0.0;

      auto p = this->plugin_objects.begin();
      for (unsigned int i=0; i<this->plugin_objects.size(); ++p, ++i)
        temperature = model_operators[i](temperature,
                                         (*p)->boundary_temperature(boundary_indicator,
                                                                    position));

      return temperature;
    }



    template <int dim>
    double
    Manager<dim>::minimal_temperature (const std::set<types::boundary_id> &fixed_boundary_ids) const
    {
      double temperature = std::numeric_limits<double>::max();

      for (const auto &p : this->plugin_objects)
        temperature = std::min(temperature,
                               p->minimal_temperature(fixed_boundary_ids));

      return temperature;
    }



    template <int dim>
    double
    Manager<dim>::maximal_temperature (const std::set<types::boundary_id> &fixed_boundary_ids) const
    {
      double temperature = 0.0;

      for (const auto &p : this->plugin_objects)
        temperature = std::max(temperature,
                               p->maximal_temperature(fixed_boundary_ids));

      return temperature;
    }



    template <int dim>
    const std::vector<std::string> &
    Manager<dim>::get_active_boundary_temperature_names () const
    {
      return this->plugin_names;
    }


    template <int dim>
    const std::list<std::unique_ptr<Interface<dim>>> &
    Manager<dim>::get_active_boundary_temperature_conditions () const
    {
      return this->plugin_objects;
    }



    template <int dim>
    const std::set<types::boundary_id> &
    Manager<dim>::get_fixed_temperature_boundary_indicators() const
    {
      return fixed_temperature_boundary_indicators;
    }



    template <int dim>
    bool
    Manager<dim>::allows_fixed_temperature_on_outflow_boundaries() const
    {
      return allow_fixed_temperature_on_outflow_boundaries;
    }



    template <int dim>
    void
    Manager<dim>::declare_parameters (ParameterHandler &prm)
    {
      // declare the entry in the parameter file
      prm.enter_subsection ("Boundary temperature model");
      {
        const std::string pattern_of_names
          = std::get<dim>(registered_plugins).get_pattern_of_names ();

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
                           "boundary temperature models and shouldn't be used. "
                           "Please use 'List of model names' instead.");

        prm.declare_entry ("Fixed temperature boundary indicators", "",
                           Patterns::List (Patterns::Anything()),
                           "A comma separated list of names denoting those boundaries "
                           "on which the temperature is fixed and described by the "
                           "boundary temperature object selected in the 'List of model names' "
                           "parameter. All boundary indicators used by the geometry "
                           "but not explicitly listed here will end up with no-flux "
                           "(insulating) boundary conditions, or, if they are listed in the "
                           "'Fixed heat flux boundary indicators', with Neumann boundary "
                           "conditions."
                           "\n\n"
                           "The names of the boundaries listed here can either be "
                           "numbers (in which case they correspond to the numerical "
                           "boundary indicators assigned by the geometry object), or they "
                           "can correspond to any of the symbolic names the geometry object "
                           "may have provided for each part of the boundary. You may want "
                           "to compare this with the documentation of the geometry model you "
                           "use in your model."
                           "\n\n"
                           "This parameter only describes which boundaries have a fixed "
                           "temperature, but not what temperature should hold on these "
                           "boundaries. The latter piece of information needs to be "
                           "implemented in a plugin in the BoundaryTemperature "
                           "group, unless an existing implementation in this group "
                           "already provides what you want.");
        prm.declare_entry ("Allow fixed temperature on outflow boundaries", "true",
                           Patterns::Bool (),
                           "When the temperature is fixed on a given boundary as determined "
                           "by the list of 'Fixed temperature boundary indicators', there "
                           "might be parts of the boundary where material flows out and "
                           "one may want to prescribe the temperature only on the parts of "
                           "the boundary where there is inflow. This parameter determines "
                           "if temperatures are only prescribed at these inflow parts of the "
                           "boundary (if false) or everywhere on a given boundary, independent "
                           "of the flow direction (if true)."
                           "Note that in this context, `fixed' refers to the fact that these "
                           "are the boundary indicators where Dirichlet boundary conditions are "
                           "applied, and does not imply that the boundary temperature is "
                           "time-independent. "
                           "\n\n"
                           "Mathematically speaking, the temperature satisfies an "
                           "advection-diffusion equation. For this type of equation, one can "
                           "prescribe the temperature even on outflow boundaries as long as the "
                           "diffusion coefficient is nonzero. This would correspond to the "
                           "``true'' setting of this parameter, which is correspondingly the "
                           "default. In practice, however, this would only make physical sense "
                           "if the diffusion coefficient is actually quite large to prevent "
                           "the creation of a boundary layer. "
                           "In addition, if there is no diffusion, one can only impose "
                           "Dirichlet boundary conditions (i.e., prescribe a fixed temperature "
                           "value at the boundary) at those boundaries where material flows in. "
                           "This would correspond to the ``false'' setting of this parameter.");
      }
      prm.leave_subsection ();

      std::get<dim>(registered_plugins).declare_parameters (prm);
    }



    template <int dim>
    void
    Manager<dim>::write_plugin_graph (std::ostream &out)
    {
      std::get<dim>(registered_plugins).write_plugin_graph ("Boundary temperature interface",
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
      std::list<internal::Plugins::PluginList<BoundaryTemperature::Interface<2>>::PluginInfo> *
      internal::Plugins::PluginList<BoundaryTemperature::Interface<2>>::plugins = nullptr;
      template <>
      std::list<internal::Plugins::PluginList<BoundaryTemperature::Interface<3>>::PluginInfo> *
      internal::Plugins::PluginList<BoundaryTemperature::Interface<3>>::plugins = nullptr;
    }
  }

  namespace BoundaryTemperature
  {
#define INSTANTIATE(dim) \
  template class Interface<dim>; \
  template class Manager<dim>;

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
