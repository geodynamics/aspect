/*
  Copyright (C) 2025 - by the authors of the ASPECT code.

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
#include <aspect/boundary_convective_heating/interface.h>

#include <aspect/utilities.h>

#include <deal.II/base/signaling_nan.h>

#include <list>
#include <tuple>


namespace aspect
{
  namespace BoundaryConvectiveHeating
  {
    // ------------------------------ Manager -----------------------------
    // -------------------------------- Deal with registering boundary_convective_heating models and automating
    // -------------------------------- their setup and selection at run time

    namespace
    {
      std::tuple
      <aspect::internal::Plugins::UnusablePluginList,
      aspect::internal::Plugins::UnusablePluginList,
      aspect::internal::Plugins::PluginList<Interface<2>>,
      aspect::internal::Plugins::PluginList<Interface<3>>> registered_plugins;

      std::tuple
      <aspect::internal::Plugins::UnusablePluginList,
      aspect::internal::Plugins::UnusablePluginList,
      aspect::internal::Plugins::PluginList<BoundaryTemperature::Interface<2>>,
      aspect::internal::Plugins::PluginList<BoundaryTemperature::Interface<3>>> registered_temperature_plugins;

      std::tuple
      <aspect::internal::Plugins::UnusablePluginList,
      aspect::internal::Plugins::UnusablePluginList,
      aspect::internal::Plugins::PluginList<BoundaryHeatFlux::Interface<2>>,
      aspect::internal::Plugins::PluginList<BoundaryHeatFlux::Interface<3>>> registered_heat_flux_plugins;
    }


    template <int dim>
    void
    Manager<dim>::register_boundary_convective_heating (const std::string &name,
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
    Manager<dim>::register_boundary_temperature (const std::string &name,
                                                 const std::string &description,
                                                 void (*declare_parameters_function) (ParameterHandler &),
                                                 std::unique_ptr<BoundaryTemperature::Interface<dim>> (*factory_function) ())
    {
      std::get<dim>(registered_temperature_plugins).register_plugin (name,
                                                                     description,
                                                                     declare_parameters_function,
                                                                     factory_function);
    }


    template <int dim>
    void
    Manager<dim>::register_boundary_heat_flux (const std::string &name,
                                               const std::string &description,
                                               void (*declare_parameters_function) (ParameterHandler &),
                                               std::unique_ptr<BoundaryHeatFlux::Interface<dim>> (*factory_function) ())
    {
      std::get<dim>(registered_heat_flux_plugins).register_plugin (name,
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
      prm.enter_subsection ("Boundary convective heating model");
      {
        this->plugin_names
          = Utilities::split_string_list(prm.get("List of heat transfer coefficient model names"));
        temperature_plugin_names
          = Utilities::split_string_list(prm.get("List of boundary temperature model names"));
        heat_flux_plugin_names
          = Utilities::split_string_list(prm.get("List of boundary heat flux model names"));

        AssertThrow(this->plugin_names.size() <= 1
                    && temperature_plugin_names.size() <= 1
                    && heat_flux_plugin_names.size() <= 1,
                    ExcMessage("The list of strings for one of the lists of model names "
                               "in the 'Boundary convective heating model' subsection "
                               "contains more than one entry. "
                               "This is not allowed. Please check your parameter file."));

        AssertThrow(Utilities::has_unique_entries(this->plugin_names)
                    && Utilities::has_unique_entries(temperature_plugin_names)
                    && Utilities::has_unique_entries(heat_flux_plugin_names),
                    ExcMessage("The list of strings for one of the lists of model names "
                               "in the 'Boundary convective heating model' subsection "
                               "contains entries more than once. "
                               "This is not allowed. Please check your parameter file."));

        try
          {
            const std::vector<types::boundary_id> x_convective_heating_boundary_indicators
              = this->get_geometry_model().translate_symbolic_boundary_names_to_ids(Utilities::split_string_list
                                                                                    (prm.get ("Convective heating boundary indicators")));
            convective_heating_boundary_indicators
              = std::set<types::boundary_id> (x_convective_heating_boundary_indicators.begin(),
                                              x_convective_heating_boundary_indicators.end());

            // If no fixed convective heating boundary indicators have been set, there should be no model_names chosen either.
            // If that not the case, raise an exception.
            if (convective_heating_boundary_indicators.size() == 0)
              {
                AssertThrow(this->plugin_names.size() == 0,
                            ExcMessage ("You have indicated that you wish to apply a boundary convective heating "
                                        "model, but the <Convective heating boundary indicators> parameter "
                                        "is empty. Please use this parameter to specify the boundaries "
                                        "on which the model(s) should be applied."));
              }
          }
        catch (const std::string &error)
          {
            AssertThrow (false, ExcMessage ("While parsing the entry <Boundary convective heating model/"
                                            "Convective heating boundary indicators>, there was an error. Specifically, "
                                            "the conversion function complained as follows:\n\n"
                                            + error));
          }
      }
      prm.leave_subsection ();

      // go through the list of boundary convective heating plugin names, create objects and let them parse
      // their own parameters
      for (auto &model_name : this->plugin_names)
        {
          // create boundary convective heating objects
          this->plugin_objects.emplace_back (std::get<dim>(registered_plugins)
                                             .create_plugin (model_name,
                                                             "Boundary convective heating model::Model names"));

          if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(this->plugin_objects.back().get()))
            sim->initialize_simulator (this->get_simulator());

          this->plugin_objects.back()->parse_parameters (prm);
          this->plugin_objects.back()->initialize ();
        }

      // go through the list of boundary temperature plugin names, create objects and let them parse
      // their own parameters
      for (auto &model_name : temperature_plugin_names)
        {
          // create boundary temperature objects
          temperature_plugin_objects.emplace_back (std::get<dim>(registered_temperature_plugins)
                                                   .create_plugin (model_name,
                                                                   "Boundary temperature::Model names"));

          if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(temperature_plugin_objects.back().get()))
            sim->initialize_simulator (this->get_simulator());

          temperature_plugin_objects.back()->parse_parameters (prm);
          temperature_plugin_objects.back()->initialize ();
        }

      // go through the list of boundary heat flux plugin names, create objects and let them parse
      // their own parameters
      for (auto &model_name : heat_flux_plugin_names)
        {
          // create boundary heat flux objects
          heat_flux_plugin_objects.emplace_back (std::get<dim>(registered_heat_flux_plugins)
                                                 .create_plugin (model_name,
                                                                 "Boundary heat flux model::Model names"));

          if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(heat_flux_plugin_objects.back().get()))
            sim->initialize_simulator (this->get_simulator());

          heat_flux_plugin_objects.back()->parse_parameters (prm);
          heat_flux_plugin_objects.back()->initialize ();
        }
    }



    template <int dim>
    std::vector<double>
    Manager<dim>::heat_transfer_coefficient (const types::boundary_id boundary_indicator,
                                             const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
                                             const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs) const
    {
      // This manager class is written to accommodate multiple plugins at least
      // in principle. However, at the moment it only works for exactly
      // one plugin because it does not yet know how to combine results
      // from multiple plugins. So assert that this is the case:
      Assert (this->plugin_objects.size() == 1, ExcNotImplemented());

      std::vector<double> coefficients;

      for (const auto &p: this->plugin_objects)
        coefficients = p->heat_transfer_coefficient(boundary_indicator,
                                                    material_model_inputs,
                                                    material_model_outputs);
      return coefficients;
    }



    template <int dim>
    double
    Manager<dim>::boundary_temperature (const types::boundary_id boundary_indicator,
                                        const Point<dim> &position) const
    {
      // This manager class is written to accommodate multiple plugins at least
      // in principle. However, at the moment it only works for exactly
      // one plugin because it does not yet know how to combine results
      // from multiple plugins. So assert that this is the case:
      Assert (temperature_plugin_objects.size() == 1, ExcNotImplemented());

      double temperature = 0.0;

      for (const auto &p: temperature_plugin_objects)
        temperature = p->boundary_temperature(boundary_indicator,
                                              position);

      return temperature;
    }



    template <int dim>
    std::vector<Tensor<1,dim>>
    Manager<dim>::heat_flux (const types::boundary_id boundary_indicator,
                             const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
                             const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
                             const std::vector<Tensor<1,dim>> &normal_vectors) const
    {
      // This manager class is written to accommodate multiple plugins at least
      // in principle. However, at the moment it only works for exactly
      // one plugin because it does not yet know how to combine results
      // from multiple plugins. So assert that this is the case:
      Assert (heat_flux_plugin_objects.size() == 1, ExcNotImplemented());

      std::vector<Tensor<1,dim>> heat_fluxes;

      for (const auto &p: heat_flux_plugin_objects)
        heat_fluxes = p->heat_flux(boundary_indicator,
                                   material_model_inputs,
                                   material_model_outputs,
                                   normal_vectors);
      return heat_fluxes;
    }



    template <int dim>
    const std::set<types::boundary_id> &
    Manager<dim>::get_convective_heating_boundary_indicators() const
    {
      return convective_heating_boundary_indicators;
    }



    template <int dim>
    void
    Manager<dim>::declare_parameters (ParameterHandler &prm)
    {
      // declare the entry in the parameter file
      prm.enter_subsection ("Boundary convective heating model");
      {
        const std::string pattern_of_names
          = std::get<dim>(registered_plugins).get_pattern_of_names ();
        const std::string pattern_of_temperature_names
          = std::get<dim>(registered_temperature_plugins).get_pattern_of_names ();
        const std::string pattern_of_heat_flux_names
          = std::get<dim>(registered_heat_flux_plugins).get_pattern_of_names ();

        prm.declare_entry("List of heat transfer coefficient model names",
                          "",
                          Patterns::MultipleSelection(pattern_of_names),
                          "A comma-separated list of boundary convective heating models that "
                          "will be used to determine the heat transfer coefficient across the boundary. "
                          "The heat transfer coefficient characterises the heat exchange "
                          "between the solid model interior and an adjacent fluid. "
                          "In the context of a Robin boundary condition, the heat transfer "
                          "coefficient governs the strength of the convective coupling: "
                          "For heat transfer coefficient --> zero, the boundary approaches "
                          "insulating (Neumann) behaviour; "
                          "For heat transfer coefficient --> infinity, the boundary approaches "
                          "a prescribed-temperature (Dirichlet) condition."
                          "The unit of the heat transfer coefficient is \\si{\\watt\\per\\meter\\squared\\per\\kelvin}."
                          "At the moment, this list can only have one entry. \n\n"
                          "The following heat transfer coefficient models are available:\n\n"
                          +
                          std::get<dim>(registered_plugins).get_description_string());
        prm.declare_entry("List of boundary temperature model names",
                          "",
                          Patterns::MultipleSelection(pattern_of_temperature_names),
                          "A comma-separated list of boundary temperature models that "
                          "will be used to determine the temperature boundary conditions. "
                          "At the moment, this list can only have one entry. \n\n"
                          "The following boundary temperature models are available:\n\n"
                          +
                          std::get<dim>(registered_temperature_plugins).get_description_string());
        prm.declare_entry("List of boundary heat flux model names",
                          "",
                          Patterns::MultipleSelection(pattern_of_heat_flux_names),
                          "A comma-separated list of boundary heat flux models that "
                          "will be used to determine the temperature boundary conditions. "
                          "At the moment, this list can only have one entry. \n\n"
                          "The following boundary heat flux models are available:\n\n"
                          +
                          std::get<dim>(registered_heat_flux_plugins).get_description_string());

        prm.declare_entry ("Convective heating boundary indicators", "",
                           Patterns::List (Patterns::Anything()),
                           "A comma separated list of names denoting those boundaries "
                           "on which convective heating should be applied, i.e., where we "
                           "want to use Robin boundary conditions. On these boundaries, we "
                           "can prescribe both a heat flux and a boundary temperature, and "
                           "the heat transfer coefficient input parameter determines the "
                           "weighting of the two conditions. The boundary temperature is "
                           "described by the boundary temperature object selected in the "
                           "'List of boundary temperature model names' parameter, and the "
                           "boundary heat flux is prescribed by the boundary heat flux object "
                           "selected in the 'List of boundary heat flux model names' parameter. "
                           "All boundary indicators used by the geometry "
                           "but not explicitly listed here will end up with no-flux "
                           "(insulating) boundary conditions, or, if they are listed in the "
                           "'Fixed heat flux boundary indicators', with Neumann boundary "
                           "conditions, or, if they are listed in the "
                           "'Fixed temperature boundary indicators', with Dirichlet boundary "
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
                           "convective flux (i.e. Robin boundarry condictions), but not what "
                           "conditions should hold on these boundaries. "
                           "The latter piece of information needs to be "
                           "implemented in a plugin in the BoundaryConvectiveHeating "
                           "group, unless an existing implementation in this group "
                           "already provides what you want.");
      }
      prm.leave_subsection ();

      std::get<dim>(registered_plugins).declare_parameters (prm);
    }



    template <int dim>
    void
    Manager<dim>::write_plugin_graph (std::ostream &out)
    {
      std::get<dim>(registered_plugins).write_plugin_graph ("Boundary convective heating interface",
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
      std::list<internal::Plugins::PluginList<BoundaryConvectiveHeating::Interface<2>>::PluginInfo> *
      internal::Plugins::PluginList<BoundaryConvectiveHeating::Interface<2>>::plugins = nullptr;
      template <>
      std::list<internal::Plugins::PluginList<BoundaryConvectiveHeating::Interface<3>>::PluginInfo> *
      internal::Plugins::PluginList<BoundaryConvectiveHeating::Interface<3>>::plugins = nullptr;
    }
  }

  namespace BoundaryConvectiveHeating
  {
#define INSTANTIATE(dim) \
  template class Interface<dim>; \
  template class Manager<dim>;

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
