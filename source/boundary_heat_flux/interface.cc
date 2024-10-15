/*
  Copyright (C) 2015 - 2024 by the authors of the ASPECT code.

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
#include <aspect/boundary_heat_flux/interface.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/exceptions.h>
#include <tuple>

#include <list>


namespace aspect
{
  namespace BoundaryHeatFlux
  {
// -------------------------------- Deal with registering models and automating
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
    register_boundary_heat_flux (const std::string &name,
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
    std::unique_ptr<Interface<dim>>
    create_boundary_heat_flux (ParameterHandler &prm)
    {
      std::string model_name;
      prm.enter_subsection ("Boundary heat flux model");
      {
        model_name = prm.get ("Model name");
      }
      prm.leave_subsection ();

      return std::get<dim>(registered_plugins).create_plugin (model_name,
                                                              "Boundary heat flux model::Model name");
    }



    template <int dim>
    void
    declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Boundary heat flux model");
      const std::string pattern_of_names
        = std::get<dim>(registered_plugins).get_pattern_of_names ();
      prm.declare_entry ("Model name", "function",
                         Patterns::Selection (pattern_of_names),
                         "Select one of the following plugins:\n\n"
                         +
                         std::get<dim>(registered_plugins).get_description_string());
      prm.leave_subsection ();

      std::get<dim>(registered_plugins).declare_parameters (prm);
    }



    template <int dim>
    void
    write_plugin_graph (std::ostream &out)
    {
      std::get<dim>(registered_plugins).write_plugin_graph ("Boundary heat flux interface",
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
      std::list<internal::Plugins::PluginList<BoundaryHeatFlux::Interface<2>>::PluginInfo> *
      internal::Plugins::PluginList<BoundaryHeatFlux::Interface<2>>::plugins = nullptr;
      template <>
      std::list<internal::Plugins::PluginList<BoundaryHeatFlux::Interface<3>>::PluginInfo> *
      internal::Plugins::PluginList<BoundaryHeatFlux::Interface<3>>::plugins = nullptr;
    }
  }

  namespace BoundaryHeatFlux
  {
#define INSTANTIATE(dim) \
  template class Interface<dim>; \
  \
  template \
  void \
  register_boundary_heat_flux<dim> (const std::string &, \
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
  write_plugin_graph<dim> (std::ostream &); \
  \
  template \
  std::unique_ptr<Interface<dim>> \
  create_boundary_heat_flux<dim> (ParameterHandler &prm);

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
