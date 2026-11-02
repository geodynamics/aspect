// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2015-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#include <aspect/global.h>
#include <aspect/boundary_fluid_pressure/interface.h>
#include <aspect/simulator_access.h>

#include <tuple>
#include <list>


namespace aspect
{
  namespace BoundaryFluidPressure
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
    register_boundary_fluid_pressure (const std::string &name,
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
    create_boundary_fluid_pressure (ParameterHandler &prm)
    {
      std::string model_name;
      prm.enter_subsection ("Boundary fluid pressure model");
      {
        model_name = prm.get ("Plugin name");
      }
      prm.leave_subsection ();

      return std::get<dim>(registered_plugins).create_plugin (model_name,
                                                              "Boundary fluid pressure model::Plugin name");
    }



    template <int dim>
    void
    declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Boundary fluid pressure model");
      const std::string pattern_of_names
        = std::get<dim>(registered_plugins).get_pattern_of_names ();
      prm.declare_entry ("Plugin name", "density",
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
      std::get<dim>(registered_plugins).write_plugin_graph ("Boundary fluid pressure interface",
                                                            out);
    }

  }
}

// explicit instantiations
namespace aspect
{
  namespace BoundaryFluidPressure
  {
#define INSTANTIATE(dim) \
  template class Interface<dim>; \
  \
  template \
  void \
  register_boundary_fluid_pressure<dim> (const std::string &, \
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
  create_boundary_fluid_pressure<dim> (ParameterHandler &prm);

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
