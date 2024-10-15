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
#include <aspect/prescribed_stokes_solution/interface.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace PrescribedStokesSolution
  {
// -------------------------------- Deal with registering prescribed_stokes_solution models and automating
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
    register_prescribed_stokes_solution_model (const std::string &name,
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
    create_prescribed_stokes_solution (ParameterHandler &prm)
    {
      std::string model_name;
      prm.enter_subsection ("Prescribed Stokes solution");
      {
        model_name = prm.get ("Model name");
      }
      prm.leave_subsection ();

      if (model_name == "unspecified")
        return nullptr;

      std::unique_ptr<Interface<dim>>plugin = std::get<dim>(registered_plugins).create_plugin (model_name,
                                               "Prescribed Stokes solution::Model name");
      return plugin;
    }



    template <int dim>
    void
    declare_parameters (ParameterHandler &prm)
    {
      // declare the entry in the parameter file
      prm.enter_subsection ("Prescribed Stokes solution");
      {
        const std::string pattern_of_names
          = std::get<dim>(registered_plugins).get_pattern_of_names ();
        prm.declare_entry ("Model name", "unspecified",
                           Patterns::Selection (pattern_of_names+"|unspecified"),
                           "Select one of the following models:\n\n"
                           +
                           std::get<dim>(registered_plugins).get_description_string());
      }
      prm.leave_subsection ();

      std::get<dim>(registered_plugins).declare_parameters (prm);
    }



    template <int dim>
    void
    write_plugin_graph (std::ostream &out)
    {
      std::get<dim>(registered_plugins).write_plugin_graph ("Prescribed Stokes solution interface",
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
      std::list<internal::Plugins::PluginList<PrescribedStokesSolution::Interface<2>>::PluginInfo> *
      internal::Plugins::PluginList<PrescribedStokesSolution::Interface<2>>::plugins = nullptr;
      template <>
      std::list<internal::Plugins::PluginList<PrescribedStokesSolution::Interface<3>>::PluginInfo> *
      internal::Plugins::PluginList<PrescribedStokesSolution::Interface<3>>::plugins = nullptr;
    }
  }

  namespace PrescribedStokesSolution
  {
#define INSTANTIATE(dim) \
  template class Interface<dim>; \
  \
  template \
  void \
  register_prescribed_stokes_solution_model<dim> (const std::string &, \
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
  create_prescribed_stokes_solution<dim> (ParameterHandler &prm);

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
