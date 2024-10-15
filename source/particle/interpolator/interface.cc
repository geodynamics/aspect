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

#include <aspect/particle/interpolator/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace Particle
  {
    namespace Interpolator
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
      register_particle_interpolator (const std::string &name,
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
      create_particle_interpolator (ParameterHandler &prm)
      {
        std::string name;
        name = prm.get ("Interpolation scheme");

        return std::get<dim>(registered_plugins).create_plugin (name,
                                                                "Particle::Interpolator name");
      }



      template <int dim>
      void
      declare_parameters (ParameterHandler &prm)
      {
        // declare the entry in the parameter file
        const std::string pattern_of_names
          = std::get<dim>(registered_plugins).get_pattern_of_names ();

        prm.declare_entry ("Interpolation scheme", "cell average",
                           Patterns::Selection (pattern_of_names),
                           "Select one of the following models:\n\n"
                           +
                           std::get<dim>(registered_plugins).get_description_string());

        std::get<dim>(registered_plugins).declare_parameters (prm);
      }



      template <int dim>
      void
      write_plugin_graph (std::ostream &out)
      {
        std::get<dim>(registered_plugins).write_plugin_graph ("Particle interpolator interface",
                                                              out);
      }
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
      std::list<internal::Plugins::PluginList<Particle::Interpolator::Interface<2>>::PluginInfo> *
      internal::Plugins::PluginList<Particle::Interpolator::Interface<2>>::plugins = nullptr;
      template <>
      std::list<internal::Plugins::PluginList<Particle::Interpolator::Interface<3>>::PluginInfo> *
      internal::Plugins::PluginList<Particle::Interpolator::Interface<3>>::plugins = nullptr;
    }
  }

  namespace Particle
  {
    namespace Interpolator
    {
#define INSTANTIATE(dim) \
  template class Interface<dim>; \
  \
  template \
  void \
  register_particle_interpolator<dim> (const std::string &, \
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
  create_particle_interpolator<dim> (ParameterHandler &prm);

      ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
    }
  }
}
