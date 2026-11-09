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
      create_particle_interpolator (const std::string &interpolator_name)
      {
        // 'bilinear least squares' is the deprecated old name of the 'linear least squares'
        // interpolator. The old name will be removed in the future.
        if (interpolator_name == "bilinear least squares")
          return std::get<dim>(registered_plugins).create_plugin ("linear least squares",
                                                                  "Particle::Interpolator name");
        else
          return std::get<dim>(registered_plugins).create_plugin (interpolator_name,
                                                                  "Particle::Interpolator name");
      }



      template <int dim>
      std::unique_ptr<Interface<dim>>
      create_particle_interpolator (ParameterHandler &prm)
      {
        const std::string name = prm.get ("Interpolation scheme");

        return create_particle_interpolator<dim> (name);
      }


      template <int dim>
      std::string
      get_valid_interpolator_names_pattern ()
      {
        return std::get<dim>(registered_plugins).get_pattern_of_names ();
      }


      template <int dim>
      void
      declare_parameters (ParameterHandler &prm)
      {
        // declare the entry in the parameter file
        const std::string pattern_of_names
          = get_valid_interpolator_names_pattern<dim>();

        // 'bilinear least squares' is the deprecated old name of the 'linear least squares'
        // interpolator. The old name will be removed in the future.
        prm.declare_entry ("Interpolation scheme", "cell average",
                           Patterns::Selection (pattern_of_names + "|bilinear least squares"),
                           "Select one of the following models:\n\n"
                           +
                           std::get<dim>(registered_plugins).get_description_string()
                           +
                           "\n\n"
                           "`bilinear least squares': "
                           "Deprecated, now an alias for `linear least squares'. "
                           "This alias will be removed in the future. ");

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
  template \
  std::string \
  get_valid_interpolator_names_pattern<dim> (); \
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
  create_particle_interpolator<dim> (const std::string &interpolator_name); \
  \
  template \
  std::unique_ptr<Interface<dim>> \
  create_particle_interpolator<dim> (ParameterHandler &prm);

      ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
    }
  }
}
