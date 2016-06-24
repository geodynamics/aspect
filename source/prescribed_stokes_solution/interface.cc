/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include <aspect/global.h>
#include <aspect/prescribed_stokes_solution/interface.h>

namespace aspect
{
  namespace PrescribedStokesSolution
  {
    template <int dim>
    Interface<dim>::~Interface ()
    {}


    template <int dim>
    void
    Interface<dim>::initialize ()
    {}


    template <int dim>
    void
    Interface<dim>::update ()
    {}


    template <int dim>
    void
    Interface<dim>::
    declare_parameters (dealii::ParameterHandler &/*prm*/)
    {}


    template <int dim>
    void
    Interface<dim>::parse_parameters (dealii::ParameterHandler &/*prm*/)
    {}


// -------------------------------- Deal with registering prescribed_stokes_solution models and automating
// -------------------------------- their setup and selection at run time

    namespace
    {
      std_cxx1x::tuple
      <void *,
      void *,
      internal::Plugins::PluginList<Interface<2> >,
      internal::Plugins::PluginList<Interface<3> > > registered_plugins;
    }



    template <int dim>
    void
    register_prescribed_stokes_solution_model (const std::string &name,
                                               const std::string &description,
                                               void (*declare_parameters_function) (ParameterHandler &),
                                               Interface<dim> *(*factory_function) ())
    {
      std_cxx1x::get<dim>(registered_plugins).register_plugin (name,
                                                               description,
                                                               declare_parameters_function,
                                                               factory_function);
    }


    template <int dim>
    Interface<dim> *
    create_prescribed_stokes_solution (ParameterHandler &prm)
    {
      std::string model_name;
      prm.enter_subsection ("Prescribed Stokes solution");
      {
        model_name = prm.get ("Model name");
      }
      prm.leave_subsection ();

      if (model_name == "unspecified")
        return NULL;

      Interface<dim> *plugin = std_cxx1x::get<dim>(registered_plugins).create_plugin (model_name,
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
          = std_cxx1x::get<dim>(registered_plugins).get_pattern_of_names ();
        prm.declare_entry ("Model name", "unspecified",
                           Patterns::Selection (pattern_of_names+"|unspecified"),
                           "Select one of the following models:\n\n"
                           +
                           std_cxx1x::get<dim>(registered_plugins).get_description_string());
      }
      prm.leave_subsection ();

      std_cxx1x::get<dim>(registered_plugins).declare_parameters (prm);
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
      std::list<internal::Plugins::PluginList<PrescribedStokesSolution::Interface<2> >::PluginInfo> *
      internal::Plugins::PluginList<PrescribedStokesSolution::Interface<2> >::plugins = 0;
      template <>
      std::list<internal::Plugins::PluginList<PrescribedStokesSolution::Interface<3> >::PluginInfo> *
      internal::Plugins::PluginList<PrescribedStokesSolution::Interface<3> >::plugins = 0;
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
                                                  Interface<dim> *( *) ()); \
  \
  template  \
  void \
  declare_parameters<dim> (ParameterHandler &); \
  \
  template \
  Interface<dim> * \
  create_prescribed_stokes_solution<dim> (ParameterHandler &prm);

    ASPECT_INSTANTIATE(INSTANTIATE)
  }
}
