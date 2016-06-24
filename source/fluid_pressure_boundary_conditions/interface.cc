/*
  Copyright (C) 2015 - 2016 by the authors of the ASPECT code.

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
#include <aspect/fluid_pressure_boundary_conditions/interface.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/std_cxx11/tuple.h>

#include <list>


namespace aspect
{
  namespace FluidPressureBoundaryConditions
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
    Interface<dim>::
    declare_parameters (dealii::ParameterHandler &/*prm*/)
    {}


    template <int dim>
    void
    Interface<dim>::parse_parameters (dealii::ParameterHandler &/*prm*/)
    {}


// -------------------------------- Deal with registering models and automating
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
    register_fluid_pressure_boundary (const std::string &name,
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
    create_fluid_pressure_boundary (ParameterHandler &prm)
    {
      std::string model_name;
      prm.enter_subsection ("Boundary fluid pressure model");
      {
        model_name = prm.get ("Plugin name");
      }
      prm.leave_subsection ();

      return std_cxx1x::get<dim>(registered_plugins).create_plugin (model_name,
                                                                    "Boundary fluid pressure model::Plugin name");
    }



    template <int dim>
    void
    declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Boundary fluid pressure model");
      const std::string pattern_of_names
        = std_cxx1x::get<dim>(registered_plugins).get_pattern_of_names ();
      prm.declare_entry ("Plugin name", "density",
                         Patterns::Selection (pattern_of_names),
                         "Select one of the following plugins:\n\n"
                         +
                         std_cxx1x::get<dim>(registered_plugins).get_description_string());
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
      std::list<internal::Plugins::PluginList<FluidPressureBoundaryConditions::Interface<2> >::PluginInfo> *
      internal::Plugins::PluginList<FluidPressureBoundaryConditions::Interface<2> >::plugins = 0;
      template <>
      std::list<internal::Plugins::PluginList<FluidPressureBoundaryConditions::Interface<3> >::PluginInfo> *
      internal::Plugins::PluginList<FluidPressureBoundaryConditions::Interface<3> >::plugins = 0;
    }
  }

  namespace FluidPressureBoundaryConditions
  {
#define INSTANTIATE(dim) \
  template class Interface<dim>; \
  \
  template \
  void \
  register_fluid_pressure_boundary<dim> (const std::string &, \
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
  create_fluid_pressure_boundary<dim> (ParameterHandler &prm);

    ASPECT_INSTANTIATE(INSTANTIATE)
  }
}
