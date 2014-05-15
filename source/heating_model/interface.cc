/*
  Copyright (C) 2011 - 2014 by the authors of the ASPECT code.

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
/*  $Id$  */


#include <aspect/global.h>
#include <aspect/heating_model/interface.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/std_cxx1x/tuple.h>

#include <list>


namespace aspect
{
  namespace HeatingModel
  {
    template <int dim>
    Interface<dim>::~Interface ()
    {}


    template <int dim>
    void
    Interface<dim>::initialize (const GeometryModel::Interface<dim> &geometry_model_)
    {
      geometry_model = &geometry_model_;
    }



    template <int dim>
    void
    Interface<dim>::set_current_time (const double time)
    {
      current_time = time;
    }



    template <int dim>
    void
    Interface<dim>::
    declare_parameters (dealii::ParameterHandler &prm)
    {}


    template <int dim>
    void
    Interface<dim>::parse_parameters (dealii::ParameterHandler &prm)
    {}


// -------------------------------- Deal with registering heating models and automating
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
    register_heating_model (const std::string &name,
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
    create_heating_model (const std::string &name,
                                         ParameterHandler &prm,
                                         const GeometryModel::Interface<dim> &geometry_model)
    {
      Interface<dim> *plugin = std_cxx1x::get<dim>(registered_plugins).create_plugin (name,
                                                                                      "Heating model",
                                                                                      prm);
      plugin->initialize (geometry_model);
      return plugin;
    }



    template <int dim>
    std::string
    get_names ()
    {
      return std_cxx1x::get<dim>(registered_plugins).get_pattern_of_names ();
    }


    template <int dim>
    void
    declare_parameters (ParameterHandler &prm)
    {
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
      std::list<internal::Plugins::PluginList<HeatingModel::Interface<2> >::PluginInfo> *
      internal::Plugins::PluginList<HeatingModel::Interface<2> >::plugins = 0;
      template <>
      std::list<internal::Plugins::PluginList<HeatingModel::Interface<3> >::PluginInfo> *
      internal::Plugins::PluginList<HeatingModel::Interface<3> >::plugins = 0;
    }
  }

  namespace HeatingModel
  {
#define INSTANTIATE(dim) \
  template class Interface<dim>; \
  \
  template \
  void \
  register_heating_model<dim> (const std::string &, \
                                                    const std::string &, \
                                                    void ( *) (ParameterHandler &), \
                                                    Interface<dim> *( *) ()); \
  \
  template  \
  void \
  declare_parameters<dim> (ParameterHandler &); \
  \
  template  \
  std::string \
  get_names<dim> (); \
  \
  template \
  Interface<dim> * \
  create_heating_model<dim> (const std::string &, \
                                            ParameterHandler &prm, \
                                            const GeometryModel::Interface<dim> &geometry_model);

    ASPECT_INSTANTIATE(INSTANTIATE)
  }
}
