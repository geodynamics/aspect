/*
  Copyright (C) 2011, 2012 by the authors of the ASPECT code.

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


#include <aspect/initial_conditions/interface.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/std_cxx1x/tuple.h>

#include <list>


namespace aspect
{
  namespace InitialConditions
  {
    template <int dim>
    Interface<dim>::~Interface ()
    {}


    template <int dim>
    void
    Interface<dim>::initialize (const GeometryModel::Interface<dim>       &geometry_model_,
                                const BoundaryTemperature::Interface<dim> &boundary_temperature_,
                                const AdiabaticConditions<dim>            &adiabatic_conditions_)
    {
      geometry_model       = &geometry_model_;
      boundary_temperature = &boundary_temperature_;
      adiabatic_conditions = &adiabatic_conditions_;
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


// -------------------------------- Deal with registering initial_conditions models and automating
// -------------------------------- their setup and selection at run time

    namespace
    {
      internal::Plugins::PluginList<Interface<deal_II_dimension> > registered_plugins;
    }



    template <int dim>
    void
    register_initial_conditions_model (const std::string &name,
                                       const std::string &description,
                                       void (*declare_parameters_function) (ParameterHandler &),
                                       Interface<dim> *(*factory_function) ())
    {
      registered_plugins.register_plugin (name,
                                          description,
                                          declare_parameters_function,
                                          factory_function);
    }


    template <int dim>
    Interface<dim> *
    create_initial_conditions (ParameterHandler &prm,
                               const GeometryModel::Interface<dim> &geometry_model,
                               const BoundaryTemperature::Interface<dim> &boundary_temperature,
                               const AdiabaticConditions<dim>      &adiabatic_conditions)
    {
      std::string model_name;
      prm.enter_subsection ("Initial conditions");
      {
        model_name = prm.get ("Model name");
      }
      prm.leave_subsection ();

      Interface<dim> *plugin = registered_plugins.create_plugin (model_name, prm);
      plugin->initialize (geometry_model,
                          boundary_temperature,
                          adiabatic_conditions);
      return plugin;
    }



    void
    declare_parameters (ParameterHandler &prm)
    {
      // declare the entry in the parameter file
      prm.enter_subsection ("Initial conditions");
      {
        const std::string pattern_of_names
          = registered_plugins.get_pattern_of_names ();
        prm.declare_entry ("Model name", "",
                           Patterns::Selection (pattern_of_names),
                           "Select one of the following models:\n\n"
                           +
                           registered_plugins.get_description_string());
      }
      prm.leave_subsection ();

      registered_plugins.declare_parameters (prm);
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
      std::list<internal::Plugins::PluginList<InitialConditions::Interface<deal_II_dimension> >::PluginInfo> *
      internal::Plugins::PluginList<InitialConditions::Interface<deal_II_dimension> >::plugins = 0;
    }
  }

  namespace InitialConditions
  {
    template class Interface<deal_II_dimension>;

    template
    void
    register_initial_conditions_model<deal_II_dimension> (const std::string &,
                                                          const std::string &,
                                                          void ( *) (ParameterHandler &),
                                                          Interface<deal_II_dimension> *( *) ());

    template
    Interface<deal_II_dimension> *
    create_initial_conditions<deal_II_dimension> (ParameterHandler &prm,
                                                  const GeometryModel::Interface<deal_II_dimension> &geometry_model,
                                                  const BoundaryTemperature::Interface<deal_II_dimension> &boundary_temperature,
                                                  const AdiabaticConditions<deal_II_dimension>      &adiabatic_conditions);
  }
}
