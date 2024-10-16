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
#include <aspect/geometry_model/initial_topography_model/interface.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/exceptions.h>
#include <tuple>
#include <list>

namespace aspect
{
  namespace InitialTopographyModel
  {
// -------------------------------- Deal with registering initial topography models and automating
// -------------------------------- their setup and selection at run time

    namespace
    {
      std::tuple
      <aspect::internal::Plugins::UnusablePluginList,
      aspect::internal::Plugins::UnusablePluginList,
      internal::Plugins::PluginList<Interface<2>>,
      internal::Plugins::PluginList<Interface<3>>> registered_plugins;
    }



    template <int dim>
    void
    register_initial_topography_model (const std::string &name,
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
    create_initial_topography_model (ParameterHandler &prm)
    {
      std::string model_name;
      prm.enter_subsection ("Geometry model");
      {
        prm.enter_subsection ("Initial topography model");
        {
          model_name = prm.get ("Model name");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();

      // If one sets the model name to an empty string in the input file,
      // ParameterHandler produces an error while reading the file. However,
      // if one omits specifying any model name at all (not even setting it to
      // the empty string) then the value we get here is the empty string. If
      // we don't catch this case here, we end up with awkward downstream
      // errors because the value obviously does not conform to the Pattern.
      AssertThrow(model_name != "unspecified",
                  ExcMessage("You need to select a Initial topography model "
                             "(`set Model name' in `subsection Initial topography model')."));

      return std::get<dim>(registered_plugins).create_plugin (model_name,
                                                              "Initial topography model::model name",
                                                              prm);
    }



    template <int dim>
    void
    declare_parameters (ParameterHandler &prm)
    {
      // declare the entry in the parameter file
      prm.enter_subsection ("Geometry model");
      {
        prm.enter_subsection ("Initial topography model");
        {
          const std::string pattern_of_names
            = std::get<dim>(registered_plugins).get_pattern_of_names ();
          prm.declare_entry ("Model name", "zero topography",
                             Patterns::Selection (pattern_of_names),
                             "Select one of the following models:\n\n"
                             +
                             std::get<dim>(registered_plugins).get_description_string());
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();

      std::get<dim>(registered_plugins).declare_parameters (prm);
    }



    template <int dim>
    void
    write_plugin_graph (std::ostream &out)
    {
      std::get<dim>(registered_plugins).write_plugin_graph ("Initial topography interface",
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
      std::list<internal::Plugins::PluginList<InitialTopographyModel::Interface<2>>::PluginInfo> *
      internal::Plugins::PluginList<InitialTopographyModel::Interface<2>>::plugins = nullptr;

      template <>
      std::list<internal::Plugins::PluginList<InitialTopographyModel::Interface<3>>::PluginInfo> *
      internal::Plugins::PluginList<InitialTopographyModel::Interface<3>>::plugins = nullptr;
    }
  }

  namespace InitialTopographyModel
  {
#define INSTANTIATE(dim) \
  template class Interface<dim>; \
  \
  template \
  void \
  register_initial_topography_model<dim> (const std::string &, \
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
  create_initial_topography_model<dim> (ParameterHandler &prm);

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
