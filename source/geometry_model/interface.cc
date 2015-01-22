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
#include <aspect/geometry_model/interface.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/std_cxx1x/tuple.h>

namespace aspect
{
  namespace GeometryModel
  {
    template <int dim>
    Interface<dim>::~Interface ()
    {}


    template <int dim>
    void
    Interface<dim>::initialize ()
    {}


    template<int dim>
    std::map<std::string,types::boundary_id>
    Interface<dim>::get_symbolic_boundary_names_map() const
    {
      //return an empty map in the base class
      return std::map<std::string,types::boundary_id>();
    }


    template<int dim>
    std::set< std::pair< std::pair<types::boundary_id, types::boundary_id>, unsigned int > >
    Interface<dim>::get_periodic_boundary_pairs() const
    {
      //return an empty set in the base class
      return std::set< std::pair< std::pair< types::boundary_id, types::boundary_id>, unsigned int > >();
    }

    template <int dim>
    bool
    Interface<dim>::has_curved_elements() const
    {
      return true;
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


    /* --------- functions to translate between symbolic and numeric boundary indicators ------ */

    namespace
    {
      types::boundary_id
      translate_boundary_indicator (const std::string &name_,
                                    const std::map<std::string,types::boundary_id> &boundary_names_mapping)
      {
        // trim whitespace on either side of the name if necessary
        std::string name = name_;
        while ((name.size() > 0) && (name[0] == ' '))
          name.erase (name.begin());
        while ((name.size() > 0) && (name[name.size()-1] == ' '))
          name.erase (name.end()-1);

        // see if the given name is a symbolic one
        if (boundary_names_mapping.find (name) != boundary_names_mapping.end())
          return boundary_names_mapping.find(name)->second;
        else
          {
            // if it wasn't a symbolic name, it better be a number. we would
            // like to use Utilities::string_to_int, but as indicated by a
            // comment in that function, as of mid-2014 the function does not
            // do any error checking, so do it by hand here. (this was fixed
            // in late July 2014, so should work in deal.II 8.2.)
            //
            // since we test for errno, we need to make sure it is zero before
            // or otherwise the conversion may succeed and strtol will just
            // leave it where it was.
            char *p;
            errno = 0;
            const long int boundary_id = std::strtol(name.c_str(), &p, 10);
            if ((errno != 0) || (name.size() == 0) || ((name.size()>0) && (*p != '\0')))
              throw std::string ("Could not convert from string <") + name + "> to a boundary indicator.";

            // seems as if the conversion worked:
            return boundary_id;
          }
      }


      std::vector<types::boundary_id>
      translate_boundary_indicators (const std::vector<std::string> &names,
                                     const std::map<std::string,types::boundary_id> &boundary_names_mapping)
      {
        std::vector<types::boundary_id> results;
        for (unsigned int i=0; i<names.size(); ++i)
          results.push_back (translate_boundary_indicator(names[i], boundary_names_mapping));

        return results;
      }
    }


    template <int dim>
    types::boundary_id
    Interface<dim>::
    translate_symbolic_boundary_name_to_id (const std::string &name) const
    {
      return translate_boundary_indicator(name, get_symbolic_boundary_names_map());
    }



    template <int dim>
    std::vector<types::boundary_id>
    Interface<dim>::
    translate_symbolic_boundary_names_to_ids (const std::vector<std::string> &names) const
    {
      return translate_boundary_indicators(names, get_symbolic_boundary_names_map());
    }


    template <int dim>
    std::string
    Interface<dim>::
    translate_id_to_symbol_name(const types::boundary_id boundary_id) const
    {
      const std::map<std::string,types::boundary_id> mapping = get_symbolic_boundary_names_map();

      // loop over all entries in the map, and set 'name' to the key if the
      // value matches the given 'boundary_id'. if 'name' has already been
      // set, then this means that we had previously already found it -- i.e.,
      // that it is in the map at least twice. produce an error in that case.
      std::string name;
      for (std::map<std::string,types::boundary_id>::const_iterator p = mapping.begin();
           p != mapping.end(); ++p)
        if (p->second == boundary_id)
          {
            Assert (name == "",
                    ExcMessage ("This material model appears to provide multiple "
                                "names for the boundary with indicator <" +
                                Utilities::int_to_string (boundary_id) + ">."));
            name = p->first;
          }

      return name;
    }


// -------------------------------- Deal with registering geometry models and automating
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
    register_geometry_model (const std::string &name,
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
    create_geometry_model (ParameterHandler &prm)
    {
      std::string model_name;
      prm.enter_subsection ("Geometry model");
      {
        model_name = prm.get ("Model name");
      }
      prm.leave_subsection ();

      return std_cxx1x::get<dim>(registered_plugins).create_plugin (model_name,
                                                                    "Geometry model::model name",
                                                                    prm);
    }



    template <int dim>
    void
    declare_parameters (ParameterHandler &prm)
    {
      // declare the entry in the parameter file
      prm.enter_subsection ("Geometry model");
      {
        const std::string pattern_of_names
          = std_cxx1x::get<dim>(registered_plugins).get_pattern_of_names ();
        try
          {
            prm.declare_entry ("Model name", "",
                               Patterns::Selection (pattern_of_names),
                               "Select one of the following models:\n\n"
                               +
                               std_cxx1x::get<dim>(registered_plugins).get_description_string());
          }
        catch (const ParameterHandler::ExcValueDoesNotMatchPattern &)
          {
            // ignore the fact that the default value for this parameter
            // does not match the pattern
          }
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
      std::list<internal::Plugins::PluginList<GeometryModel::Interface<2> >::PluginInfo> *
      internal::Plugins::PluginList<GeometryModel::Interface<2> >::plugins = 0;

      template <>
      std::list<internal::Plugins::PluginList<GeometryModel::Interface<3> >::PluginInfo> *
      internal::Plugins::PluginList<GeometryModel::Interface<3> >::plugins = 0;
    }
  }

  namespace GeometryModel
  {
#define INSTANTIATE(dim) \
  template class Interface<dim>; \
  \
  template \
  void \
  register_geometry_model<dim> (const std::string &, \
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
  create_geometry_model<dim> (ParameterHandler &prm);

    ASPECT_INSTANTIATE(INSTANTIATE)
  }
}
