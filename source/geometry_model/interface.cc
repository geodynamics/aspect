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
#include <aspect/geometry_model/interface.h>
#include <aspect/simulator_access.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/base/utilities.h>

#include <tuple>

namespace aspect
{
  namespace GeometryModel
  {
    template <int dim>
    std::map<std::string,types::boundary_id>
    Interface<dim>::get_symbolic_boundary_names_map() const
    {
      // return an empty map in the base class
      return {};
    }



    template <int dim>
    std::set<std::pair<std::pair<types::boundary_id, types::boundary_id>, unsigned int>>
    Interface<dim>::get_periodic_boundary_pairs() const
    {
      // return an empty set in the base class
      return {};
    }



    template <int dim>
    void
    Interface<dim>::adjust_positions_for_periodicity (Point<dim> &/*position*/,
                                                      const ArrayView<Point<dim>> &/*connected_positions*/,
                                                      const ArrayView<Tensor<1, dim>> &/*connected_velocities*/) const
    {
      AssertThrow(false,
                  ExcMessage("Positions cannot be adjusted for periodicity in the chosen geometry model."));
      return;
    }



    template <int dim>
    bool
    Interface<dim>::has_curved_elements() const
    {
      return true;
    }



    template <int dim>
    std::array<double,dim>
    Interface<dim>::cartesian_to_natural_coordinates(const Point<dim> &) const
    {
      Assert (false,
              ExcMessage ("The cartesian_to_natural_coordinates function has "
                          "not been implemented in this geometry model."));
      return {};
    }


    template <int dim>
    Utilities::NaturalCoordinate<dim>
    Interface<dim>::cartesian_to_other_coordinates(const Point<dim> &position,
                                                   const Utilities::Coordinates::CoordinateSystem &coordinate_system) const
    {
      std::array<double, dim> other_coord;
      switch (coordinate_system)
        {
          case Utilities::Coordinates::cartesian:
            other_coord = Utilities::convert_point_to_array(position);
            break;

          case Utilities::Coordinates::spherical:
            other_coord = Utilities::Coordinates::cartesian_to_spherical_coordinates(position);
            break;

          case Utilities::Coordinates::depth:
            other_coord[0] = depth(position);
            break;

          case Utilities::Coordinates::ellipsoidal:
          default:
            AssertThrow(false, ExcNotImplemented());
        }

      return Utilities::NaturalCoordinate<dim>(other_coord, coordinate_system);
    }

    template <int dim>
    Point<dim>
    Interface<dim>::natural_to_cartesian_coordinates(const std::array<double,dim> &) const
    {
      Assert (false,
              ExcMessage ("The natural_to_cartesian_coordinates function has "
                          "not been implemented in this geometry model."));
      return Point<dim>();
    }


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

        // backwards compatibility (rename boundaries to all use "top" and "bottom"
        if (name == "surface" || name == "outer")
          name = "top";
        else if (name == "inner")
          name = "bottom";

        // see if the given name is a symbolic one
        if (boundary_names_mapping.find (name) != boundary_names_mapping.end())
          return boundary_names_mapping.find(name)->second;
        else
          return dealii::Utilities::string_to_int(name);
      }


      std::vector<types::boundary_id>
      translate_boundary_indicators (const std::vector<std::string> &names,
                                     const std::map<std::string,types::boundary_id> &boundary_names_mapping)
      {
        std::vector<types::boundary_id> results;
        results.reserve(names.size());
        for (const auto &name : names)
          results.push_back (translate_boundary_indicator(name, boundary_names_mapping));

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
      for (const auto &p : mapping)
        if (p.second == boundary_id)
          {
            Assert (name == "",
                    ExcMessage ("This geometry model appears to provide multiple "
                                "names for the boundary with indicator <" +
                                Utilities::int_to_string (boundary_id) + ">."));
            name = p.first;
          }

      return name;
    }


// -------------------------------- Deal with registering geometry models and automating
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
    register_geometry_model (const std::string &name,
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
    create_geometry_model (ParameterHandler &prm)
    {
      std::string model_name;
      prm.enter_subsection ("Geometry model");
      {
        model_name = prm.get ("Model name");
      }
      prm.leave_subsection ();

      // If one sets the model name to an empty string in the input file,
      // ParameterHandler produces an error while reading the file. However,
      // if one omits specifying any model name at all (not even setting it to
      // the empty string) then the value we get here is the empty string. If
      // we don't catch this case here, we end up with awkward downstream
      // errors because the value obviously does not conform to the Pattern.
      AssertThrow(model_name != "unspecified",
                  ExcMessage("You need to select a Geometry model "
                             "(`set Model name' in `subsection Geometry model')."));

      return std::get<dim>(registered_plugins).create_plugin (model_name,
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
      std::get<dim>(registered_plugins).write_plugin_graph ("Geometry model interface",
                                                            out);
    }



    template <int dim>
    void
    Interface<dim>::make_periodicity_constraints(const DoFHandler<dim> &dof_handler,
                                                 AffineConstraints<double> &constraints) const
    {
      using periodic_boundary_set
        = std::set<std::pair<std::pair<types::boundary_id, types::boundary_id>, unsigned int>>;
      periodic_boundary_set pbs = get_periodic_boundary_pairs();

      for (const auto &pb : pbs)
        {
          DoFTools::make_periodicity_constraints(dof_handler,
                                                 pb.first.first,  // first boundary id
                                                 pb.first.second, // second boundary id
                                                 pb.second,       // cartesian direction for translational symmetry
                                                 constraints);
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
      std::list<internal::Plugins::PluginList<GeometryModel::Interface<2>>::PluginInfo> *
      internal::Plugins::PluginList<GeometryModel::Interface<2>>::plugins = nullptr;

      template <>
      std::list<internal::Plugins::PluginList<GeometryModel::Interface<3>>::PluginInfo> *
      internal::Plugins::PluginList<GeometryModel::Interface<3>>::plugins = nullptr;
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
  create_geometry_model<dim> (ParameterHandler &prm);

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
