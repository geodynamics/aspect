/*
  Copyright (C) 2011 - 2020 by the authors of the ASPECT code.

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
#include <aspect/utilities.h>
#include <aspect/simulator_access.h>

#ifdef ASPECT_WITH_LIBDAP
#include <D4Connect.h>
#include <Connect.h>
#include <Response.h>
#include <Array.h>
#endif

#include <array>
#include <deal.II/base/point.h>

#include <deal.II/base/mpi.h>
#include <deal.II/base/table.h>
#include <deal.II/base/table_indices.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/patterns.h>


#include <aspect/geometry_model/box.h>
#include <aspect/geometry_model/two_merged_boxes.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/geometry_model/sphere.h>
#include <aspect/geometry_model/chunk.h>
#include <aspect/geometry_model/initial_topography_model/ascii_data.h>

#include <fstream>
#include <string>
#include <locale>
#include <dirent.h>
#include <sys/stat.h>
#include <errno.h>

#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <boost/lexical_cast.hpp>

namespace aspect
{
  /**
   * A namespace for utility functions that might be used in many different
   * places to prevent code duplication.
   */
  namespace Utilities
  {
    template <typename T>
    Table<2,T>
    parse_input_table (const std::string &input_string,
                       const unsigned int n_rows,
                       const unsigned int n_columns,
                       const std::string &property_name)
    {
      Table<2,T> input_table(n_rows,n_columns);

      const std::vector<std::string> rows = Utilities::possibly_extend_from_1_to_N(Utilities::split_string_list(input_string,';'),
                                                                                   n_rows,
                                                                                   property_name);

      for (unsigned int i=0; i<rows.size(); ++i)
        {
          std::vector<std::string> current_columns = Utilities::possibly_extend_from_1_to_N(Utilities::split_string_list(rows[i]),
                                                     n_columns,
                                                     property_name);

          for (unsigned int j=0; j<current_columns.size(); ++j)
            {
              // get rid of surrounding whitespace
              trim(current_columns[j]);

              input_table[i][j] = boost::lexical_cast<T>(current_columns[j]);
            }
        }

      return input_table;
    }

    namespace
    {
      // This is a helper function used in parse_map_to_double_array below.
      // It takes an input_string that is expected to follow the input format
      // explained in the documentation of the parse_map_to_double_array function
      // and parses it into a multimap, only performing rudimentary error checking
      // for correct formatting.
      std::multimap<std::string, double>
      parse_string_to_map (const std::string &input_string,
                           const std::vector<std::string> &list_of_keys,
                           const std::string &property_name)
      {
        std::multimap<std::string, double> parsed_map;

        if (Patterns::Map(Patterns::Anything(),
                          Patterns::List(Patterns::Double(),
                                         0,
                                         std::numeric_limits<unsigned int>::max(),
                                         "|")).match(input_string))
          {
            // Split the list by comma delimited components.
            const std::vector<std::string> field_entries = dealii::Utilities::split_string_list(input_string, ',');

            for (const auto &field_entry : field_entries)
              {
                // Split each entry into string and value ( <id> : <value>)
                std::vector<std::string> key_and_value = Utilities::split_string_list (field_entry, ':');

                // Ensure that each entry has the correct form.
                AssertThrow (key_and_value.size() == 2,
                             ExcMessage ("The format for mapped "
                                         + property_name
                                         + "requires that each entry has the "
                                         "form `<key> : <value>' "
                                         ", but the entry <"
                                         + field_entry
                                         + "> does not appear to follow this pattern."));

                // Handle special key "all", which must be the only entry if found
                if (key_and_value[0] == "all")
                  {
                    AssertThrow (field_entries.size() == 1,
                                 ExcMessage ("The keyword `all' in the property "
                                             + property_name
                                             + " is only allowed if there is no other "
                                             "keyword."));

                    const std::vector<std::string> values = dealii::Utilities::split_string_list(key_and_value[1], '|');

                    // Assign all the values to all fields
                    for (const std::string &key: list_of_keys)
                      for (const std::string &value : values)
                        {
                          parsed_map.emplace(key, Utilities::string_to_double(value));
                        }
                  }
                // Handle lists of multiple unique entries
                else
                  {
                    AssertThrow (parsed_map.find(key_and_value[0]) == parsed_map.end(),
                                 ExcMessage ("The keyword <"
                                             + key_and_value[0]
                                             + "> in "
                                             + property_name
                                             + " is listed multiple times. "
                                             "Check that you have only one value for "
                                             "each field id in your list."));

                    const std::vector<std::string> values = dealii::Utilities::split_string_list(key_and_value[1], '|');

                    for (const auto &value : values)
                      {
                        parsed_map.emplace(key_and_value[0],Utilities::string_to_double(value));
                      }
                  }
              }
          }
        else if (Patterns::List(Patterns::Double(),1,list_of_keys.size()).match(input_string))
          {
            // Handle the format of a comma separated list of doubles, with no keywords
            const std::vector<double> values = possibly_extend_from_1_to_N (dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(input_string)),
                                                                            list_of_keys.size(),
                                                                            property_name);

            for (unsigned int i=0; i<values.size(); ++i)
              {
                // list_of_keys and values have the same length, which is guaranteed by the
                // call to possibly_extend_from_1_to_N() above
                parsed_map.emplace(list_of_keys[i],values[i]);
              }
          }
        else
          {
            // No Patterns matches were found!
            AssertThrow (false,
                         ExcMessage ("The required format for property <"
                                     + property_name
                                     + "> was not found. Specify a comma separated "
                                     + "list of `<double>' or `<key1> : <double>|<double>|..., "
                                     + "<key2> : <double>|... , ... '."));
          }

        return parsed_map;
      }
    }


    std::vector<double>
    parse_map_to_double_array (const std::string &input_string,
                               const std::vector<std::string> &list_of_keys,
                               const bool expects_background_field,
                               const std::string &property_name,
                               const bool allow_multiple_values_per_key,
                               const std::shared_ptr<std::vector<unsigned int> > &n_values_per_key,
                               const bool allow_missing_keys)
    {
      std::vector<std::string> field_names = list_of_keys;
      if (expects_background_field)
        field_names.insert(field_names.begin(),"background");
      const unsigned int n_fields = field_names.size();

      // First: parse the string into a map depending on what Pattern we are dealing with
      std::multimap<std::string, double> parsed_map = parse_string_to_map(input_string,
                                                                          field_names,
                                                                          property_name);

      // Second: Now check that the structure of the map is as expected
      {
        const bool check_structure = (n_values_per_key && n_values_per_key->size() != 0);
        const bool store_structure = (n_values_per_key && n_values_per_key->size() == 0);
        std::vector<unsigned int> values_per_key(n_fields, 0);

        if (check_structure)
          AssertThrow(n_values_per_key->size() == n_fields,
                      ExcMessage("When providing an expected structure for input parameter " + property_name + " you need to provide "
                                 + "as many entries in the structure vector as there are input field names (+1 if there is a background field). "
                                 + "The current structure vector has " + std::to_string(n_values_per_key->size()) + " entries, but there are "
                                 + std::to_string(n_fields) + " field names." ));

        for (const std::pair<const std::string, double> &key_and_value: parsed_map)
          {
            const std::vector<std::string>::iterator field_name =
              std::find(field_names.begin(),field_names.end(),key_and_value.first);

            // Ensure that each key is in the list of field names
            AssertThrow (field_name != field_names.end(),
                         ExcMessage ("The keyword <" + key_and_value.first + "> in "
                                     + property_name + " does not match any entries "
                                     "from the list of field names"
                                     + ((expects_background_field)
                                        ?
                                        " (plus `background' for the background field). "
                                        :
                                        ". ")
                                     + "Check that you only use valid names.\n\n"
                                     "One example of where to check this is if "
                                     "Compositional fields are used, "
                                     "then check the id list "
                                     "from `set Names of fields' in the "
                                     "Compositional fields subsection. "
                                     "Alternatively, if `set Names of fields' "
                                     "is not set, the default names are "
                                     "C_1, C_2, ..., C_n."));

            const unsigned int field_index = std::distance(field_names.begin(), field_name);
            values_per_key[field_index] += 1;
          }

        if (store_structure)
          *n_values_per_key = values_per_key;

        unsigned int field_index = 0;
        for (const unsigned int &n_values: values_per_key)
          {
            if (allow_multiple_values_per_key == false)
              AssertThrow (n_values <= 1,
                           ExcMessage ("The keyword <"
                                       + field_names[field_index]
                                       + "> in "
                                       + property_name
                                       + " has multiple values, which is unexpected. "
                                       "Check that you have only one value for "
                                       "each field id in your list."));

            if (allow_missing_keys == false)
              AssertThrow (n_values > 0,
                           ExcMessage ("The keyword <"
                                       + field_names[field_index]
                                       + "> in "
                                       + property_name
                                       + " is not listed, although it is expected. "
                                       "Check that you have at least one value for "
                                       "each field id in your list (possibly plus "
                                       "`background` if a background field is expected "
                                       "for this property)."));

            if (check_structure)
              {
                AssertThrow(((*n_values_per_key)[field_index] == n_values || n_values == 1),
                            ExcMessage("The key <" + field_names[field_index] + "> in <"+ property_name + "> does not have "
                                       + "the expected number of values. It expects " + std::to_string((*n_values_per_key)[field_index])
                                       + "or 1 values, but we found " + std::to_string(n_values) + " values."));
                if (n_values == 1)
                  {
                    const std::string field_name = field_names[field_index];
                    const double field_value = parsed_map.find(field_name)->second;
                    for (unsigned int i=1; i<(*n_values_per_key)[field_index]; ++i)
                      parsed_map.emplace(field_name, field_value);
                  }
              }

            ++field_index;
          }
      }

      // Finally: Convert the map into a vector of doubles, sorted in the order
      // of the field_names input parameter
      std::vector<double> return_values;
      for (const std::string &field_name: field_names)
        {
          const std::pair<std::multimap<std::string, double>::const_iterator,
                std::multimap<std::string, double>::const_iterator> entry_range = parsed_map.equal_range(field_name);

          for (auto entry = entry_range.first; entry != entry_range.second; ++entry)
            return_values.push_back(entry->second);
        }
      return return_values;
    }


#if !DEAL_II_VERSION_GTE(9,2,0)
    /**
     * Split the set of DoFs (typically locally owned or relevant) in @p whole_set into blocks
     * given by the @p dofs_per_block structure.
     */
    void split_by_block (const std::vector<types::global_dof_index> &dofs_per_block,
                         const IndexSet &whole_set,
                         std::vector<IndexSet> &partitioned)
    {
      const unsigned int n_blocks = dofs_per_block.size();
      partitioned.clear();

      partitioned.resize(n_blocks);
      types::global_dof_index start = 0;
      for (unsigned int i=0; i<n_blocks; ++i)
        {
          partitioned[i] = whole_set.get_view(start, start + dofs_per_block[i]);
          start += dofs_per_block[i];
        }
    }
#endif

    template <int dim>
    std::vector<std::string>
    expand_dimensional_variable_names (const std::vector<std::string> &var_declarations)
    {
      std::string dim_names[3] = {"x", "y", "z"};
      char fn_split = '(', fn_end = ')';
      std::vector<std::string> var_name_list;

      for (const auto &var_decl : var_declarations)
        {
          if (var_decl.find(fn_split) != std::string::npos && var_decl[var_decl.length()-1]==fn_end)
            {
              const std::string fn_name = var_decl.substr(0, var_decl.find(fn_split));

              // Cannot be const because will be manipulated to strip whitespace
              std::string var_name = var_decl.substr(var_decl.find(fn_split)+1, var_decl.length()-2-var_decl.find(fn_split));
              while ((var_name.length() != 0) && (var_name[0] == ' '))
                var_name.erase(0, 1);
              while ((var_name.length() != 0) && (var_name[var_name.length()-1] == ' '))
                var_name.erase(var_name.length()-1, 1);
              if (fn_name == "vector")
                {
                  for (int i=0; i<dim; ++i)
                    var_name_list.push_back(var_name+"_"+dim_names[i]);
                }
              else if (fn_name == "tensor")
                {
                  for (int i=0; i<dim; ++i)
                    for (int j=0; j< dim; ++j)
                      var_name_list.push_back(var_name+"_"+dim_names[i]+dim_names[j]);
                }
              else
                {
                  var_name_list.push_back(var_decl);
                }
            }
          else
            {
              var_name_list.push_back(var_decl);
            }
        }
      return var_name_list;
    }


    /**
    * This is an internal deal.II function stolen from dof_tools.cc
    *
    * Return an array that for each dof on the reference cell lists the
    * corresponding vector component.
    */
    template <int dim, int spacedim>
    std::vector<unsigned char>
    get_local_component_association (const FiniteElement<dim,spacedim>  &fe,
                                     const ComponentMask        & /*component_mask*/)
    {
      const unsigned char invalid = static_cast<unsigned char>(-1);
      std::vector<unsigned char> local_component_association (fe.dofs_per_cell, invalid);

      // compute the component each local dof belongs to.
      // if the shape function is primitive, then this
      // is simple and we can just associate it with
      // what system_to_component_index gives us
      for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
        {
          // see the deal.II version if we ever need non-primitive FEs.
          Assert (fe.is_primitive(i), ExcNotImplemented());
          local_component_association[i] =
            fe.system_to_component_index(i).first;
        }

      Assert (std::find (local_component_association.begin(),
                         local_component_association.end(),
                         invalid)
              ==
              local_component_association.end(),
              ExcInternalError());

      return local_component_association;
    }


    template <int dim>
    IndexSet extract_locally_active_dofs_with_component(const DoFHandler<dim> &dof_handler,
                                                        const ComponentMask &component_mask)
    {
      std::vector<unsigned char> local_asoc =
        get_local_component_association (dof_handler.get_fe(),
                                         ComponentMask(dof_handler.get_fe().n_components(), true));

      IndexSet ret(dof_handler.n_dofs());

      unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;
      std::vector<types::global_dof_index> indices(dofs_per_cell);
      for (const auto &cell : dof_handler.active_cell_iterators())
        if (cell->is_locally_owned())
          {
            cell->get_dof_indices(indices);
            for (unsigned int i=0; i<dofs_per_cell; ++i)
              if (component_mask[local_asoc[i]])
                ret.add_index(indices[i]);
          }

      return ret;
    }



    namespace Coordinates
    {

      template <int dim>
      std::array<double,dim>
      WGS84_coordinates(const Point<dim> &position)
      {
        std::array<double,dim> ecoord;

        // Define WGS84 ellipsoid constants.
        const double radius = 6378137.;
        const double ellipticity = 8.1819190842622e-2;
        const double b = std::sqrt(radius * radius
                                   * (1 - ellipticity * ellipticity));
        const double ep = std::sqrt((radius * radius - b * b) / (b * b));
        const double p = std::sqrt(position(0) * position(0) + position(1) * position(1));
        const double th = std::atan2(radius * position(2), b * p);
        ecoord[2] = std::atan2((position(2) + ep * ep * b * std::sin(th)
                                * std::sin(th) * std::sin(th)),
                               (p - (ellipticity * ellipticity * radius * (std::cos(th)
                                                                           * std::cos(th) * std::cos(th)))))
                    * (180. / numbers::PI);

        if (dim == 3)
          {
            ecoord[1] = std::atan2(position(1), position(0))
                        * (180. / numbers::PI);

            /* Set all longitudes between [0,360]. */
            if (ecoord[1] < 0.)
              ecoord[1] += 360.;
            else if (ecoord[1] > 360.)
              ecoord[1] -= 360.;
          }
        else
          ecoord[1] = 0.0;


        ecoord[0] = radius/std::sqrt(1- ellipticity * ellipticity
                                     * std::sin(numbers::PI * ecoord[2]/180)
                                     * std::sin(numbers::PI * ecoord[2]/180));
        return ecoord;
      }

      template <int dim>
      std::array<double,dim>
      cartesian_to_spherical_coordinates(const Point<dim> &position)
      {
        std::array<double,dim> scoord;

        scoord[0] = position.norm(); // R
        scoord[1] = std::atan2(position(1),position(0)); // Phi
        if (scoord[1] < 0.0)
          scoord[1] += 2.0*numbers::PI; // correct phi to [0,2*pi]
        if (dim==3)
          {
            if (scoord[0] > std::numeric_limits<double>::min())
              scoord[2] = std::acos(position(2)/scoord[0]);
            else
              scoord[2] = 0.0;
          }
        return scoord;
      }

      template <int dim>
      Point<dim>
      spherical_to_cartesian_coordinates(const std::array<double,dim> &scoord)
      {
        Point<dim> ccoord;

        switch (dim)
          {
            case 2:
            {
              ccoord[0] = scoord[0] * std::cos(scoord[1]); // X
              ccoord[1] = scoord[0] * std::sin(scoord[1]); // Y
              break;
            }
            case 3:
            {
              ccoord[0] = scoord[0] * std::sin(scoord[2]) * std::cos(scoord[1]); // X
              ccoord[1] = scoord[0] * std::sin(scoord[2]) * std::sin(scoord[1]); // Y
              ccoord[2] = scoord[0] * std::cos(scoord[2]); // Z
              break;
            }
            default:
              Assert (false, ExcNotImplemented());
              break;
          }

        return ccoord;
      }

      template <int dim>
      std::array<double,3>
      cartesian_to_ellipsoidal_coordinates(const Point<3> &x,
                                           const double semi_major_axis_a,
                                           const double eccentricity)
      {
        const double R    = semi_major_axis_a;
        const double b      = std::sqrt(R * R * (1 - eccentricity * eccentricity));
        const double ep     = std::sqrt((R * R - b * b) / (b * b));
        const double p      = std::sqrt(x(0) * x(0) + x(1) * x(1));
        const double th     = std::atan2(R * x(2), b * p);
        const double phi    = std::atan2(x(1), x(0));
        const double theta  = std::atan2(x(2) + ep * ep * b * std::pow(std::sin(th),3),
                                         (p - (eccentricity * eccentricity * R  * std::pow(std::cos(th),3))));
        const double R_bar = R / (std::sqrt(1 - eccentricity * eccentricity * std::sin(theta) * std::sin(theta)));
        const double R_plus_d = p / std::cos(theta);

        std::array<double,3> phi_theta_d;
        phi_theta_d[0] = phi;

        phi_theta_d[1] = theta;
        phi_theta_d[2] = R_plus_d - R_bar;
        return phi_theta_d;
      }

      template <int dim>
      Point<3>
      ellipsoidal_to_cartesian_coordinates(const std::array<double,3> &phi_theta_d,
                                           const double semi_major_axis_a,
                                           const double eccentricity)
      {
        const double phi   = phi_theta_d[0];
        const double theta = phi_theta_d[1];
        const double d     = phi_theta_d[2];

        const double R_bar = semi_major_axis_a / std::sqrt(1 - (eccentricity * eccentricity *
                                                                std::sin(theta) * std::sin(theta)));

        return Point<3> ((R_bar + d) * std::cos(phi) * std::cos(theta),
                         (R_bar + d) * std::sin(phi) * std::cos(theta),
                         ((1 - eccentricity * eccentricity) * R_bar + d) * std::sin(theta));

      }


      template <int dim>
      Tensor<1, dim>
      spherical_to_cartesian_vector(const Tensor<1, dim> &spherical_vector,
                                    const Point<dim> &position)
      {
        Tensor<1, dim> cartesian_vector;

        const std::array<double, dim> r_phi_theta = cartesian_to_spherical_coordinates(position);

        switch (dim)
          {
            case 2:
            {
              const double phi = r_phi_theta[1];

              const double u_r   = spherical_vector[0];
              const double u_phi = spherical_vector[1];

              cartesian_vector[0] = std::cos(phi)*u_r
                                    - std::sin(phi)*u_phi; // X
              cartesian_vector[1] = std::sin(phi)*u_r
                                    + std::cos(phi)*u_phi; // Y

              break;
            }
            case 3:
            {
              const double phi   = r_phi_theta[1];
              const double theta = r_phi_theta[2];

              const double u_r     = spherical_vector[0];
              const double u_phi   = spherical_vector[1];
              const double u_theta = spherical_vector[2];

              cartesian_vector[0] = std::cos(phi)*std::sin(theta)*u_r
                                    - std::sin(phi)*u_phi
                                    + std::cos(phi)*std::cos(theta)*u_theta; // X
              cartesian_vector[1] = std::sin(phi)*std::sin(theta)*u_r
                                    + std::cos(phi)*u_phi
                                    + std::sin(phi)*std::cos(theta)*u_theta; // Y
              cartesian_vector[2] = std::cos(theta)*u_r
                                    - std::sin(theta)*u_theta; // Z
              break;
            }

            default:
              Assert (false, ExcNotImplemented());
              break;
          }

        return cartesian_vector;
      }



      CoordinateSystem
      string_to_coordinate_system(const std::string &coordinate_system)
      {
        if (coordinate_system == "cartesian")
          return cartesian;
        else if (coordinate_system == "spherical")
          return spherical;
        else if (coordinate_system == "depth")
          return Coordinates::depth;
        else
          AssertThrow(false, ExcNotImplemented());

        return Coordinates::invalid;
      }


    }

    template <int dim>
    bool
    polygon_contains_point(const std::vector<Point<2> > &point_list,
                           const dealii::Point<2> &point)
    {
      /**
       * This code has been based on http://geomalgorithms.com/a03-_inclusion.html,
       * and therefore requires the following copyright notice:
       *
       * Copyright 2000 softSurfer, 2012 Dan Sunday
       * This code may be freely used and modified for any purpose
       * providing that this copyright notice is included with it.
       * SoftSurfer makes no warranty for this code, and cannot be held
       * liable for any real or imagined damage resulting from its use.
       * Users of this code must verify correctness for their application.
       *
       * The main functional difference between the original code and this
       * code is that all the boundaries are considered to be inside the
       * polygon. One should of course realize that with floating point
       * arithmetic no guarantees can be made for the borders, but for
       * exact arithmetic this algorithm would work (also see polygon
       * in point test).
       */
      int pointNo = point_list.size();
      int    wn = 0;    // the  winding number counter
      int   j=pointNo-1;

      // loop through all edges of the polygon
      for (int i=0; i<pointNo; i++)
        {
          // edge from V[i] to  V[i+1]
          if (point_list[j][1] <= point[1])
            {
              // start y <= P.y
              if (point_list[i][1] >= point[1])      // an upward crossing
                {
                  const double is_left = (point_list[i][0] - point_list[j][0]) * (point[1] - point_list[j][1])
                                         - (point[0] -  point_list[j][0]) * (point_list[i][1] - point_list[j][1]);

                  if ( is_left > 0 && point_list[i][1] > point[1])
                    {
                      // P left of  edge
                      ++wn;            // have  a valid up intersect
                    }
                  else if ( is_left == 0)
                    {
                      // The point is exactly on the infinite line.
                      // determine if it is on the segment
                      const double dot_product = (point - point_list[j])*(point_list[i] - point_list[j]);

                      if (dot_product >= 0)
                        {
                          const double squaredlength = (point_list[i] - point_list[j]).norm_square();

                          if (dot_product <= squaredlength)
                            {
                              return true;
                            }
                        }
                    }
                }
            }
          else
            {
              // start y > P.y (no test needed)
              if (point_list[i][1]  <= point[1])     // a downward crossing
                {
                  const double is_left = (point_list[i][0] - point_list[j][0]) * (point[1] - point_list[j][1])
                                         - (point[0] -  point_list[j][0]) * (point_list[i][1] - point_list[j][1]);

                  if ( is_left < 0)
                    {
                      // P right of  edge
                      --wn;            // have  a valid down intersect
                    }
                  else if ( is_left == 0)
                    {
                      // This code is to make sure that the boundaries are included in the polygon.
                      // The point is exactly on the infinite line.
                      // determine if it is on the segment
                      const double dot_product = (point - point_list[j])*(point_list[i] - point_list[j]);

                      if (dot_product >= 0)
                        {
                          const double squaredlength = (point_list[i] - point_list[j]).norm_square();

                          if (dot_product <= squaredlength)
                            {
                              return true;
                            }
                        }
                    }
                }
            }
          j=i;
        }

      return (wn != 0);
    }

    template <int dim>
    double
    signed_distance_to_polygon(const std::vector<Point<2> > &point_list,
                               const dealii::Point<2> &point)
    {
      // If the point lies outside polygon, we give it a negative sign,
      // inside a positive sign.
      const double sign = polygon_contains_point<dim>(point_list, point) ? 1.0 : -1.0;

      /**
       * This code is based on http://geomalgorithms.com/a02-_lines.html#Distance-to-Infinite-Line,
       * and therefore requires the following copyright notice:
       *
       * Copyright 2000 softSurfer, 2012 Dan Sunday
       * This code may be freely used and modified for any purpose
       * providing that this copyright notice is included with it.
       * SoftSurfer makes no warranty for this code, and cannot be held
       * liable for any real or imagined damage resulting from its use.
       * Users of this code must verify correctness for their application.
       *
       */

      const unsigned int n_poly_points = point_list.size();
      AssertThrow(n_poly_points >= 3, ExcMessage("Not enough polygon points were specified."));

      // Initialize a vector of distances for each point of the polygon with a very large distance
      std::vector<double> distances(n_poly_points, 1e23);

      // Create another polygon but with all points shifted 1 position to the right
      std::vector<Point<2> > shifted_point_list(n_poly_points);
      shifted_point_list[0] = point_list[n_poly_points-1];

      for (unsigned int i = 0; i < n_poly_points-1; ++i)
        shifted_point_list[i+1] = point_list[i];

      for (unsigned int i = 0; i < n_poly_points; ++i)
        {
          const std::array<Point<2>,2 > list = {{point_list[i], shifted_point_list[i]}};
          distances[i] = distance_to_line(list, point);
        }

      // Return the minimum of the distances of the point to all polygon segments
      return *std::min_element(distances.begin(),distances.end()) * sign;
    }

    double
    distance_to_line(const std::array<dealii::Point<2>,2 > &point_list,
                     const dealii::Point<2> &point)
    {

      /**
       * This code is based on http://geomalgorithms.com/a02-_lines.html#Distance-to-Infinite-Line,
       * and therefore requires the following copyright notice:
       *
       * Copyright 2000 softSurfer, 2012 Dan Sunday
       * This code may be freely used and modified for any purpose
       * providing that this copyright notice is included with it.
       * SoftSurfer makes no warranty for this code, and cannot be held
       * liable for any real or imagined damage resulting from its use.
       * Users of this code must verify correctness for their application.
       *
       */

      const unsigned int n_poly_points = point_list.size();
      AssertThrow(n_poly_points == 2, ExcMessage("A list of points for a line segment should consist of 2 points."));

      // Create vector along the polygon line segment P0 to P1
      const Tensor<1,2> vector_segment = point_list[1] - point_list[0];
      // Create vector from point P to the second segment point
      const Tensor<1,2> vector_point_segment = point - point_list[0];

      // Compute dot products to get angles
      const double c1 = vector_point_segment * vector_segment;

      // Point P's perpendicular base line lies outside segment, before P0.
      // Return distance between points P and P0.
      if (c1 <= 0.0)
        return (Tensor<1,2> (point_list[0] - point)).norm();

      const double c2 = vector_segment * vector_segment;

      // Point P's perpendicular base line lies outside segment, after P1.
      // Return distance between points P and P1.
      if (c2 <= c1)
        return (Tensor<1,2> (point_list[1] - point)).norm();

      // Point P's perpendicular base line lies on the line segment.
      // Return distance between point P and the base point.
      const Point<2> point_on_segment = point_list[0] + (c1/c2) * vector_segment;
      return (Tensor<1,2> (point - point_on_segment)).norm();
    }

    template <int dim>
    std::array<Tensor<1,dim>,dim-1>
    orthogonal_vectors (const Tensor<1,dim> &v)
    {
      Assert (v.norm() > 0,
              ExcMessage ("This function can not be called with a zero "
                          "input vector."));

      std::array<Tensor<1,dim>,dim-1> return_value;
      switch (dim)
        {
          case 2:
          {
            // create a direction by swapping the two coordinates and
            // flipping one sign; this is orthogonal to 'v' and has
            // the same length already
            return_value[0][0] = v[1];
            return_value[0][1] = -v[0];
            break;
          }

          case 3:
          {
            // In 3d, we can get two other vectors in a 3-step procedure:
            // - compute a 'w' that is definitely not collinear with 'v'
            // - compute the first direction as u[1] = v \times w,
            //   normalize it
            // - compute the second direction as u[2] = v \times u[1],
            //   normalize it
            //
            // For the first step, use a procedure suggested by Luca Heltai:
            // Set d to the index of the largest component of v. Make a vector
            // with this component equal to zero, and with the other two
            // equal to the norm of the point. Call this 'w'.
            unsigned int max_component = 0;
            for (unsigned int d=1; d<dim; ++d)
              if (std::fabs(v[d]) > std::fabs(v[max_component]))
                max_component = d;
            Tensor<1,dim> w = v;
            w[max_component] = 0;
            for (unsigned int d=1; d<dim; ++d)
              if (d != max_component)
                w[d] = v.norm();

            return_value[0] = cross_product_3d(v, w);
            return_value[0] *= v.norm() / return_value[0].norm();

            return_value[1] = cross_product_3d(v, return_value[0]);
            return_value[1] *= v.norm() / return_value[1].norm();

            break;
          }

          default:
            Assert (false, ExcNotImplemented());
        }

      return return_value;
    }


    std::pair<double,double> real_spherical_harmonic( const unsigned int l,
                                                      const unsigned int m,
                                                      const double theta,
                                                      const double phi)
    {
      const double sqrt_2 = numbers::SQRT2;
      const std::complex<double> sph_harm_val = boost::math::spherical_harmonic( l, m, theta, phi );
      if ( m == 0 )
        return std::make_pair( sph_harm_val.real(), 0.0 );
      else
        return std::make_pair( sqrt_2 * sph_harm_val.real(), sqrt_2 * sph_harm_val.imag() );
    }


    bool
    fexists(const std::string &filename)
    {
      std::ifstream ifile(filename.c_str());

      // return whether construction of the input file has succeeded;
      // success requires the file to exist and to be readable
      return static_cast<bool>(ifile);
    }


    bool
    filename_is_url(const std::string &filename)
    {
      if (filename.find("http://") == 0 || filename.find("https://") == 0 || filename.find("file://") == 0)
        return true;
      else
        return false;
    }


    std::string
    read_and_distribute_file_content(const std::string &filename,
                                     const MPI_Comm &comm)
    {
      std::string data_string;

      if (Utilities::MPI::this_mpi_process(comm) == 0)
        {
          // set file size to an invalid size (signaling an error if we can not read it)
          unsigned int filesize = numbers::invalid_unsigned_int;

          // Check to see if the prm file will be reading data from disk or
          // from a provided URL
          if (filename_is_url(filename))
            {
#ifdef ASPECT_WITH_LIBDAP
              libdap::Connect *url = new libdap::Connect(filename);
              libdap::BaseTypeFactory factory;
              libdap::DataDDS dds(&factory);
              libdap::DAS das;

              url->request_data(dds, "");
              url->request_das(das);


              // Temporary vector that will hold the different arrays stored in urlArray
              std::vector<std::string> tmp;
              // Vector that will hold the arrays (columns) and the values within those arrays
              std::vector<std::vector<std::string>> columns;

              // Check dds values to make sure the arrays are of the same length and of type string
              for (libdap::DDS::Vars_iter i = dds.var_begin(); i != dds.var_end(); ++i)
                {
                  libdap::BaseType *btp = *i;
                  if ((*i)->type() == libdap::dods_array_c)
                    {
                      // Array to store the url data
                      libdap::Array *urlArray;
                      urlArray = static_cast <libdap::Array *>(btp);
                      if (urlArray->var() != nullptr && urlArray->var()->type() == libdap::dods_str_c)
                        {
                          // The url Array contains a separate array for each column of data.
                          // This will put each of these individual arrays into its own vector.
                          urlArray->value(tmp);
                          columns.push_back(tmp);
                        }
                      else
                        {
                          AssertThrow (false,
                                       ExcMessage (std::string("Error when reading from url: ") + filename +
                                                   " Check your connection to the server and make sure the server "
                                                   "delivers correct data."));
                        }

                    }
                  else
                    {
                      AssertThrow (false,
                                   ExcMessage (std::string("Error when reading from url: ") + filename +
                                               " Check your connection to the server and make sure the server "
                                               "delivers correct data."));
                    }
                }

              // Add the POINTS data that is required and found at the top of the data file.
              // The POINTS values are set as attributes inside a table.
              // Loop through the Attribute table to locate the points values within
              std::vector<std::string> points;
              for (libdap::AttrTable::Attr_iter i = das.var_begin(); i != das.var_end(); i++)
                {
                  libdap::AttrTable *table;

                  table = das.get_table(i);
                  if (table->get_attr("POINTS") != "")
                    points.push_back(table->get_attr("POINTS"));
                }

              std::stringstream urlString;

              // Append the gathered POINTS in the proper format:
              // "# POINTS: <val1> <val2> <val3>"
              urlString << "# POINTS:";
              for (unsigned int i = 0; i < points.size(); i++)
                {
                  urlString << " " << points[i];
                }
              urlString << "\n";

              // Add the values from the arrays into the stringstream. The values are passed in
              // per row with a character return added at the end of each row.
              // TODO: Add a check to make sure that each column is the same size before writing
              //     to the stringstream
              for (unsigned int i = 0; i < tmp.size(); i++)
                {
                  for (unsigned int j = 0; j < columns.size(); j++)
                    {
                      urlString << columns[j][i];
                      urlString << " ";
                    }
                  urlString << "\n";
                }

              data_string = urlString.str();
              filesize = data_string.size();
              delete url;
#else // ASPECT_WITH_LIBDAP

              // broadcast failure state, then throw
              const int ierr = MPI_Bcast(&filesize, 1, MPI_UNSIGNED, 0, comm);
              AssertThrowMPI(ierr);
              AssertThrow(false,
                          ExcMessage(std::string("Reading of file ") + filename + " failed. " +
                                     "Make sure you have the dependencies for reading a url " +
                                     "(run cmake with -DASPECT_WITH_LIBDAP=ON)"));

#endif // ASPECT_WITH_LIBDAP
            }
          else
            {
              std::ifstream filestream(filename.c_str());

              if (!filestream)
                {
                  // broadcast failure state, then throw
                  const int ierr = MPI_Bcast(&filesize,1,MPI_UNSIGNED,0,comm);
                  AssertThrowMPI(ierr);
                  AssertThrow (false,
                               ExcMessage (std::string("Could not open file <") + filename + ">."));
                  return data_string; // never reached
                }


              // Read data from disk
              std::stringstream datastream;
              filestream >> datastream.rdbuf();

              if (!filestream.eof())
                {
                  // broadcast failure state, then throw
                  const int ierr = MPI_Bcast(&filesize,1,MPI_UNSIGNED,0,comm);
                  AssertThrowMPI(ierr);
                  AssertThrow (false,
                               ExcMessage (std::string("Reading of file ") + filename + " finished " +
                                           "before the end of file was reached. Is the file corrupted or"
                                           "too large for the input buffer?"));
                  return data_string; // never reached
                }

              data_string = datastream.str();
              filesize = data_string.size();
            }

          // Distribute data_size and data across processes
          int ierr = MPI_Bcast(&filesize,1,MPI_UNSIGNED,0,comm);
          AssertThrowMPI(ierr);
          ierr = MPI_Bcast(&data_string[0],filesize,MPI_CHAR,0,comm);
          AssertThrowMPI(ierr);
        }
      else
        {
          // Prepare for receiving data
          unsigned int filesize;
          int ierr = MPI_Bcast(&filesize,1,MPI_UNSIGNED,0,comm);
          AssertThrowMPI(ierr);
          if (filesize == numbers::invalid_unsigned_int)
            throw QuietException();

          data_string.resize(filesize);

          // Receive and store data
          ierr = MPI_Bcast(&data_string[0],filesize,MPI_CHAR,0,comm);
          AssertThrowMPI(ierr);
        }

      return data_string;
    }

    int
    mkdirp(std::string pathname,const mode_t mode)
    {
      // force trailing / so we can handle everything in loop
      if (pathname[pathname.size()-1] != '/')
        {
          pathname += '/';
        }

      size_t pre = 0;
      size_t pos;

      while ((pos = pathname.find_first_of('/',pre)) != std::string::npos)
        {
          const std::string subdir = pathname.substr(0,pos++);
          pre = pos;

          // if leading '/', first string is 0 length
          if (subdir.size() == 0)
            continue;

          int mkdir_return_value;
          if ((mkdir_return_value = mkdir(subdir.c_str(),mode)) && (errno != EEXIST))
            return mkdir_return_value;

        }

      return 0;
    }

    void create_directory(const std::string &pathname,
                          const MPI_Comm &comm,
                          bool silent)
    {
      // verify that the output directory actually exists. if it doesn't, create
      // it on processor zero
      int error;

      if ((Utilities::MPI::this_mpi_process(comm) == 0))
        {
          DIR *output_directory = opendir(pathname.c_str());
          if (output_directory == nullptr)
            {
              if (!silent)
                std::cout << "\n"
                          << "-----------------------------------------------------------------------------\n"
                          << "The output directory <" << pathname
                          << "> provided in the input file appears not to exist.\n"
                          << "ASPECT will create it for you.\n"
                          << "-----------------------------------------------------------------------------\n\n"
                          << std::endl;

              error = Utilities::mkdirp(pathname, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);

            }
          else
            {
              error = closedir(output_directory);
            }
          // Broadcast error code
          const int ierr = MPI_Bcast (&error, 1, MPI_INT, 0, comm);
          AssertThrowMPI(ierr);
          AssertThrow (error == 0,
                       ExcMessage (std::string("Can't create the output directory at <") + pathname + ">"));
        }
      else
        {
          // Wait to receive error code, and throw QuietException if directory
          // creation has failed
          const int ierr = MPI_Bcast (&error, 1, MPI_INT, 0, comm);
          AssertThrowMPI(ierr);
          if (error!=0)
            throw aspect::QuietException();
        }
    }

// tk does the cubic spline interpolation that can be used between different spherical layers in the mantle.
// This interpolation is based on the script spline.h, which was downloaded from
// http://kluge.in-chemnitz.de/opensource/spline/spline.h   //
// copyright (C) 2011, 2014 Tino Kluge (ttk448 at gmail.com)
    namespace tk
    {
      /**
       * Band matrix solver for a banded, square matrix of size @p dim. Used
       * by spline interpolation in tk::spline.
       */
      class band_matrix
      {
        public:
          /**
           * Constructor, see resize()
           */
          band_matrix(int dim, int n_u, int n_l);
          /**
           * Resize to a @p dim by @dim matrix with given number
           * of off-diagonals.
           */
          void resize(int dim, int n_u, int n_l);

          /**
           * Return the dimension of the matrix
           */
          int dim() const;

          /**
           * Number of off-diagonals above.
           */
          int num_upper() const
          {
            return m_upper.size()-1;
          }

          /**
           * Number of off-diagonals below.
           */
          int num_lower() const
          {
            return m_lower.size()-1;
          }

          /**
           * Writeable access to element A(i,j), indices going from
           * i=0,...,dim()-1
           */
          double &operator () (int i, int j);
          /**
           * Read-only access
           */
          double operator () (int i, int j) const;

          /**
           * second diagonal (used in LU decomposition), saved in m_lower[0]
           */
          double &saved_diag(int i);

          /**
           * second diagonal (used in LU decomposition), saved in m_lower[0]
           */
          double saved_diag(int i) const;

          /**
           * LU-Decomposition of a band matrix
           */
          void lu_decompose();

          /**
           * solves Ux=y
           */
          std::vector<double> r_solve(const std::vector<double> &b) const;

          /**
           * solves Ly=b
           */
          std::vector<double> l_solve(const std::vector<double> &b) const;

          /**
           * Solve Ax=b and builds LU decomposition using lu_decompose()
           * if @p is_lu_decomposed is false.
           */
          std::vector<double> lu_solve(const std::vector<double> &b,
                                       bool is_lu_decomposed=false);
        private:
          /**
           * diagonal and off-diagonals above
           */
          std::vector< std::vector<double> > m_upper;
          /**
           * diagonals below the diagonal
           */
          std::vector< std::vector<double> > m_lower;
      };

      band_matrix::band_matrix(int dim, int n_u, int n_l)
      {
        resize(dim, n_u, n_l);
      }

      void band_matrix::resize(int dim, int n_u, int n_l)
      {
        assert(dim > 0);
        assert(n_u >= 0);
        assert(n_l >= 0);
        m_upper.resize(n_u+1);
        m_lower.resize(n_l+1);
        for (size_t i=0; i<m_upper.size(); i++)
          {
            m_upper[i].resize(dim);
          }
        for (size_t i=0; i<m_lower.size(); i++)
          {
            m_lower[i].resize(dim);
          }
      }

      int band_matrix::dim() const
      {
        if (m_upper.size()>0)
          {
            return m_upper[0].size();
          }
        else
          {
            return 0;
          }
      }

      double &band_matrix::operator () (int i, int j)
      {
        int k = j - i;       // what band is the entry
        assert( (i >= 0) && (i<dim()) && (j >= 0) && (j < dim()) );
        assert( (-num_lower() <= k) && (k <= num_upper()) );
        // k=0 -> diagonal, k<0 lower left part, k>0 upper right part
        if (k >= 0)
          return m_upper[k][i];
        else
          return m_lower[-k][i];
      }

      double band_matrix::operator () (int i, int j) const
      {
        int k=j-i;       // what band is the entry
        assert( (i >= 0) && (i < dim()) && (j >= 0) && (j < dim()) );
        assert( (-num_lower() <= k) && (k <= num_upper()) );
        // k=0 -> diagonal, k<0 lower left part, k>0 upper right part
        if (k >= 0)
          return m_upper[k][i];
        else
          return m_lower[-k][i];
      }

      double band_matrix::saved_diag(int i) const
      {
        assert( (i >= 0) && (i < dim()) );
        return m_lower[0][i];
      }

      double &band_matrix::saved_diag(int i)
      {
        assert( (i >= 0) && (i < dim()) );
        return m_lower[0][i];
      }

      void band_matrix::lu_decompose()
      {
        int i_max,j_max;
        int j_min;
        double x;

        // preconditioning
        // normalize column i so that a_ii=1
        for (int i = 0; i < this->dim(); i++)
          {
            assert(this->operator()(i,i) != 0.0);
            this->saved_diag(i) = 1.0/this->operator()(i,i);
            j_min = std::max(0,i-this->num_lower());
            j_max = std::min(this->dim()-1,i+this->num_upper());
            for (int j = j_min; j <= j_max; j++)
              {
                this->operator()(i,j) *= this->saved_diag(i);
              }
            this->operator()(i,i) = 1.0;          // prevents rounding errors
          }

        // Gauss LR-Decomposition
        for (int k = 0; k < this->dim(); k++)
          {
            i_max = std::min(this->dim()-1,k+this->num_lower());  // num_lower not a mistake!
            for (int i = k+1; i <= i_max; i++)
              {
                assert(this->operator()(k,k) != 0.0);
                x = -this->operator()(i,k)/this->operator()(k,k);
                this->operator()(i,k) = -x;                         // assembly part of L
                j_max = std::min(this->dim()-1, k + this->num_upper());
                for (int j = k+1; j <= j_max; j++)
                  {
                    // assembly part of R
                    this->operator()(i,j) = this->operator()(i,j)+x*this->operator()(k,j);
                  }
              }
          }
      }

      std::vector<double> band_matrix::l_solve(const std::vector<double> &b) const
      {
        assert( this->dim() == static_cast<int>(b.size()) );
        std::vector<double> x(this->dim());
        int j_start;
        double sum;
        for (int i = 0; i < this->dim(); i++)
          {
            sum = 0;
            j_start = std::max(0,i-this->num_lower());
            for (int j = j_start; j < i; j++) sum += this->operator()(i,j)*x[j];
            x[i] = (b[i]*this->saved_diag(i)) - sum;
          }
        return x;
      }


      std::vector<double> band_matrix::r_solve(const std::vector<double> &b) const
      {
        assert( this->dim() == static_cast<int>(b.size()) );
        std::vector<double> x(this->dim());
        int j_stop;
        double sum;
        for (int i = this->dim()-1; i >= 0; i--)
          {
            sum = 0;
            j_stop = std::min(this->dim()-1, i + this->num_upper());
            for (int j = i+1; j <= j_stop; j++) sum += this->operator()(i,j)*x[j];
            x[i] = (b[i] - sum) / this->operator()(i,i);
          }
        return x;
      }

      std::vector<double> band_matrix::lu_solve(const std::vector<double> &b,
                                                bool is_lu_decomposed)
      {
        assert( this->dim() == static_cast<int>(b.size()) );
        std::vector<double>  x,y;
        // TODO: this is completely unsafe because you rely on the user
        // if the function is called more than once.
        if (is_lu_decomposed == false)
          {
            this->lu_decompose();
          }
        y = this->l_solve(b);
        x = this->r_solve(y);
        return x;
      }


      void spline::set_points(const std::vector<double> &x,
                              const std::vector<double> &y,
                              bool cubic_spline,
                              bool monotone_spline)
      {
        assert(x.size() == y.size());
        m_x = x;
        m_y = y;
        const unsigned int n = x.size();
        for (unsigned int i = 0; i < n-1; i++)
          {
            assert(m_x[i] < m_x[i+1]);
          }

        if (cubic_spline == true)  // cubic spline interpolation
          {
            if (monotone_spline == true)
              {
                /**
                 * This monotone spline algorithm is based on the javascript version
                 * at https://en.wikipedia.org/wiki/Monotone_cubic_interpolation. The
                 * parameters from this algorithm prevent overshooting in the
                 * interpolation spline.
                 */
                std::vector<double> dys(n-1), dxs(n-1), ms(n-1);
                for (unsigned int i=0; i < n-1; i++)
                  {
                    dxs[i] = x[i+1]-x[i];
                    dys[i] = y[i+1]-y[i];
                    ms[i] = dys[i]/dxs[i];
                  }

                // get m_a parameter
                m_c.resize(n);
                m_c[0] = 0;

                for (unsigned int i = 0; i < n-2; i++)
                  {
                    const double m0 = ms[i];
                    const double m1 = ms[i+1];

                    if (m0 * m1 <= 0)
                      {
                        m_c[i+1] = 0;
                      }
                    else
                      {
                        const double dx0 = dxs[i];
                        const double dx1 = dxs[i+1];
                        const double common = dx0 + dx1;
                        m_c[i+1] = 3*common/((common + dx0)/m0 + (common + dx1)/m1);
                      }
                  }
                m_c[n-1] = ms[n-2];

                // Get b and c coefficients
                m_a.resize(n);
                m_b.resize(n);
                for (unsigned int i = 0; i < m_c.size()-1; i++)
                  {
                    const double c1 = m_c[i];
                    const double m0 = ms[i];

                    const double invDx = 1/dxs[i];
                    const double common0 = c1 + m_c[i+1] - m0 - m0;
                    m_b[i] = (m0 - c1 - common0) * invDx;
                    m_a[i] = common0 * invDx * invDx;
                  }
              }
            else
              {
                // setting up the matrix and right hand side of the equation system
                // for the parameters b[]
                band_matrix A(n,1,1);
                std::vector<double>  rhs(n);
                for (unsigned int i = 1; i<n-1; i++)
                  {
                    A(i,i-1) = 1.0/3.0*(x[i]-x[i-1]);
                    A(i,i) = 2.0/3.0*(x[i+1]-x[i-1]);
                    A(i,i+1) = 1.0/3.0*(x[i+1]-x[i]);
                    rhs[i] = (y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
                  }
                // boundary conditions, zero curvature b[0]=b[n-1]=0
                A(0,0) = 2.0;
                A(0,1) = 0.0;
                rhs[0] = 0.0;
                A(n-1,n-1) = 2.0;
                A(n-1,n-2) = 0.0;
                rhs[n-1] = 0.0;

                // solve the equation system to obtain the parameters b[]
                m_b = A.lu_solve(rhs);

                // calculate parameters a[] and c[] based on b[]
                m_a.resize(n);
                m_c.resize(n);
                for (unsigned int i = 0; i<n-1; i++)
                  {
                    m_a[i] = 1.0/3.0*(m_b[i+1]-m_b[i])/(x[i+1]-x[i]);
                    m_c[i] = (y[i+1]-y[i])/(x[i+1]-x[i])
                             - 1.0/3.0*(2.0*m_b[i]+m_b[i+1])*(x[i+1]-x[i]);
                  }
              }
          }
        else     // linear interpolation
          {
            m_a.resize(n);
            m_b.resize(n);
            m_c.resize(n);
            for (unsigned int i = 0; i<n-1; i++)
              {
                m_a[i] = 0.0;
                m_b[i] = 0.0;
                m_c[i] = (m_y[i+1]-m_y[i])/(m_x[i+1]-m_x[i]);
              }
          }

        // for the right boundary we define
        // f_{n-1}(x) = b*(x-x_{n-1})^2 + c*(x-x_{n-1}) + y_{n-1}
        double h = x[n-1]-x[n-2];
        // m_b[n-1] is determined by the boundary condition
        if (!monotone_spline)
          {
            m_a[n-1] = 0.0;
            m_c[n-1] = 3.0*m_a[n-2]*h*h+2.0*m_b[n-2]*h+m_c[n-2];   // = f'_{n-2}(x_{n-1})
          }
      }

      double spline::operator() (double x) const
      {
        size_t n = m_x.size();
        // find the closest point m_x[idx] < x, idx=0 even if x<m_x[0]
        std::vector<double>::const_iterator it;
        it = std::lower_bound(m_x.begin(),m_x.end(),x);
        int idx = std::max( int(it-m_x.begin())-1, 0);

        double h = x-m_x[idx];
        double interpol;
        if (x<m_x[0])
          {
            // extrapolation to the left
            interpol = ((m_b[0])*h + m_c[0])*h + m_y[0];
          }
        else if (x>m_x[n-1])
          {
            // extrapolation to the right
            interpol = ((m_b[n-1])*h + m_c[n-1])*h + m_y[n-1];
          }
        else
          {
            // interpolation
            interpol = ((m_a[idx]*h + m_b[idx])*h + m_c[idx])*h + m_y[idx];
          }
        return interpol;
      }

    } // namespace tk


    std::string
    expand_ASPECT_SOURCE_DIR (const std::string &location)
    {
      return Utilities::replace_in_string(location,
                                          "$ASPECT_SOURCE_DIR",
                                          ASPECT_SOURCE_DIR);
    }

    std::string parenthesize_if_nonempty (const std::string &s)
    {
      if (s.size() > 0)
        return " (\"" + s + "\")";
      else
        return "";
    }

    bool
    has_unique_entries (const std::vector<std::string> &strings)
    {
      const std::set<std::string> set_of_strings(strings.begin(),strings.end());
      return (set_of_strings.size() == strings.size());
    }

    template <int dim>
    AsciiDataLookup<dim>::AsciiDataLookup(const unsigned int components,
                                          const double scale_factor)
      :
      components(components),
      data(components),
      maximum_component_value(components),
      scale_factor(scale_factor),
      coordinate_values_are_equidistant(false)
    {}



    template <int dim>
    AsciiDataLookup<dim>::AsciiDataLookup(const double scale_factor)
      :
      components(numbers::invalid_unsigned_int),
      data(),
      maximum_component_value(),
      scale_factor(scale_factor),
      coordinate_values_are_equidistant(false)
    {}



    template <int dim>
    std::vector<std::string>
    AsciiDataLookup<dim>::get_column_names() const
    {
      return data_component_names;
    }



    template <int dim>
    bool
    AsciiDataLookup<dim>::has_equidistant_coordinates() const
    {
      return coordinate_values_are_equidistant;
    }



    template <int dim>
    const std::vector<double> &
    AsciiDataLookup<dim>::get_coordinates(const unsigned int dimension) const
    {
      AssertThrow(dimension < dim,
                  ExcMessage("There is no spatial dimension number " + std::to_string(dimension)
                             + " in the current data file."));

      return coordinate_values[dimension];
    }



    template <int dim>
    unsigned int
    AsciiDataLookup<dim>::get_column_index_from_name(const std::string &column_name) const
    {
      const std::vector<std::string>::const_iterator column_position =
        std::find(data_component_names.begin(),data_component_names.end(),column_name);

      AssertThrow(column_position != data_component_names.end(),
                  ExcMessage("There is no data column named " + column_name
                             + " in the current data file. Please check the name and the "
                             "first line not starting with '#' of your data file."));

      return std::distance(data_component_names.begin(),column_position);
    }

    template <int dim>
    std::string
    AsciiDataLookup<dim>::get_column_name_from_index(const unsigned int column_index) const
    {
      AssertThrow(data_component_names.size() > column_index,
                  ExcMessage("There is no data column number " + Utilities::to_string(column_index)
                             + " in the current data file. The data file only contains "
                             + Utilities::to_string(data_component_names.size()) + " named columns."));

      return data_component_names[column_index];
    }

    template <int dim>
    double
    AsciiDataLookup<dim>::get_maximum_component_value(const unsigned int component) const
    {
      return maximum_component_value[component];
    }



    template <int dim>
    void
    AsciiDataLookup<dim>::reinit(const std::vector<std::string> &column_names,
                                 const std::vector<std::vector<double>> &coordinate_values_,
                                 const std::vector<Table<dim,double> > &raw_data)
    {
      Assert(coordinate_values_.size()==dim, ExcMessage("Invalid size of coordinate_values."));
      for (unsigned int d=0; d<dim; ++d)
        {
          coordinate_values[d] = coordinate_values_[d];
          AssertThrow(coordinate_values[d].size()>1,
                      ExcMessage("Error: At least 2 entries per coordinate direction are required."));
          table_points[d] = coordinate_values_[d].size();
        }

      components = column_names.size();
      data_component_names = column_names;
      Assert(raw_data.size() == components,
             ExcMessage("Error: Incorrect number of columns specified."));

      // compute maximum_component_value for each component:
      maximum_component_value = std::vector<double>(components,-std::numeric_limits<double>::max());
      for (unsigned int c=0; c<components; ++c)
        {
          Assert(raw_data[c].size() == table_points,
                 ExcMessage("Error: One of the data tables has an incorrect size."));

          const unsigned int n_elements = raw_data[c].n_elements();
          for (unsigned int idx=0; idx<n_elements; ++idx)
            maximum_component_value[c] = std::max(maximum_component_value[c], raw_data[c](
                                                    compute_table_indices(table_points, idx)));
        }

      // In case the data is specified on a grid that is equidistant
      // in each coordinate direction, we only need to store
      // (besides the data) the number of intervals in each direction and
      // the begin- and endpoints of the coordinates.
      // In case the grid is not equidistant, we need to keep
      // all the coordinates in each direction, which is more costly.
      std::array<unsigned int,dim> table_intervals;

      coordinate_values_are_equidistant = true;
      for (unsigned int d=0; d<dim; ++d)
        {
          table_intervals[d] = table_points[d]-1;

          // The minimum and maximum coordinate values:
          grid_extent[d].first = coordinate_values[d][0];
          grid_extent[d].second = coordinate_values[d][table_points[d]-1];

          const double grid_spacing = coordinate_values[d][1] - coordinate_values[d][0];

          for (unsigned int n = 1; n < table_points[d]; ++n)
            {
              const double current_grid_spacing = coordinate_values[d][n] - coordinate_values_[d][n-1];

              AssertThrow(current_grid_spacing > 0.0,
                          ExcMessage ("Coordinates in dimension "
                                      + Utilities::int_to_string(d)
                                      + " are not strictly ascending."));

              // If spacing between coordinates changed (with a relative
              // tolerance), keep track of that information.  Note that we do
              // not break out of this loop in this case but run through the
              // whole array, so that the AssertThrow above is executed for
              // each entry to ensure increasing coordinate values.
              if (std::abs(current_grid_spacing - grid_spacing) > 0.005*(current_grid_spacing+grid_spacing))
                coordinate_values_are_equidistant = false;
            }
        }

      // For each data component, set up a GridData,
      // its type depending on the read-in grid.
      data.resize(components);
      for (unsigned int c = 0; c < components; ++c)
        {
          if (coordinate_values_are_equidistant)
            data[c]
              = std_cxx14::make_unique<Functions::InterpolatedUniformGridData<dim>> (grid_extent,
                                                                                     table_intervals,
                                                                                     raw_data[c]);
          else
            data[c]
              = std_cxx14::make_unique<Functions::InterpolatedTensorProductGridData<dim>> (coordinate_values,
                                                                                           raw_data[c]);
        }
    }



    template <int dim>
    void
    AsciiDataLookup<dim>::load_file(const std::string &filename,
                                    const MPI_Comm &comm)
    {
      // Grab the values already stored in this class (if they exist), this way we can
      // check if somebody changes the size of the table over time and error out (see below)
      TableIndices<dim> new_table_points = this->table_points;
      std::vector<std::string> column_names;

      // Read data from disk and distribute among processes
      std::stringstream in(read_and_distribute_file_content(filename, comm));

      // Read header lines and table size
      while (in.peek() == '#')
        {
          std::string line;
          std::getline(in,line);
          std::stringstream linestream(line);
          std::string word;
          while (linestream >> word)
            if (word == "POINTS:")
              for (unsigned int i = 0; i < dim; i++)
                {
                  unsigned int temp_index;
                  linestream >> temp_index;

                  if (new_table_points[i] == 0)
                    new_table_points[i] = temp_index;
                  else
                    AssertThrow (new_table_points[i] == temp_index,
                                 ExcMessage("The file grid must not change over model runtime. "
                                            "Either you prescribed a conflicting number of points in "
                                            "the input file, or the POINTS comment in your data files "
                                            "is changing between following files."));
                }
        }

      for (unsigned int i = 0; i < dim; i++)
        {
          AssertThrow(new_table_points[i] != 0,
                      ExcMessage("Could not successfully read in the file header of the "
                                 "ascii data file <" + filename + ">. One header line has to "
                                 "be of the format: '#POINTS: N1 [N2] [N3]', where N1 and "
                                 "potentially N2 and N3 have to be the number of data points "
                                 "in their respective dimension. Check for typos in this line "
                                 "(e.g. a missing space character)."));
        }

      // Read column lines if present
      unsigned int name_column_index = 0;
      double temp_data;

      while (true)
        {
          AssertThrow (name_column_index < 100,
                       ExcMessage("The program found more than 100 columns in the first line of the data file. "
                                  "This is unlikely intentional. Check your data file and make sure the data can be "
                                  "interpreted as floating point numbers. If you do want to read a data file with more "
                                  "than 100 columns, please remove this assertion."));

          std::string column_name_or_data;
          in >> column_name_or_data;
          try
            {
              // If the data field contains a name this will throw an exception
              temp_data = boost::lexical_cast<double>(column_name_or_data);

              // If there was no exception we have left the line containing names
              // and have read the first data field. Save number of components, and
              // make sure there is no contradiction if the components were already given to
              // the constructor of this class.
              if (components == numbers::invalid_unsigned_int)
                components = name_column_index - dim;
              else if (name_column_index != 0)
                AssertThrow (components == name_column_index,
                             ExcMessage("The number of expected data columns and the "
                                        "list of column names at the beginning of the data file "
                                        + filename + " do not match. The file should contain "
                                        "one column name per column (one for each dimension "
                                        "and one per data column)."));

              break;
            }
          catch (const boost::bad_lexical_cast &e)
            {
              // The first dim columns are coordinates and contain no data
              if (name_column_index >= dim)
                {
                  // Transform name to lower case to prevent confusion with capital letters
                  // Note: only ASCII characters allowed
                  std::transform(column_name_or_data.begin(), column_name_or_data.end(), column_name_or_data.begin(), ::tolower);

                  AssertThrow(std::find(column_names.begin(),column_names.end(),column_name_or_data)
                              == column_names.end(),
                              ExcMessage("There are multiple fields named " + column_name_or_data +
                                         " in the data file " + filename + ". Please remove duplication to "
                                         "allow for unique association between column and name."));

                  column_names.push_back(column_name_or_data);
                }
              ++name_column_index;
            }
        }

      // Create table for the data. This peculiar reinit is necessary, because
      // there is no constructor for Table, which takes TableIndices as
      // argument.
      Table<dim,double> data_table;
      data_table.TableBase<dim,double>::reinit(new_table_points);
      std::vector<Table<dim,double> > data_tables(components, data_table);

      std::vector<std::vector<double>> coordinate_values(dim);
      for (unsigned int d=0; d<dim; ++d)
        coordinate_values[d].resize(new_table_points[d]);

      if (column_names.size()==0)
        {
          // set default column names:
          for (unsigned int c=0; c<components; ++c)
            column_names.push_back("column " + Utilities::int_to_string(c,2));
        }

      // Finally read data lines:
      unsigned int read_data_entries = 0;
      do
        {
          // what row and column of the file are we in?
          const unsigned int column_num = read_data_entries%(components+dim);
          const unsigned int row_num = read_data_entries/(components+dim);
          TableIndices<dim> idx = compute_table_indices(new_table_points, row_num);

          if (column_num < dim)
            {
              // This is a coordinate. Store (and check that they are consistent)
              const double old_value = coordinate_values[column_num][idx[column_num]];

              AssertThrow(old_value == 0. ||
                          (std::abs(old_value-temp_data) < 1e-8*std::abs(old_value)),
                          ExcMessage("Invalid coordinate "
                                     + Utilities::int_to_string(column_num) + " in row "
                                     + Utilities::int_to_string(row_num)
                                     + " in file " + filename +
                                     "\nThis class expects the coordinates to be structured, meaning "
                                     "the coordinate values in each coordinate direction repeat exactly "
                                     "each time."));

              coordinate_values[column_num][idx[column_num]] = temp_data;
            }
          else
            {
              // This is a data value, so scale and store:
              const unsigned int component = column_num - dim;
              data_tables[component](idx) = temp_data * scale_factor;
            }

          ++read_data_entries;
        }
      while (in >> temp_data);

      AssertThrow(in.eof(),
                  ExcMessage ("While reading the data file '" + filename + "' the ascii data "
                              "plugin has encountered an error before the end of the file. "
                              "Please check for malformed data values (e.g. NaN) or superfluous "
                              "lines at the end of the data file."));

      const unsigned int n_expected_data_entries = (components + dim) * data_table.n_elements();
      AssertThrow(read_data_entries == n_expected_data_entries,
                  ExcMessage ("While reading the data file '" + filename + "' the ascii data "
                              "plugin has reached the end of the file, but has not found the "
                              "expected number of data values considering the spatial dimension, "
                              "data columns, and number of lines prescribed by the POINTS header "
                              "of the file. Please check the number of data "
                              "lines against the POINTS header in the file."));

      // finally create the data:
      this->reinit(column_names, coordinate_values, data_tables);
    }


    template <int dim>
    double
    AsciiDataLookup<dim>::get_data(const Point<dim> &position,
                                   const unsigned int component) const
    {
      Assert(component<components, ExcMessage("Invalid component index"));
      return data[component]->value(position);
    }

    template <int dim>
    Tensor<1,dim>
    AsciiDataLookup<dim>::get_gradients(const Point<dim> &position,
                                        const unsigned int component)
    {
      return data[component]->gradient(position,0);
    }


    template <int dim>
    TableIndices<dim>
    AsciiDataLookup<dim>::compute_table_indices(const TableIndices<dim> &sizes, const unsigned int i) const
    {
      TableIndices<dim> idx;
      idx[0] = i % sizes[0];
      if (dim >= 2)
        idx[1] = (i / sizes[0]) % sizes[1];
      if (dim == 3)
        idx[2] = i / (sizes[0] * sizes[1]);

      return idx;
    }



    template <int dim>
    AsciiDataBase<dim>::AsciiDataBase ()
    {}


    template <int dim>
    void
    AsciiDataBase<dim>::declare_parameters (ParameterHandler  &prm,
                                            const std::string &default_directory,
                                            const std::string &default_filename,
                                            const std::string &subsection_name)
    {
      prm.enter_subsection (subsection_name);
      {
        prm.declare_entry ("Data directory",
                           default_directory,
                           Patterns::DirectoryName (),
                           "The name of a directory that contains the model data. This path "
                           "may either be absolute (if starting with a `/') or relative to "
                           "the current directory. The path may also include the special "
                           "text `$ASPECT_SOURCE_DIR' which will be interpreted as the path "
                           "in which the ASPECT source files were located when ASPECT was "
                           "compiled. This interpretation allows, for example, to reference "
                           "files located in the `data/' subdirectory of ASPECT.");
        prm.declare_entry ("Data file name",
                           default_filename,
                           Patterns::Anything (),
                           "The file name of the model data.");
        prm.declare_entry ("Scale factor", "1.",
                           Patterns::Double (),
                           "Scalar factor, which is applied to the model data. "
                           "You might want to use this to scale the input to a "
                           "reference model. Another way to use this factor is to "
                           "convert units of the input files. For instance, if you "
                           "provide velocities in cm/yr set this factor to 0.01.");
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    AsciiDataBase<dim>::parse_parameters (ParameterHandler &prm,
                                          const std::string &subsection_name)
    {
      prm.enter_subsection (subsection_name);
      {
        // Get the path to the data files. If it contains a reference
        // to $ASPECT_SOURCE_DIR, replace it by what CMake has given us
        // as a #define
        data_directory = Utilities::expand_ASPECT_SOURCE_DIR(prm.get ("Data directory"));
        data_file_name    = prm.get ("Data file name");
        scale_factor      = prm.get_double ("Scale factor");
      }
      prm.leave_subsection();
    }



    template <int dim>
    AsciiDataBoundary<dim>::AsciiDataBoundary ()
      :
      current_file_number(0),
      first_data_file_model_time(0.0),
      first_data_file_number(0),
      decreasing_file_order(false),
      data_file_time_step(0.0),
      time_weight(0.0),
      time_dependent(true),
      lookups(),
      old_lookups()
    {}



    template <int dim>
    void
    AsciiDataBoundary<dim>::initialize(const std::set<types::boundary_id> &boundary_ids,
                                       const unsigned int components)
    {
      AssertThrow ((Plugins::plugin_type_matches<const GeometryModel::SphericalShell<dim>> (this->get_geometry_model()))
                   || (Plugins::plugin_type_matches<const GeometryModel::Chunk<dim>> (this->get_geometry_model()))
                   || (Plugins::plugin_type_matches<const GeometryModel::Sphere<dim>> (this->get_geometry_model()))
                   || (Plugins::plugin_type_matches<const GeometryModel::Box<dim>> (this->get_geometry_model()))
                   || (Plugins::plugin_type_matches<const GeometryModel::TwoMergedBoxes<dim>> (this->get_geometry_model())),
                   ExcMessage ("This ascii data plugin can only be used when using "
                               "a spherical shell, chunk, box or two merged boxes geometry."));


      for (const auto &boundary_id : boundary_ids)
        {
          lookups.insert(std::make_pair(boundary_id,
                                        std_cxx14::make_unique<Utilities::AsciiDataLookup<dim-1>>
                                        (components,
                                         this->scale_factor)));

          old_lookups.insert(std::make_pair(boundary_id,
                                            std_cxx14::make_unique<Utilities::AsciiDataLookup<dim-1>>
                                            (components,
                                             this->scale_factor)));

          // Set the first file number and load the first files
          current_file_number = first_data_file_number;

          const int next_file_number =
            (decreasing_file_order) ?
            current_file_number - 1
            :
            current_file_number + 1;

          const std::string filename (create_filename (current_file_number, boundary_id));

          this->get_pcout() << std::endl << "   Loading Ascii data boundary file "
                            << filename << "." << std::endl << std::endl;


          AssertThrow(Utilities::fexists(filename) || filename_is_url(filename),
                      ExcMessage (std::string("Ascii data file <")
                                  +
                                  filename
                                  +
                                  "> not found!"));
          lookups.find(boundary_id)->second->load_file(filename,this->get_mpi_communicator());

          // If the boundary condition is constant, switch off time_dependence
          // immediately. If not, also load the second file for interpolation.
          // This catches the case that many files are present, but the
          // parameter file requests a single file.
          if (filename == create_filename (current_file_number+1, boundary_id))
            {
              end_time_dependence ();
            }
          else
            {
              const std::string filename (create_filename (next_file_number, boundary_id));
              this->get_pcout() << std::endl << "   Loading Ascii data boundary file "
                                << filename << "." << std::endl << std::endl;
              if (Utilities::fexists(filename))
                {
                  lookups.find(boundary_id)->second.swap(old_lookups.find(boundary_id)->second);
                  lookups.find(boundary_id)->second->load_file(filename, this->get_mpi_communicator());
                }
              else
                end_time_dependence ();
            }
        }
    }



    template <int dim>
    std::array<unsigned int,dim-1>
    AsciiDataBoundary<dim>::get_boundary_dimensions (const types::boundary_id boundary_id) const
    {
      std::array<unsigned int,dim-1> boundary_dimensions;

      switch (dim)
        {
          case 2:
            if ((boundary_id == 2) || (boundary_id == 3))
              {
                boundary_dimensions[0] = 0;
              }
            else if ((boundary_id == 0) || (boundary_id == 1))
              {
                boundary_dimensions[0] = 1;
              }
            else
              {
                boundary_dimensions[0] = numbers::invalid_unsigned_int;
                AssertThrow(false,ExcNotImplemented());
              }

            break;

          case 3:
            if ((boundary_id == 4) || (boundary_id == 5))
              {
                boundary_dimensions[0] = 0;
                boundary_dimensions[1] = 1;
              }
            else if ((boundary_id == 0) || (boundary_id == 1))
              {
                boundary_dimensions[0] = 1;
                boundary_dimensions[1] = 2;
              }
            else if ((boundary_id == 2) || (boundary_id == 3))
              {
                boundary_dimensions[0] = 0;
                boundary_dimensions[1] = 2;
              }
            else
              {
                boundary_dimensions[0] = numbers::invalid_unsigned_int;
                boundary_dimensions[1] = numbers::invalid_unsigned_int;
                AssertThrow(false,ExcNotImplemented());
              }

            break;

          default:
            for (unsigned int d=0; d<dim-1; ++d)
              boundary_dimensions[d] = numbers::invalid_unsigned_int;
            AssertThrow(false,ExcNotImplemented());
        }
      return boundary_dimensions;
    }

    namespace
    {
      /**
       * Given a string @p filename_and_path that contains exactly one
       * <code>%s</code> and one <code>%d</code> code (possibly modified
       * by flag, field, and length modifiers as discussed in the man
       * pages of the <code>printf()</code> family of functions),
       * return the expanded string where the <code>%s</code> code is
       * replaced by @p boundary_name, and <code>%d</code> is replaced
       * by @p filenumber.
       */
      std::string replace_placeholders(const std::string &filename_and_path,
                                       const std::string &boundary_name,
                                       const int filenumber)
      {
        const int maxsize = filename_and_path.length() + 256;
        char *filename = static_cast<char *>(malloc (maxsize * sizeof(char)));
        int ret = snprintf (filename,
                            maxsize,
                            filename_and_path.c_str(),
                            boundary_name.c_str(),
                            filenumber);

        AssertThrow(ret >= 0, ExcMessage("Invalid string placeholder in filename detected."));
        AssertThrow(ret< maxsize, ExcInternalError("snprintf string overflow detected."));
        const std::string str_result (filename);
        free (filename);
        return str_result;
      }

    }

    template <int dim>
    std::string
    AsciiDataBoundary<dim>::create_filename (const int filenumber,
                                             const types::boundary_id boundary_id) const
    {
      std::string templ = this->data_directory + this->data_file_name;

      const std::string boundary_name = this->get_geometry_model().translate_id_to_symbol_name(boundary_id);

      const std::string result = replace_placeholders(templ, boundary_name, filenumber);
      if (fexists(result))
        return result;

      // Backwards compatibility check: people might still be using the old
      // names of the top/bottom boundary. If they do, print a warning but
      // accept those files.
      std::string compatible_result;
      if (boundary_name == "top")
        {
          compatible_result = replace_placeholders(templ, "surface", filenumber);
          if (!fexists(compatible_result))
            compatible_result = replace_placeholders(templ, "outer", filenumber);
        }
      else if (boundary_name == "bottom")
        compatible_result = replace_placeholders(templ, "inner", filenumber);

      if (!fexists(result) && fexists(compatible_result))
        {
          this->get_pcout() << "WARNING: Filename convention concerning geometry boundary "
                            "names changed. Please rename '" << compatible_result << "'"
                            << " to '" << result << "'"
                            << std::endl;
          return compatible_result;
        }

      return result;
    }


    template <int dim>
    void
    AsciiDataBoundary<dim>::update ()
    {
      if (time_dependent && (this->get_time() - first_data_file_model_time >= 0.0))
        {
          const double time_steps_since_start = (this->get_time() - first_data_file_model_time)
                                                / data_file_time_step;
          // whether we need to update our data files. This looks so complicated
          // because we need to catch increasing and decreasing file orders and all
          // possible first_data_file_model_times and first_data_file_numbers.
          const bool need_update =
            static_cast<int> (time_steps_since_start)
            > std::abs(current_file_number - first_data_file_number);

          if (need_update)
            {
              // The last file, which was tried to be loaded was
              // number current_file_number +/- 1, because current_file_number
              // is the file older than the current model time
              const int old_file_number =
                (decreasing_file_order) ?
                current_file_number - 1
                :
                current_file_number + 1;

              // Calculate new file_number
              current_file_number =
                (decreasing_file_order) ?
                first_data_file_number
                - static_cast<unsigned int> (time_steps_since_start)
                :
                first_data_file_number
                + static_cast<unsigned int> (time_steps_since_start);

              const bool load_both_files = std::abs(current_file_number - old_file_number) >= 1;

              for (const auto &boundary_id : lookups)
                update_data(boundary_id.first, load_both_files);
            }

          time_weight = time_steps_since_start
                        - std::abs(current_file_number - first_data_file_number);

          Assert ((0 <= time_weight) && (time_weight <= 1),
                  ExcMessage (
                    "Error in set_current_time. Time_weight has to be in [0,1]"));
        }
    }

    template <int dim>
    void
    AsciiDataBoundary<dim>::update_data (const types::boundary_id boundary_id,
                                         const bool load_both_files)
    {
      // If the time step was large enough to move forward more
      // then one data file we need to load both current files
      // to stay accurate in interpolation
      if (load_both_files)
        {
          const std::string filename (create_filename (current_file_number,boundary_id));
          this->get_pcout() << std::endl << "   Loading Ascii data boundary file "
                            << filename << "." << std::endl << std::endl;
          if (Utilities::fexists(filename))
            {
              lookups.find(boundary_id)->second.swap(old_lookups.find(boundary_id)->second);
              lookups.find(boundary_id)->second->load_file(filename,this->get_mpi_communicator());
            }

          // If loading current_time_step failed, end time dependent part with old_file_number.
          else
            end_time_dependence ();
        }

      // Now load the next data file. This part is the main purpose of this function.
      const int next_file_number =
        (decreasing_file_order) ?
        current_file_number - 1
        :
        current_file_number + 1;

      const std::string filename (create_filename (next_file_number,boundary_id));
      this->get_pcout() << std::endl << "   Loading Ascii data boundary file "
                        << filename << "." << std::endl << std::endl;
      if (Utilities::fexists(filename))
        {
          lookups.find(boundary_id)->second.swap(old_lookups.find(boundary_id)->second);
          lookups.find(boundary_id)->second->load_file(filename,this->get_mpi_communicator());
        }

      // If next file does not exist, end time dependent part with current_time_step.
      else
        end_time_dependence ();
    }

    template <int dim>
    void
    AsciiDataBoundary<dim>::end_time_dependence ()
    {
      // no longer consider the problem time dependent from here on out
      // this cancels all attempts to read files at the next time steps
      time_dependent = false;
      // Give warning if first processor
      this->get_pcout() << std::endl
                        << "   From this timestep onwards, ASPECT will not attempt to load new Ascii data files." << std::endl
                        << "   This is either because ASPECT has already read all the files necessary to impose" << std::endl
                        << "   the requested boundary condition, or that the last available file has been read." << std::endl
                        << "   If the Ascii data represented a time-dependent boundary condition," << std::endl
                        << "   that time-dependence ends at this timestep  (i.e. the boundary condition" << std::endl
                        << "   will continue unchanged from the last known state into the future)." << std::endl << std::endl;
    }


    template <int dim>
    double
    AsciiDataBoundary<dim>::
    get_data_component (const types::boundary_id             boundary_indicator,
                        const Point<dim>                    &position,
                        const unsigned int                   component) const
    {
      // For initial ascii data topography, we need access to the data before get_time() is set,
      // as this is when the grid including topography is constructed for the chunk geometry.
      if ( (dynamic_cast<const GeometryModel::Chunk<dim>*>(&this->get_geometry_model()) != nullptr &&
            dynamic_cast<const InitialTopographyModel::AsciiData<dim>*>(&this->get_initial_topography_model()) != nullptr &&
            this->get_timestep_number() == numbers::invalid_unsigned_int) ||
           this->get_time() - first_data_file_model_time >= 0.0)
        {
          const std::array<double,dim> natural_position = this->get_geometry_model().cartesian_to_natural_coordinates(position);

          Point<dim> internal_position;
          for (unsigned int i = 0; i < dim; i++)
            internal_position[i] = natural_position[i];

          // The chunk model has latitude as natural coordinate. We need to convert this to colatitude
          if (Plugins::plugin_type_matches<const GeometryModel::Chunk<dim>> (this->get_geometry_model()) && dim == 3)
            {
              internal_position[2] = numbers::PI/2. - internal_position[2];
            }

          const std::array<unsigned int,dim-1> boundary_dimensions =
            get_boundary_dimensions(boundary_indicator);

          Point<dim-1> data_position;
          for (unsigned int i = 0; i < dim-1; i++)
            data_position[i] = internal_position[boundary_dimensions[i]];

          const double data = lookups.find(boundary_indicator)->second->get_data(data_position,component);

          if (!time_dependent)
            return data;

          const double old_data = old_lookups.find(boundary_indicator)->second->get_data(data_position,component);

          return time_weight * data + (1 - time_weight) * old_data;
        }
      else
        return 0.0;
    }


    template <int dim>
    Tensor<1,dim-1>
    AsciiDataBoundary<dim>::vector_gradient (const types::boundary_id             boundary_indicator,
                                             const Point<dim>                    &position,
                                             const unsigned int                   component) const
    {
      // For initial ascii data topography, we need access to the data before get_time() is set,
      // as this is when the grid including topography is constructed for the chunk geometry.
      if ((dynamic_cast<const GeometryModel::Chunk<dim>*>(&this->get_geometry_model()) != nullptr &&
           dynamic_cast<const InitialTopographyModel::AsciiData<dim>*>(&this->get_initial_topography_model()) != nullptr &&
           this->get_timestep_number() == numbers::invalid_unsigned_int) ||
          this->get_time() - first_data_file_model_time >= 0.0 )
        {
          const std::array<double,dim> natural_position = this->get_geometry_model().cartesian_to_natural_coordinates(position);

          Point<dim> internal_position;
          for (unsigned int i = 0; i < dim; i++)
            internal_position[i] = natural_position[i];

          // The chunk model has latitude as natural coordinate. We need to convert this to colatitude
          if (dynamic_cast<const GeometryModel::Chunk<dim>*> (&this->get_geometry_model()) != nullptr && dim == 3)
            {
              internal_position[2] = numbers::PI/2. - internal_position[2];
            }

          const std::array<unsigned int,dim-1> boundary_dimensions =
            get_boundary_dimensions(boundary_indicator);

          Point<dim-1> data_position;
          for (unsigned int i = 0; i < dim-1; ++i)
            data_position[i] = internal_position[boundary_dimensions[i]];

          const Tensor<1,dim-1>  gradients = lookups.find(boundary_indicator)->second->get_gradients(data_position,component);

          if (!time_dependent)
            return gradients;

          const Tensor<1,dim-1> old_gradients = old_lookups.find(boundary_indicator)->second->get_gradients(data_position,component);

          return time_weight * gradients + (1 - time_weight) * old_gradients;
        }
      else
        return Tensor<1,dim-1>();
    }


    template <int dim>
    double
    AsciiDataBoundary<dim>::get_maximum_component_value (const types::boundary_id boundary_indicator, const unsigned int component) const
    {
      return lookups.find(boundary_indicator)->second->get_maximum_component_value(component);
    }


    template <int dim>
    void
    AsciiDataBoundary<dim>::declare_parameters (ParameterHandler  &prm,
                                                const std::string &default_directory,
                                                const std::string &default_filename,
                                                const std::string &subsection_name)
    {
      Utilities::AsciiDataBase<dim>::declare_parameters(prm,
                                                        default_directory,
                                                        default_filename,
                                                        subsection_name);

      prm.enter_subsection (subsection_name);
      {
        prm.declare_entry ("Data file name",
                           default_filename,
                           Patterns::Anything (),
                           "The file name of the model data. Provide file in format: "
                           "(File name).\\%s\\%d, where \\%s is a string specifying "
                           "the boundary of the model according to the names of the boundary "
                           "indicators (of the chosen geometry model), and \\%d is any sprintf "
                           "integer qualifier specifying the format of the current file number.");
        prm.declare_entry ("Data file time step", "1e6",
                           Patterns::Double (0.),
                           "Time step between following data files. "
                           "Depending on the setting of the global `Use years in output instead of seconds' flag "
                           "in the input file, this number is either interpreted as seconds or as years. "
                           "The default is one million, i.e., either one million seconds or one million years.");
        prm.declare_entry ("First data file model time", "0",
                           Patterns::Double (0.),
                           "Time from which on the data file with number `First data "
                           "file number' is used as boundary condition. Until this "
                           "time, a boundary condition equal to zero everywhere is assumed. "
                           "Depending on the setting of the global `Use years in output instead of seconds' flag "
                           "in the input file, this number is either interpreted as seconds or as years.");
        prm.declare_entry ("First data file number", "0",
                           Patterns::Integer (),
                           "Number of the first velocity file to be loaded when the model time "
                           "is larger than `First velocity file model time'.");
        prm.declare_entry ("Decreasing file order", "false",
                           Patterns::Bool (),
                           "In some cases the boundary files are not numbered in increasing "
                           "but in decreasing order (e.g. `Ma BP'). If this flag is set to "
                           "`True' the plugin will first load the file with the number "
                           "`First data file number' and decrease the file number during "
                           "the model run.");
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    AsciiDataBoundary<dim>::parse_parameters (ParameterHandler &prm,
                                              const std::string &subsection_name)
    {
      Utilities::AsciiDataBase<dim>::parse_parameters(prm,
                                                      subsection_name);

      prm.enter_subsection(subsection_name);
      {
        data_file_time_step             = prm.get_double ("Data file time step");
        first_data_file_model_time      = prm.get_double ("First data file model time");
        first_data_file_number          = prm.get_integer("First data file number");
        decreasing_file_order           = prm.get_bool   ("Decreasing file order");

        if (this->convert_output_to_years() == true)
          {
            data_file_time_step        *= year_in_seconds;
            first_data_file_model_time *= year_in_seconds;
          }
      }
      prm.leave_subsection();
    }



    template <int dim>
    AsciiDataLayered<dim>::AsciiDataLayered ()
    {}



    template <int dim>
    void
    AsciiDataLayered<dim>::initialize(const unsigned int components)
    {
      AssertThrow ((Plugins::plugin_type_matches<GeometryModel::SphericalShell<dim> >(this->get_geometry_model()) ||
                    Plugins::plugin_type_matches<GeometryModel::Chunk<dim> >(this->get_geometry_model()) ||
                    Plugins::plugin_type_matches<GeometryModel::Sphere<dim> >(this->get_geometry_model()) ||
                    Plugins::plugin_type_matches<GeometryModel::Box<dim> >(this->get_geometry_model())) ||
                   Plugins::plugin_type_matches<GeometryModel::TwoMergedBoxes<dim>> (this->get_geometry_model()),
                   ExcMessage ("This ascii data plugin can only be used when using "
                               "a spherical shell, chunk, sphere, box or two merged boxes geometry."));

      // Create the lookups for each file
      number_of_layer_boundaries = data_file_names.size();
      for (unsigned int i=0; i<number_of_layer_boundaries; ++i)
        {
          const std::string filename = data_directory + data_file_names[i];
          AssertThrow(Utilities::fexists(filename) || filename_is_url(filename),
                      ExcMessage (std::string("Ascii data file <")
                                  +
                                  filename
                                  +
                                  "> not found!"));

          lookups.push_back(std_cxx14::make_unique<Utilities::AsciiDataLookup<dim-1>> (components,
                                                                                       this->scale_factor));
          lookups[i]->load_file(filename,this->get_mpi_communicator());
        }
    }



    template <int dim>
    double
    AsciiDataLayered<dim>::
    get_data_component (const Point<dim> &position,
                        const unsigned int component) const
    {
      // Get the location of the component in the coordinate system of the ascii data input
      const std::array<double,dim> natural_position = this->get_geometry_model().cartesian_to_natural_coordinates(position);

      Point<dim> internal_position;
      for (unsigned int i = 0; i < dim; i++)
        internal_position[i] = natural_position[i];

      // The chunk model has latitude as natural coordinate. We need to convert this to colatitude
      if (Plugins::plugin_type_matches<GeometryModel::Chunk<dim> >(this->get_geometry_model()) && dim == 3)
        {
          internal_position[2] = numbers::PI/2. - internal_position[2];
        }

      double vertical_position;
      Point<dim-1> horizontal_position;
      if (this->get_geometry_model().natural_coordinate_system() == Utilities::Coordinates::CoordinateSystem::cartesian)
        {
          // in cartesian coordinates, the vertical component comes last
          vertical_position = internal_position[dim-1];
          for (unsigned int i = 0; i < dim-1; i++)
            horizontal_position[i] = internal_position[i];
        }
      else
        {
          // in spherical coordinates, the vertical component comes first
          vertical_position = internal_position[0];
          for (unsigned int i = 0; i < dim-1; i++)
            horizontal_position[i] = internal_position[i+1];
        }

      // Find which layer we're in
      unsigned int layer_boundary_index=0;
      // if position < position of the first boundary layer, stop
      double old_difference_in_vertical_position = vertical_position - lookups[layer_boundary_index]->get_data(horizontal_position,0);
      double difference_in_vertical_position = old_difference_in_vertical_position;
      while (difference_in_vertical_position > 0. && layer_boundary_index < number_of_layer_boundaries-1)
        {
          ++layer_boundary_index;
          old_difference_in_vertical_position = difference_in_vertical_position;
          difference_in_vertical_position = vertical_position - lookups[layer_boundary_index]->get_data(horizontal_position,0);
        }

      if (interpolation_scheme == "piecewise constant")
        {
          return lookups[layer_boundary_index]->get_data(horizontal_position,component); // takes value from layer above
        }
      else if (interpolation_scheme == "linear")
        {
          if (difference_in_vertical_position > 0 || layer_boundary_index == 0) // if the point is above the first layer or below the last
            {
              return lookups[layer_boundary_index]->get_data(horizontal_position,component);
            }
          else
            {
              const double f = difference_in_vertical_position/(difference_in_vertical_position-old_difference_in_vertical_position);
              return ((1.-f)*lookups[layer_boundary_index]->get_data(horizontal_position,component) +
                      f*lookups[layer_boundary_index-1]->get_data(horizontal_position,component));
            }
        }
      return 0;
    }


    template <int dim>
    void
    AsciiDataLayered<dim>::declare_parameters (ParameterHandler  &prm,
                                               const std::string &default_directory,
                                               const std::string &default_filename,
                                               const std::string &subsection_name)
    {
      Utilities::AsciiDataBase<dim>::declare_parameters(prm,
                                                        default_directory,
                                                        default_filename,
                                                        subsection_name);

      prm.enter_subsection (subsection_name);
      {
        prm.declare_entry ("Data directory",
                           default_directory,
                           Patterns::DirectoryName (),
                           "The name of a directory that contains the model data. This path "
                           "may either be absolute (if starting with a `/') or relative to "
                           "the current directory. The path may also include the special "
                           "text `$ASPECT_SOURCE_DIR' which will be interpreted as the path "
                           "in which the ASPECT source files were located when ASPECT was "
                           "compiled. This interpretation allows, for example, to reference "
                           "files located in the `data/' subdirectory of ASPECT. ");

        prm.declare_entry ("Data file names",
                           default_filename,
                           Patterns::List (Patterns::Anything()),
                           "The file names of the model data (comma separated). ");

        prm.declare_entry ("Interpolation scheme", "linear",
                           Patterns::Selection("piecewise constant|linear"),
                           "Method to interpolate between layer boundaries. Select from "
                           "piecewise constant or linear. Piecewise constant takes the "
                           "value from the nearest layer boundary above the data point. "
                           "The linear option interpolates linearly between layer boundaries. "
                           "Above and below the domain given by the layer boundaries, the values are"
                           "given by the top and bottom layer boundary.");

      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    AsciiDataLayered<dim>::parse_parameters (ParameterHandler &prm,
                                             const std::string &subsection_name)
    {
      Utilities::AsciiDataBase<dim>::parse_parameters(prm, subsection_name);

      prm.enter_subsection(subsection_name);
      {
        data_directory = Utilities::expand_ASPECT_SOURCE_DIR(prm.get ("Data directory"));
        data_file_names = Utilities::split_string_list(prm.get ("Data file names"), ',');
        interpolation_scheme = prm.get("Interpolation scheme");
      }
      prm.leave_subsection();
    }



    template <int dim>
    AsciiDataInitial<dim>::AsciiDataInitial ()
    {}



    template <int dim>
    void
    AsciiDataInitial<dim>::initialize (const unsigned int components)
    {
      AssertThrow ((Plugins::plugin_type_matches<const GeometryModel::SphericalShell<dim>> (this->get_geometry_model()))
                   || (Plugins::plugin_type_matches<const GeometryModel::Chunk<dim>> (this->get_geometry_model()))
                   || (Plugins::plugin_type_matches<const GeometryModel::Sphere<dim>> (this->get_geometry_model()))
                   || (Plugins::plugin_type_matches<const GeometryModel::Box<dim>> (this->get_geometry_model()))
                   || (Plugins::plugin_type_matches<const GeometryModel::TwoMergedBoxes<dim>> (this->get_geometry_model())),
                   ExcMessage ("This ascii data plugin can only be used when using "
                               "a spherical shell, chunk, or box geometry."));

      lookup = std_cxx14::make_unique<Utilities::AsciiDataLookup<dim>> (components,
                                                                        this->scale_factor);

      const std::string filename = this->data_directory + this->data_file_name;

      this->get_pcout() << std::endl << "   Loading Ascii data initial file "
                        << filename << "." << std::endl << std::endl;


      AssertThrow(Utilities::fexists(filename) || filename_is_url(filename),
                  ExcMessage (std::string("Ascii data file <")
                              +
                              filename
                              +
                              "> not found!"));
      lookup->load_file(filename, this->get_mpi_communicator());
    }



    template <int dim>
    double
    AsciiDataInitial<dim>::
    get_data_component (const Point<dim>                    &position,
                        const unsigned int                   component) const
    {
      Point<dim> internal_position = position;

      if (Plugins::plugin_type_matches<const GeometryModel::SphericalShell<dim>> (this->get_geometry_model())
          || (Plugins::plugin_type_matches<const GeometryModel::Chunk<dim>> (this->get_geometry_model())))
        {
          const std::array<double,dim> spherical_position =
            Utilities::Coordinates::cartesian_to_spherical_coordinates(position);

          for (unsigned int i = 0; i < dim; i++)
            internal_position[i] = spherical_position[i];
        }
      return lookup->get_data(internal_position,component);
    }



    template <int dim>
    AsciiDataProfile<dim>::AsciiDataProfile ()
    {}



    template <int dim>
    void
    AsciiDataProfile<dim>::initialize (const MPI_Comm &communicator)
    {
      lookup = std_cxx14::make_unique<Utilities::AsciiDataLookup<1>> (this->scale_factor);

      const std::string filename = this->data_directory + this->data_file_name;

      AssertThrow(Utilities::fexists(filename) || filename_is_url(filename),
                  ExcMessage (std::string("Ascii data file <")
                              +
                              filename
                              +
                              "> not found!"));
      lookup->load_file(filename,communicator);
    }



    template <int dim>
    std::vector<std::string>
    AsciiDataProfile<dim>::get_column_names() const
    {
      return lookup->get_column_names();
    }



    template <int dim>
    unsigned int
    AsciiDataProfile<dim>::get_column_index_from_name(const std::string &column_name) const
    {
      return lookup->get_column_index_from_name(column_name);
    }



    template <int dim>
    unsigned int
    AsciiDataProfile<dim>::maybe_get_column_index_from_name(const std::string &column_name) const
    {
      try
        {
          // read the entries in if they exist
          return lookup->get_column_index_from_name(column_name);
        }
      catch (...)
        {
          // return an invalid unsigned int entry if the column does not exist
          return numbers::invalid_unsigned_int;
        }
    }



    template <int dim>
    std::string
    AsciiDataProfile<dim>::get_column_name_from_index(const unsigned int column_index) const
    {
      return lookup->get_column_name_from_index(column_index);
    }



    template <int dim>
    double
    AsciiDataProfile<dim>::
    get_data_component (const Point<1>                      &position,
                        const unsigned int                   component) const
    {
      return lookup->get_data(position,component);
    }



    template <int dim>
    const std::vector<double> &
    AsciiDataProfile<dim>::get_coordinates() const
    {
      return lookup->get_coordinates(0);
    }



    double
    weighted_p_norm_average ( const std::vector<double> &weights,
                              const std::vector<double> &values,
                              const double p)
    {
      // TODO: prevent division by zero for all
      double averaged_parameter = 0.0;

      // first look at the special cases which can be done faster
      if (p <= -1000)
        {
          // Minimum
          double min_value = 0;
          unsigned int first_element_with_nonzero_weight = 0;
          for (; first_element_with_nonzero_weight < weights.size(); ++first_element_with_nonzero_weight)
            if (weights[first_element_with_nonzero_weight] > 0)
              {
                min_value = values[first_element_with_nonzero_weight];
                break;
              }
          Assert (first_element_with_nonzero_weight < weights.size(),
                  ExcMessage ("There are only zero (or smaller) weights in the weights vector."));

          for (unsigned int i=first_element_with_nonzero_weight+1; i < weights.size(); ++i)
            if (weights[i] != 0)
              if (values[i] < min_value)
                min_value = values[i];

          return min_value;
        }
      else if (p == -1)
        {
          // Harmonic average
          for (unsigned int i=0; i< weights.size(); ++i)
            {
              /**
               * if the value is zero, we get a division by zero. To prevent this
               * we look at what should happen in this case. When a value is zero,
               * and the correspondent weight is non-zero, this corresponds to no
               * resistance in a parallel system. This means that this will dominate,
               * and we should return zero. If the value is zero and the weight is
               * zero, we just ignore it.
               */
              if (values[i] == 0 && weights[i] > 0)
                return 0;
              else if (values[i] != 0)
                averaged_parameter += weights[i]/values[i];
            }

          Assert (averaged_parameter > 0, ExcMessage ("The sum of the weights/values may not be smaller or equal to zero."));
          const double sum_of_weights = std::accumulate(weights.begin(), weights.end(), 0.0);
          return sum_of_weights/averaged_parameter;
        }
      else if (p == 0)
        {
          // Geometric average
          for (unsigned int i=0; i < weights.size(); ++i)
            averaged_parameter += weights[i]*std::log(values[i]);

          const double sum_of_weights = std::accumulate(weights.begin(), weights.end(), 0.0);
          Assert (sum_of_weights != 0,
                  ExcMessage ("The sum of the weights may not be equal to zero, because we need to divide through it."));
          return std::exp(averaged_parameter/sum_of_weights);
        }
      else if (p == 1)
        {
          // Arithmetic average
          for (unsigned int i=0; i< weights.size(); ++i)
            averaged_parameter += weights[i]*values[i];

          const double sum_of_weights = std::accumulate(weights.begin(), weights.end(), 0.0);
          Assert (sum_of_weights != 0,
                  ExcMessage ("The sum of the weights may not be equal to zero, because we need to divide through it."));
          return averaged_parameter/sum_of_weights;
        }
      else if (p == 2)
        {
          // Quadratic average (RMS)
          for (unsigned int i=0; i< weights.size(); ++i)
            averaged_parameter += weights[i]*values[i]*values[i];

          const double sum_of_weights = std::accumulate(weights.begin(), weights.end(), 0.0);
          Assert (sum_of_weights != 0,
                  ExcMessage ("The sum of the weights may not be equal to zero, because we need to divide through it."));
          Assert (averaged_parameter/sum_of_weights > 0, ExcMessage ("The sum of the weights is smaller or equal to zero."));
          return std::sqrt(averaged_parameter/sum_of_weights);
        }
      else if (p == 3)
        {
          // Cubic average
          for (unsigned int i=0; i< weights.size(); ++i)
            averaged_parameter += weights[i]*values[i]*values[i]*values[i];

          const double sum_of_weights = std::accumulate(weights.begin(), weights.end(), 0.0);
          Assert (sum_of_weights != 0,
                  ExcMessage ("The sum of the weights may not be equal to zero, because we need to divide through it."));
          return cbrt(averaged_parameter/sum_of_weights);
        }
      else if (p >= 1000)
        {
          // Maximum
          double max_value = 0;
          unsigned int first_element_with_nonzero_weight = 0;
          for (; first_element_with_nonzero_weight < weights.size(); ++first_element_with_nonzero_weight)
            if (weights[first_element_with_nonzero_weight] > 0)
              {
                max_value = values[first_element_with_nonzero_weight];
                break;
              }
          Assert (first_element_with_nonzero_weight < weights.size(),
                  ExcMessage ("There are only zero (or smaller) weights in the weights vector."));

          for (unsigned int i=first_element_with_nonzero_weight+1; i < weights.size(); ++i)
            if (weights[i] != 0)
              if (values[i] > max_value)
                max_value = values[i];

          return max_value;
        }
      else
        {
          for (unsigned int i=0; i< weights.size(); ++i)
            {
              /**
               * When a value is zero or smaller, the exponent is smaller then one and the
               * correspondent  weight is non-zero, this corresponds to no resistance in a
               * parallel system.  This means that this 'path' will be followed, and we
               * return zero.
               */
              if (values[i] <= 0 && p < 0)
                return 0;
              averaged_parameter += weights[i] * std::pow(values[i],p);
            }

          const double sum_of_weights = std::accumulate(weights.begin(), weights.end(), 0.0);

          Assert (sum_of_weights > 0, ExcMessage ("The sum of the weights may not be smaller or equal to zero."));
          Assert (averaged_parameter > 0,
                  ExcMessage ("The sum of the weights times the values to the power p may not be smaller or equal to zero."));
          return std::pow(averaged_parameter/sum_of_weights, 1/p);
        }
    }



    template <typename T>
    T
    derivative_of_weighted_p_norm_average (const double /*averaged_parameter*/,
                                           const std::vector<double> &weights,
                                           const std::vector<double> &values,
                                           const std::vector<T> &derivatives,
                                           const double p)
    {
      // TODO: use averaged_parameter to speed up computation?
      // TODO: add special cases p = 2 and p = 3
      double averaged_parameter_derivative_part_1 = 0.0;
      T averaged_parameter_derivative_part_2 = T();

      // first look at the special cases which can be done faster
      if (p <= -1000)
        {
          // Minimum
          double min_value = 0;
          unsigned int element_with_minimum_value = 0;
          unsigned int first_element_with_nonzero_weight = 0;
          for (; first_element_with_nonzero_weight < weights.size(); ++first_element_with_nonzero_weight)
            if (weights[first_element_with_nonzero_weight] > 0)
              {
                min_value = values[first_element_with_nonzero_weight];
                element_with_minimum_value = first_element_with_nonzero_weight;
                break;
              }
          Assert (first_element_with_nonzero_weight < weights.size(),
                  ExcMessage ("There are only zero (or smaller) weights in the weights vector."));

          for (unsigned int i=first_element_with_nonzero_weight+1; i < weights.size(); ++i)
            if (weights[i] != 0)
              if (values[i] < min_value)
                {
                  min_value = values[i];
                  element_with_minimum_value = i;
                }
          return derivatives[element_with_minimum_value];
        }
      else if (p == -1)
        {
          // Harmonic average
          for (unsigned int i=0; i< weights.size(); ++i)
            {
              /**
               * if the value is zero, we get a division by zero. To prevent this
               * we look at what should happen in this case. When a value is zero,
               * and the correspondent weight is non-zero, this corresponds to no
               * resistance in a parallel system. This means that this will dominate,
               * and we should return this derivative. If the value is zero and the
               * weight is zero, we just ignore it.
               */
              if (values[i] == 0 && weights[i] > 0)
                return derivatives[i];
              else if (values[i] != 0)
                {
                  averaged_parameter_derivative_part_1 += weights[i] / values[i];
                  averaged_parameter_derivative_part_2 += weights[i] * (1/(values[i] * values[i])) * derivatives[i];
                }
            }
          const double sum_of_weights = std::accumulate(weights.begin(), weights.end(), 0.0);
          Assert (sum_of_weights > 0, ExcMessage ("The sum of the weights may not be smaller or equal to zero."));
          return std::pow(averaged_parameter_derivative_part_1/sum_of_weights,-2) * averaged_parameter_derivative_part_2/sum_of_weights;
        }
      else if (p == 0)
        {
          // Geometric average
          for (unsigned int i=0; i < weights.size(); ++i)
            {
              averaged_parameter_derivative_part_1 += weights[i]*std::log(values[i]);
              averaged_parameter_derivative_part_2 += weights[i]*(1/values[i])*derivatives[i];
            }

          const double sum_of_weights = std::accumulate(weights.begin(), weights.end(), 0.0);
          Assert (sum_of_weights != 0,
                  ExcMessage ("The sum of the weights may not be equal to zero, because we need to divide through it."));
          return std::exp(averaged_parameter_derivative_part_1/sum_of_weights) * averaged_parameter_derivative_part_2/sum_of_weights;
        }
      else if (p == 1)
        {
          // Arithmetic average
          for (unsigned int i=0; i< weights.size(); ++i)
            averaged_parameter_derivative_part_2 += weights[i]*derivatives[i];

          const double sum_of_weights = std::accumulate(weights.begin(), weights.end(), 0.0);
          Assert (sum_of_weights != 0,
                  ExcMessage ("The sum of the weights may not be equal to zero, because we need to divide through it."));
          return averaged_parameter_derivative_part_2/sum_of_weights;
        }
      else if (p >= 1000)
        {
          // Maximum
          double max_value = 0;
          unsigned int element_with_maximum_value = 0;
          unsigned int first_element_with_nonzero_weight = 0;
          for (; first_element_with_nonzero_weight < weights.size(); ++first_element_with_nonzero_weight)
            if (weights[first_element_with_nonzero_weight] > 0)
              {
                max_value = values[first_element_with_nonzero_weight];
                element_with_maximum_value = first_element_with_nonzero_weight;
                break;
              }
          Assert (first_element_with_nonzero_weight < weights.size(),
                  ExcMessage ("There are only zero (or smaller) weights in the weights vector."));

          for (unsigned int i=first_element_with_nonzero_weight+1; i < weights.size(); ++i)
            if (weights[i] != 0)
              if (values[i] > max_value)
                {
                  max_value = values[i];
                  element_with_maximum_value = i;
                }

          return derivatives[element_with_maximum_value];
        }
      else
        {
          // The general case: We can simplify the equation by stating that (1/p) * p = 1
          // TODO: This can probably be optimized by using:
          // averaged_parameter_derivative_part_2 += weights[i]*values_p[i]*(1/values[i])*derivatives[i]; and
          // averaged_parameter_derivative = averaged_parameter * (1/averaged_parameter_derivative_part_1) * averaged_parameter_derivative_part_2;
          for (unsigned int i=0; i< weights.size(); ++i)
            {
              /**
               * When a value is zero or smaller, the exponent is smaller then one and the
               * correspondent  weight is non-zero, this corresponds to no resistance in a
               * parallel system. This means that this 'path' will be followed, and we
               * return that derivative.
               */
              if (values[i] <= 0 && p < 0)
                return derivatives[i];
              averaged_parameter_derivative_part_1 += weights[i] * std::pow(values[i],p);
              averaged_parameter_derivative_part_2 += weights[i] * std::pow(values[i],p-1) * derivatives[i];
            }
          const double sum_of_weights = std::accumulate(weights.begin(), weights.end(), 0.0);
          Assert (sum_of_weights > 0, ExcMessage ("The sum of the weights may not be smaller or equal to zero."));
          Assert (averaged_parameter_derivative_part_1/sum_of_weights > 0,
                  ExcMessage ("The sum of the weights times the values to the power p may not be smaller or equal to zero."));
          return std::pow(averaged_parameter_derivative_part_1/sum_of_weights,(1/p)-1) * averaged_parameter_derivative_part_2/sum_of_weights;
          // TODO: find a way to check if value is finite for any type? Or just leave this kind of checking up to the user?
        }
    }



    template <int dim>
    double compute_spd_factor(const double eta,
                              const SymmetricTensor<2,dim> &strain_rate,
                              const SymmetricTensor<2,dim> &dviscosities_dstrain_rate,
                              const double SPD_safety_factor)
    {
      // if the strain rate is zero, or the derivative is zero, then
      // the exact choice of alpha factor does not matter because the
      // factor that it multiplies is zero -- so return the best value
      // (i.e., one)
      if ((strain_rate.norm() == 0) || (dviscosities_dstrain_rate.norm() == 0))
        return 1;

      const double norm_a_b = std::sqrt((strain_rate*strain_rate)*(dviscosities_dstrain_rate*dviscosities_dstrain_rate));;//std::sqrt((deviator(strain_rate)*deviator(strain_rate))*(dviscosities_dstrain_rate*dviscosities_dstrain_rate));
      const double contract_b_a = (dviscosities_dstrain_rate*strain_rate);
      const double one_minus_part = 1 - (contract_b_a / norm_a_b);
      const double denom = one_minus_part * one_minus_part * norm_a_b;

      // the case denom == 0 (smallest eigenvalue is zero), should return one,
      // and it does here, because C_safety * 2.0 * eta is always larger then zero.
      if (denom <= SPD_safety_factor * 2.0 * eta)
        return 1.0;
      else
        return std::max(0.0, SPD_safety_factor * ((2.0 * eta) / denom));
    }



    template <int dim>
    Point<dim> convert_array_to_point(const std::array<double,dim> &array)
    {
      Point<dim> point;
      for (unsigned int i = 0; i < dim; i++)
        point[i] = array[i];

      return point;
    }



    template <int dim>
    std::array<double,dim> convert_point_to_array(const Point<dim> &point)
    {
      std::array<double,dim> array;
      for (unsigned int i = 0; i < dim; i++)
        array[i] = point[i];

      return array;
    }



    Operator::Operator()
      :
      op(uninitialized)
    {}



    Operator::Operator(const operation _op)
      :
      op(_op)
    {}



    double
    Operator::operator() (const double x, const double y) const
    {
      switch (op)
        {
          case Utilities::Operator::add:
          {
            return x + y;
          }
          case Utilities::Operator::subtract:
          {
            return x - y;
          }
          case Utilities::Operator::minimum:
          {
            return std::min(x,y);
          }
          case Utilities::Operator::maximum:
          {
            return std::max(x,y);
          }
          case Utilities::Operator::replace_if_valid:
          {
            if (std::isnan(y))
              return x;
            else
              return y;
          }
          default:
          {
            Assert (false, ExcInternalError());
          }
        }
      return numbers::signaling_nan<double>();
    }



    bool
    Operator::operator== (const operation other_op) const
    {
      return other_op == op;
    }



    std::vector<Operator> create_model_operator_list(const std::vector<std::string> &operator_names)
    {
      std::vector<Operator> operator_list(operator_names.size());
      for (unsigned int i=0; i<operator_names.size(); ++i)
        {
          // create operator list
          if (operator_names[i] == "add")
            operator_list[i] = Operator(Operator::add);
          else if (operator_names[i] == "subtract")
            operator_list[i] = Operator(Operator::subtract);
          else if (operator_names[i] == "minimum")
            operator_list[i] = Operator(Operator::minimum);
          else if (operator_names[i] == "maximum")
            operator_list[i] = Operator(Operator::maximum);
          else if (operator_names[i] == "replace if valid")
            operator_list[i] = Operator(Operator::replace_if_valid);
          else
            AssertThrow(false,
                        ExcMessage ("ASPECT only accepts the following operators: "
                                    "add, subtract, minimum, maximum, and replace if valid. But your parameter file "
                                    "contains: " + operator_names[i] + ". Please check your parameter file.") );
        }

      return operator_list;
    }



    const std::string get_model_operator_options()
    {
      return "add|subtract|minimum|maximum|replace if valid";
    }



    template <int dim>
    SymmetricTensor<2,dim>
    nth_basis_for_symmetric_tensors (const unsigned int k)
    {
      Assert((k < SymmetricTensor<2,dim>::n_independent_components),
             ExcMessage("The component is larger then the amount of independent components in the matrix.") );

      const TableIndices<2> indices_ij = SymmetricTensor<2,dim>::unrolled_to_component_indices (k);

      Tensor<2,dim> result;
      result[indices_ij] = 1;

      return symmetrize(result);
    }

    template <int dim>
    NaturalCoordinate<dim>::NaturalCoordinate(Point<dim> &position,
                                              const GeometryModel::Interface<dim> &geometry_model)
    {
      coordinate_system = geometry_model.natural_coordinate_system();
      coordinates = geometry_model.cartesian_to_natural_coordinates(position);
    }

    template <int dim>
    NaturalCoordinate<dim>::NaturalCoordinate(const std::array<double, dim> &coord,
                                              const Utilities::Coordinates::CoordinateSystem &coord_system) :
      coordinate_system (coord_system), coordinates (coord)
    {}



    template <int dim>
    std::array<double,dim> &NaturalCoordinate<dim>::get_coordinates()
    {
      return coordinates;
    }

    template <>
    std::array<double,1> NaturalCoordinate<2>::get_surface_coordinates() const
    {
      std::array<double,1> coordinate;

      switch (coordinate_system)
        {
          case Coordinates::CoordinateSystem::cartesian:
            coordinate[0] = coordinates[0];
            break;

          case Coordinates::CoordinateSystem::spherical:
            coordinate[0] = coordinates[1];
            break;

          case Coordinates::CoordinateSystem::ellipsoidal:
            coordinate[0] = coordinates[1];
            break;

          default:
            coordinate[0] = 0;
            Assert (false, ExcNotImplemented());
            break;
        }

      return coordinate;
    }

    template <>
    std::array<double,2> NaturalCoordinate<3>::get_surface_coordinates() const
    {
      std::array<double,2> coordinate;

      switch (coordinate_system)
        {
          case Coordinates::CoordinateSystem::cartesian:
            coordinate[0] = coordinates[0];
            coordinate[1] = coordinates[1];
            break;

          case Coordinates::CoordinateSystem::spherical:
            coordinate[0] = coordinates[1];
            coordinate[1] = coordinates[2];
            break;

          case Coordinates::CoordinateSystem::ellipsoidal:
            coordinate[0] = coordinates[1];
            coordinate[1] = coordinates[2];
            break;

          default:
            Assert (false, ExcNotImplemented());
        }

      return coordinate;
    }

    template <int dim>
    double NaturalCoordinate<dim>::get_depth_coordinate() const
    {
      switch (coordinate_system)
        {
          case Coordinates::CoordinateSystem::cartesian:
            return coordinates[dim-1];

          case Coordinates::CoordinateSystem::spherical:
            return coordinates[0];

          case Coordinates::CoordinateSystem::ellipsoidal:
            return coordinates[0];

          default:
            Assert (false, ExcNotImplemented());
        }

      return 0;
    }


    template <int dim, typename VectorType>
    void
    project_cellwise(const Mapping<dim>                                        &mapping,
                     const DoFHandler<dim>                                     &dof_handler,
                     const unsigned int                                         component_index,
                     const Quadrature<dim>                                     &quadrature,
                     const std::function<void(
                       const typename DoFHandler<dim>::active_cell_iterator &,
                       const std::vector<Point<dim> > &,
                       std::vector<double> &)>                                 &function,
                     VectorType                                                &vec_result)
    {
      const FEValuesExtractors::Scalar extractor(component_index);

      UpdateFlags update_flags = UpdateFlags(update_values   |
                                             update_quadrature_points |
                                             update_JxW_values);

      FEValues<dim> fe_values (mapping, dof_handler.get_fe(), quadrature, update_flags);

      const unsigned int
      dofs_per_cell = fe_values.dofs_per_cell,
      n_q_points    = fe_values.n_quadrature_points;

      std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
      Vector<double> cell_vector (dofs_per_cell);
      Vector<double> local_projection (dofs_per_cell);
      FullMatrix<double> local_mass_matrix (dofs_per_cell, dofs_per_cell);

      std::vector<double> rhs_values(n_q_points);

      for (const auto &cell : dof_handler.active_cell_iterators())
        if (cell->is_locally_owned())
          {
            // For each cell, create a local mass matrix and rhs.
            cell->get_dof_indices (local_dof_indices);
            fe_values.reinit(cell);

            function(cell, fe_values.get_quadrature_points(), rhs_values);

            cell_vector = 0;
            local_mass_matrix = 0;
            for (unsigned int point=0; point<n_q_points; ++point)
              for (unsigned int i=0; i<dofs_per_cell; ++i)
                {
                  if (dof_handler.get_fe().system_to_component_index(i).first == component_index)
                    cell_vector(i) +=
                      rhs_values[point] *
                      fe_values[extractor].value(i,point) *
                      fe_values.JxW(point);

                  for (unsigned int j=0; j<dofs_per_cell; ++j)
                    if ((dof_handler.get_fe().system_to_component_index(i).first ==
                         component_index)
                        &&
                        (dof_handler.get_fe().system_to_component_index(j).first ==
                         component_index))
                      local_mass_matrix(j,i) += (fe_values[extractor].value(i,point) * fe_values[extractor].value(j,point) *
                                                 fe_values.JxW(point));
                    else if (i == j)
                      local_mass_matrix(i,j) = 1.;
                }

            // now invert the local mass matrix and multiply it with the rhs
            local_mass_matrix.gauss_jordan();
            local_mass_matrix.vmult (local_projection, cell_vector);

            // then set the global solution vector to the values just computed
            cell->set_dof_values (local_projection, vec_result);
          }
    }




    template <int dim>
    VectorFunctionFromVelocityFunctionObject<dim>::
    VectorFunctionFromVelocityFunctionObject
    (const unsigned int n_components,
     const std::function<Tensor<1,dim> (const Point<dim> &)> &function_object)
      :
      Function<dim>(n_components),
      function_object (function_object)
    {
    }



    template <int dim>
    double
    VectorFunctionFromVelocityFunctionObject<dim>::value (const Point<dim> &p,
                                                          const unsigned int component) const
    {
      Assert (component < this->n_components,
              ExcIndexRange (component, 0, this->n_components));

      if (component < dim)
        {
          const Tensor<1,dim> v = function_object(p);
          return v[component];
        }
      else
        return 0;
    }



    template <int dim>
    void
    VectorFunctionFromVelocityFunctionObject<dim>::
    vector_value (const Point<dim>   &p,
                  Vector<double>     &values) const
    {
      AssertDimension(values.size(), this->n_components);

      // set everything to zero, and then the right components to their correct values
      values = 0;

      const Tensor<1,dim> v = function_object(p);
      for (unsigned int d=0; d<dim; ++d)
        values(d) = v[d];
    }

// Explicit instantiations

#define INSTANTIATE(dim) \
  template \
  IndexSet extract_locally_active_dofs_with_component(const DoFHandler<dim> &, \
                                                      const ComponentMask &); \
  template \
  std::vector<std::string> \
  expand_dimensional_variable_names<dim> (const std::vector<std::string> &var_declarations); \
  \
  template \
  class VectorFunctionFromVelocityFunctionObject<dim>;

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE



    template class AsciiDataLookup<1>;
    template class AsciiDataLookup<2>;
    template class AsciiDataLookup<3>;
    template class AsciiDataBase<2>;
    template class AsciiDataBase<3>;
    template class AsciiDataBoundary<2>;
    template class AsciiDataBoundary<3>;
    template class AsciiDataLayered<2>;
    template class AsciiDataLayered<3>;
    template class AsciiDataInitial<2>;
    template class AsciiDataInitial<3>;
    template class AsciiDataProfile<1>;
    template class AsciiDataProfile<2>;
    template class AsciiDataProfile<3>;

    template Point<2> Coordinates::spherical_to_cartesian_coordinates<2>(const std::array<double,2> &scoord);
    template Point<3> Coordinates::spherical_to_cartesian_coordinates<3>(const std::array<double,3> &scoord);

    template std::array<double,2> Coordinates::cartesian_to_spherical_coordinates<2>(const Point<2> &position);
    template std::array<double,3> Coordinates::cartesian_to_spherical_coordinates<3>(const Point<3> &position);

    template Tensor<1,2> Coordinates::spherical_to_cartesian_vector<2>(const Tensor<1,2> &spherical_vector,
                                                                       const Point<2> &position);
    template Tensor<1,3> Coordinates::spherical_to_cartesian_vector<3>(const Tensor<1,3> &spherical_vector,
                                                                       const Point<3> &position);


    template std::array<double,2> Coordinates::WGS84_coordinates<2>(const Point<2> &position);
    template std::array<double,3> Coordinates::WGS84_coordinates<3>(const Point<3> &position);

    template bool polygon_contains_point<2>(const std::vector<Point<2> > &pointList, const dealii::Point<2> &point);
    template bool polygon_contains_point<3>(const std::vector<Point<2> > &pointList, const dealii::Point<2> &point);

    template double signed_distance_to_polygon<2>(const std::vector<Point<2> > &pointList, const dealii::Point<2> &point);
    template double signed_distance_to_polygon<3>(const std::vector<Point<2> > &pointList, const dealii::Point<2> &point);


    template std::array<Tensor<1,2>,1> orthogonal_vectors (const Tensor<1,2> &v);
    template std::array<Tensor<1,3>,2> orthogonal_vectors (const Tensor<1,3> &v);

    template double
    derivative_of_weighted_p_norm_average (const double averaged_parameter,
                                           const std::vector<double> &weights,
                                           const std::vector<double> &values,
                                           const std::vector<double> &derivatives,
                                           const double p);

    template dealii::SymmetricTensor<2, 2, double>
    derivative_of_weighted_p_norm_average (const double averaged_parameter,
                                           const std::vector<double> &weights,
                                           const std::vector<double> &values,
                                           const std::vector<dealii::SymmetricTensor<2, 2, double> > &derivatives,
                                           const double p);

    template dealii::SymmetricTensor<2, 3, double>
    derivative_of_weighted_p_norm_average (const double averaged_parameter,
                                           const std::vector<double> &weights,
                                           const std::vector<double> &values,
                                           const std::vector<dealii::SymmetricTensor<2, 3, double> > &derivatives,
                                           const double p);

    template double compute_spd_factor(const double eta,
                                       const SymmetricTensor<2,2> &strain_rate,
                                       const SymmetricTensor<2,2> &dviscosities_dstrain_rate,
                                       const double safety_factor);

    template double compute_spd_factor(const double eta,
                                       const SymmetricTensor<2,3> &strain_rate,
                                       const SymmetricTensor<2,3> &dviscosities_dstrain_rate,
                                       const double safety_factor);

    template Point<2> convert_array_to_point<2>(const std::array<double,2> &array);
    template Point<3> convert_array_to_point<3>(const std::array<double,3> &array);

    template std::array<double,2> convert_point_to_array<2>(const Point<2> &point);
    template std::array<double,3> convert_point_to_array<3>(const Point<3> &point);

    template SymmetricTensor<2,2> nth_basis_for_symmetric_tensors (const unsigned int k);
    template SymmetricTensor<2,3> nth_basis_for_symmetric_tensors (const unsigned int k);

    template class NaturalCoordinate<2>;
    template class NaturalCoordinate<3>;


    template Table<2,double> parse_input_table(const std::string &input_string,
                                               const unsigned int n_rows,
                                               const unsigned int n_columns,
                                               const std::string &property_name);

    template
    void
    project_cellwise(const Mapping<2> &mapping,
                     const DoFHandler<2> &dof_handler,
                     const unsigned int component_index,
                     const Quadrature<2> &quadrature,
                     const std::function<void(
                       const DoFHandler<2>::active_cell_iterator &,
                       const std::vector<Point<2> > &,
                       std::vector<double> &)> &function,
                     dealii::LinearAlgebra::distributed::Vector<double> &vec_result);
    template
    void
    project_cellwise(const Mapping<3> &mapping,
                     const DoFHandler<3> &dof_handler,
                     const unsigned int component_index,
                     const Quadrature<3> &quadrature,
                     const std::function<void(
                       const DoFHandler<3>::active_cell_iterator &,
                       const std::vector<Point<3> > &,
                       std::vector<double> &)> &function,
                     dealii::LinearAlgebra::distributed::Vector<double> &vec_result);
    template
    void
    project_cellwise(const Mapping<2> &mapping,
                     const DoFHandler<2> &dof_handler,
                     const unsigned int component_index,
                     const Quadrature<2> &quadrature,
                     const std::function<void(
                       const DoFHandler<2>::active_cell_iterator &,
                       const std::vector<Point<2> > &,
                       std::vector<double> &)> &function,
                     LinearAlgebra::BlockVector &vec_result);
    template
    void
    project_cellwise(const Mapping<3> &mapping,
                     const DoFHandler<3> &dof_handler,
                     const unsigned int component_index,
                     const Quadrature<3> &quadrature,
                     const std::function<void(
                       const DoFHandler<3>::active_cell_iterator &,
                       const std::vector<Point<3> > &,
                       std::vector<double> &)> &function,
                     LinearAlgebra::BlockVector &vec_result);
  }
}
