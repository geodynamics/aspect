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
#include <aspect/utilities.h>
#include <aspect/simulator_access.h>
#include <aspect/geometry_model/interface.h>

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
#include <deal.II/grid/grid_tools.h>

#include <cerrno>
#include <dirent.h>
#include <fstream>
#include <string>
#include <sys/stat.h>
#include <iostream>
#include <regex>

#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/copy.hpp>

namespace aspect
{
  /**
   * A namespace for utility functions that might be used in many different
   * places to prevent code duplication.
   */
  namespace Utilities
  {
    namespace internal
    {
      namespace MPI
      {
        // --------------------------------------------------------------------
        // The following is copied from deal.II's mpi.templates.h file.
        // We should instead import it from deal.II's header files directly
        // if that information is made available via one of the existing .h
        // files.
        // --------------------------------------------------------------------
#ifdef DEAL_II_WITH_MPI
        /**
         * Return the corresponding MPI data type id for the argument given.
         */
        inline MPI_Datatype
        mpi_type_id(const bool *)
        {
          return MPI_CXX_BOOL;
        }



        inline MPI_Datatype
        mpi_type_id(const char *)
        {
          return MPI_CHAR;
        }



        inline MPI_Datatype
        mpi_type_id(const signed char *)
        {
          return MPI_SIGNED_CHAR;
        }



        inline MPI_Datatype
        mpi_type_id(const short *)
        {
          return MPI_SHORT;
        }



        inline MPI_Datatype
        mpi_type_id(const int *)
        {
          return MPI_INT;
        }



        inline MPI_Datatype
        mpi_type_id(const long int *)
        {
          return MPI_LONG;
        }



        inline MPI_Datatype
        mpi_type_id(const unsigned char *)
        {
          return MPI_UNSIGNED_CHAR;
        }



        inline MPI_Datatype
        mpi_type_id(const unsigned short *)
        {
          return MPI_UNSIGNED_SHORT;
        }



        inline MPI_Datatype
        mpi_type_id(const unsigned int *)
        {
          return MPI_UNSIGNED;
        }



        inline MPI_Datatype
        mpi_type_id(const unsigned long int *)
        {
          return MPI_UNSIGNED_LONG;
        }



        inline MPI_Datatype
        mpi_type_id(const unsigned long long int *)
        {
          return MPI_UNSIGNED_LONG_LONG;
        }



        inline MPI_Datatype
        mpi_type_id(const float *)
        {
          return MPI_FLOAT;
        }



        inline MPI_Datatype
        mpi_type_id(const double *)
        {
          return MPI_DOUBLE;
        }



        inline MPI_Datatype
        mpi_type_id(const long double *)
        {
          return MPI_LONG_DOUBLE;
        }



        inline MPI_Datatype
        mpi_type_id(const std::complex<float> *)
        {
          return MPI_COMPLEX;
        }



        inline MPI_Datatype
        mpi_type_id(const std::complex<double> *)
        {
          return MPI_DOUBLE_COMPLEX;
        }
#endif
      }
    }






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



    namespace MapParsing
    {
      namespace
      {
        // This is a helper function used in parse_map_to_double_array below.
        // It takes an input_string that is expected to follow the input format
        // explained in the documentation of the parse_map_to_double_array function
        // and parses it into a multimap, only performing rudimentary error checking
        // for correct formatting.
        std::multimap<std::string, double>
        parse_string_to_map (const std::string &input_string,
                             const Options &options)
        {
          std::multimap<std::string, double> parsed_map;

          // Parse the input string, if it follows the structure of
          // 'key1:value1, key2:value2', or 'key1:value1|value2, ...'.
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
                                           + options.property_name
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
                                               + options.property_name
                                               + " is only allowed if there is no other "
                                               "keyword."));

                      const std::vector<std::string> values = dealii::Utilities::split_string_list(key_and_value[1], '|');

                      // Assign all the values to all fields
                      for (const std::string &key: options.list_of_required_keys)
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
                                               + options.property_name
                                               + " is listed multiple times. "
                                               "Check that you have only one value for "
                                               "each field id in your list."));

                      const std::vector<std::string> values = dealii::Utilities::split_string_list(key_and_value[1], '|');

                      for (const auto &value : values)
                        parsed_map.emplace(key_and_value[0],Utilities::string_to_double(value));
                    }
                }
            }
          // Parse the input string, if it follows the structure of
          // 'value1, value2, value3' with as many entries as allowed keys.
          else if (Patterns::List(Patterns::Double(),options.list_of_allowed_keys.size(),options.list_of_allowed_keys.size()).match(input_string))
            {
              const std::vector<double> values = possibly_extend_from_1_to_N(dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(input_string)),
                                                                             options.list_of_allowed_keys.size(),
                                                                             options.property_name);

              for (unsigned int i=0; i<values.size(); ++i)
                {
                  // list_of_keys and values have the same length, which is guaranteed by the
                  // call to possibly_extend_from_1_to_N() above
                  parsed_map.emplace(options.list_of_allowed_keys[i],values[i]);
                }
            }
          // Parse the input string, if it follows the structure of
          // 'value1, value2, value3' with one entry or as many entries as required keys.
          else if (Patterns::List(Patterns::Double(),1,options.list_of_required_keys.size()).match(input_string))
            {
              const std::vector<double> values = possibly_extend_from_1_to_N (dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(input_string)),
                                                                              options.list_of_required_keys.size(),
                                                                              options.property_name);

              for (unsigned int i=0; i<values.size(); ++i)
                {
                  // list_of_keys and values have the same length, which is guaranteed by the
                  // call to possibly_extend_from_1_to_N() above
                  parsed_map.emplace(options.list_of_required_keys[i],values[i]);
                }
            }
          else
            {
              // No Patterns matches were found!
              AssertThrow (false,
                           ExcMessage ("The string for property <"
                                       + options.property_name
                                       + "> does not have the expected format. "
                                       + "Check that the string is either a "
                                       + "comma separated list of `<double>' or "
                                       + "`<key1> : <double>|<double>|..., "
                                       + "<key2> : <double>|... , ... '. "
                                       + "If the string looks correct, "
                                       + "it is likely that the length of the "
                                       + "list of keys passed to "
                                       + "parse_map_to_double_array does not "
                                       + "match the length of the "
                                       + "comma separated property list."));
            }

          // Now remove all allowed but not requested keys from the map
          // If no keys are specifically requested, do not delete any.
          if (options.list_of_allowed_keys.size() > options.list_of_required_keys.size())
            {
              for (const auto &allowed_key: options.list_of_allowed_keys)
                if (std::find(options.list_of_required_keys.begin(),
                              options.list_of_required_keys.end(),
                              allowed_key) == options.list_of_required_keys.end())
                  parsed_map.erase(allowed_key);
            }

          return parsed_map;
        }

        // This is a helper function used in parse_map_to_double_array below.
        // It takes an input multimap for example generated by parse_string_to_map
        // above and flattens it into a
        // vector of plain doubles (in the order of the given keys).
        std::vector<double>
        flatten_map_to_vector (const std::multimap<std::string, double> &map,
                               const std::vector<std::string> &keys)
        {
          std::vector<double> values;
          values.reserve(map.size());

          for (const std::string &field_name: keys)
            {
              const std::pair<std::multimap<std::string, double>::const_iterator,
                    std::multimap<std::string, double>::const_iterator> entry_range = map.equal_range(field_name);

              for (auto entry = entry_range.first; entry != entry_range.second; ++entry)
                values.push_back(entry->second);
            }
          return values;
        }

        void
        read_or_check_map_structure (std::multimap<std::string, double> &map,
                                     Options &options)
        {
          const unsigned int n_fields = options.list_of_required_keys.size();

          std::vector<unsigned int> values_per_key(n_fields, 0);

          if (options.check_values_per_key)
            AssertThrow(options.n_values_per_key.size() == n_fields,
                        ExcMessage("When providing an expected structure for input parameter " + options.property_name + " you need to provide "
                                   + "as many entries in the structure vector as there are input field names (+1 if there is a background field). "
                                   + "The current structure vector has " + std::to_string(options.n_values_per_key.size()) + " entries, but there are "
                                   + std::to_string(n_fields) + " field names." ));

          for (const std::pair<const std::string, double> &key_and_value: map)
            {
              const std::vector<std::string>::const_iterator field_name =
                std::find(options.list_of_required_keys.begin(),
                          options.list_of_required_keys.end(),
                          key_and_value.first);

              // Ensure that each key is in the list of field names
              AssertThrow (field_name != options.list_of_required_keys.end(),
                           ExcMessage ("The keyword <" + key_and_value.first + "> in "
                                       + options.property_name + " does not match any entries "
                                       "from the list of requested field names."
                                       "Check that you only use valid names.\n\n"
                                       "One example of where to check this is if "
                                       "Compositional fields are used, "
                                       "then check the id list "
                                       "from `set Names of fields' in the "
                                       "Compositional fields subsection. "
                                       "Alternatively, if `set Names of fields' "
                                       "is not set, the default names are "
                                       "C_1, C_2, ..., C_n."));

              const unsigned int field_index = std::distance(options.list_of_required_keys.cbegin(), field_name);
              values_per_key[field_index] += 1;
            }

          if (options.store_values_per_key)
            options.n_values_per_key = values_per_key;

          unsigned int field_index = 0;
          for (const unsigned int &n_values: values_per_key)
            {
              if (options.allow_multiple_values_per_key == false)
                AssertThrow (n_values <= 1,
                             ExcMessage ("The keyword <"
                                         + options.list_of_required_keys[field_index]
                                         + "> in "
                                         + options.property_name
                                         + " has multiple values, which is unexpected. "
                                         "Check that you have only one value for "
                                         "each field id in your list."));

              if (options.allow_missing_keys == false)
                AssertThrow (n_values > 0,
                             ExcMessage ("The keyword <"
                                         + options.list_of_required_keys[field_index]
                                         + "> in "
                                         + options.property_name
                                         + " is not listed, although it is expected. "
                                         "Check that you have at least one value for "
                                         "each field id in your list (possibly plus "
                                         "`background` if a background field is expected "
                                         "for this property)."));

              if (options.check_values_per_key)
                {
                  const unsigned int n_expected_values = options.n_values_per_key[field_index];
                  const std::string field_name = options.list_of_required_keys[field_index];

                  AssertThrow((n_expected_values == n_values || n_values == 1),
                              ExcMessage("The key <" + field_name + "> in <"+ options.property_name + "> does not have "
                                         + "the expected number of values. It expects " + std::to_string(n_expected_values)
                                         + " or 1 values, but we found " + std::to_string(n_values) + " values."));

                  // If we expect multiple values for a key, but found exactly one: assume
                  // the one value stands for every expected value. This allows
                  // for short and simpler input if all values for a key are the same.
                  if (n_values == 1)
                    {
                      const double field_value = map.find(field_name)->second;
                      for (unsigned int i=1; i<n_expected_values; ++i)
                        map.emplace(field_name, field_value);
                    }
                }

              ++field_index;
            }
        }
      }

      std::vector<double>
      parse_map_to_double_array(const std::string &input_string,
                                Options &options)
      {
        // Check options for consistency
        AssertThrow (options.property_name != "",
                     ExcMessage("parse_map_to_double_array needs a property name to be able to properly report parsing errors."));
        AssertThrow (options.list_of_required_keys.size() != 0,
                     ExcMessage("parse_map_to_double_array needs at least one required key name for property "
                                + options.property_name
                                + "."));
        AssertThrow (options.check_values_per_key == false ||
                     options.store_values_per_key == false,
                     ExcMessage("parse_map_to_double_array can not simultaneously store the structure "
                                "of the parsed map for "
                                + options.property_name
                                + " and check that structure against a given structure."));
        AssertThrow (options.check_values_per_key == false ||
                     options.n_values_per_key.size() == options.list_of_required_keys.size(),
                     ExcMessage("parse_map_to_double_array can only check the structure "
                                "of the parsed map for "
                                + options.property_name
                                + " if an expected number of values for each key is given."));

        // First: parse the string into a map depending on what Pattern we are dealing with
        std::multimap<std::string, double> parsed_map = parse_string_to_map(input_string,
                                                                            options);

        // Second: Now check that the structure of the map is as expected
        read_or_check_map_structure(parsed_map,
                                    options);

        // Finally: Convert the map into a vector of doubles, sorted in the order
        // of the list_of_required_keys option
        return flatten_map_to_vector(parsed_map, options.list_of_required_keys);
      }
    }



    std::vector<double>
    parse_map_to_double_array (const std::string &input_string,
                               const std::vector<std::string> &list_of_keys,
                               const bool expects_background_field,
                               const std::string &property_name,
                               const bool allow_multiple_values_per_key,
                               const std::unique_ptr<std::vector<unsigned int>> &n_values_per_key,
                               const bool allow_missing_keys)
    {
      std::vector<std::string> input_field_names = list_of_keys;

      if (expects_background_field)
        input_field_names.insert(input_field_names.begin(),"background");

      MapParsing::Options options(input_field_names, property_name);
      options.allow_multiple_values_per_key = allow_multiple_values_per_key;
      options.allow_missing_keys = allow_missing_keys;

      if (n_values_per_key)
        {
          options.n_values_per_key = *n_values_per_key;
          options.check_values_per_key = (n_values_per_key->size() != 0);
          options.store_values_per_key = (n_values_per_key->size() == 0);
        }

      const auto parsed_map = MapParsing::parse_map_to_double_array(input_string, options);

      if (n_values_per_key)
        *n_values_per_key = options.n_values_per_key;

      return parsed_map;
    }



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



    template <int dim>
    std::vector<Point<dim>> get_unit_support_points(const SimulatorAccess<dim> &simulator_access)
    {
      if ( !simulator_access.get_parameters().use_locally_conservative_discretization )
        {
          return simulator_access.get_fe().get_unit_support_points();
        }
      else
        {
          //special case for discontinuous pressure elements, which lack unit support points
          const unsigned int dofs_per_cell = simulator_access.get_fe().dofs_per_cell;
          std::vector<Point<dim>> unit_support_points;
          unit_support_points.reserve(dofs_per_cell);

          for (unsigned int dof=0; dof < dofs_per_cell; ++dof)
            {
              // base will hold element, base_index holds node/shape function within that element
              const unsigned int base       = simulator_access.get_fe().system_to_base_index(dof).first.first;
              const unsigned int base_index = simulator_access.get_fe().system_to_base_index(dof).second;
              // get the unit support points for the relevant element
              const std::vector<Point<dim>> &my_support_points = simulator_access.get_fe().base_element(base).get_unit_support_points();
              if ( my_support_points.size() == 0 )
                {
                  //manufacture a support point, arbitrarily at cell center
                  if (dim==2)
                    unit_support_points.push_back(Point<dim> (0.5,0.5));
                  if (dim==3)
                    unit_support_points.push_back(Point<dim> (0.5,0.5,0.5));
                }
              else
                {
                  unit_support_points.push_back(my_support_points[base_index]);
                }
            }
          return unit_support_points;
        }
    }



    template <int dim>
    bool
    point_is_in_triangulation(const Mapping<dim> &mapping,
                              const parallel::distributed::Triangulation<dim> &triangulation,
                              const Point<dim> &point,
                              const MPI_Comm mpi_communicator)
    {
      // Try to find the cell around the given point.
      bool cell_found = false;
      std::pair<const typename parallel::distributed::Triangulation<dim>::active_cell_iterator,
          Point<dim>> it =
            GridTools::find_active_cell_around_point<>(mapping, triangulation, point);

      // If we found the correct cell on this MPI process, we have found the right cell.
      if (it.first.state() == IteratorState::valid && it.first->is_locally_owned())
        cell_found = true;

      // Compute how many processes found the cell.
      const int n_procs_cell_found = Utilities::MPI::sum(cell_found ? 1 : 0, mpi_communicator);
      // If at least one process found the cell, the point is in the triangulation.
      if (n_procs_cell_found > 0)
        return true;
      else
        return false;
    }



    namespace Coordinates
    {

      template <int dim>
      std::array<double,dim>
      WGS84_coordinates(const Point<dim> &position)
      {
        Assert (dim==3, ExcNotImplemented());

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
                    * constants::radians_to_degree;

        ecoord[1] = std::atan2(position(1), position(0))
                    * constants::radians_to_degree;

        // Set all longitudes between [0,360]:
        if (ecoord[1] < 0.)
          ecoord[1] += 360.;
        else if (ecoord[1] > 360.)
          ecoord[1] -= 360.;


        ecoord[0] = radius/std::sqrt(1- ellipticity * ellipticity
                                     * std::sin(constants::degree_to_radians * ecoord[2])
                                     * std::sin(constants::degree_to_radians * ecoord[2]));
        return ecoord;
      }



      template <int dim>
      std::array<double,dim>
      cartesian_to_spherical_coordinates(const Point<dim> &position)
      {
        std::array<double,dim> scoord;

        scoord[0] = position.norm(); // R

        // Compute the longitude phi. Note that atan2 is documented to return
        // its result as a value between -pi and +pi, whereas we use the
        // convention that we consider eastern longitude between 0 and 2pi.
        // As a consequence, we correct where necessary.
        scoord[1] = std::atan2(position(1),position(0));
        if (scoord[1] < 0.0)
          scoord[1] += 2.0*numbers::PI; // correct phi to [0,2*pi]

        // In 3d also compute the polar angle (=colatitude)
        if (dim==3)
          {
            if (/* R= */scoord[0] > std::numeric_limits<double>::min())
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
        const double theta  = std::atan2(x(2) + ep * ep * b * Utilities::fixed_power<3>(std::sin(th)),
                                         (p - (eccentricity * eccentricity * R  * Utilities::fixed_power<3>(std::cos(th)))));
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
    polygon_contains_point(const std::vector<Point<2>> &point_list,
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
      for (int i=0; i<pointNo; ++i)
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
    signed_distance_to_polygon(const std::vector<Point<2>> &point_list,
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
      std::vector<Point<2>> shifted_point_list(n_poly_points);
      shifted_point_list[0] = point_list[n_poly_points-1];

      for (unsigned int i = 0; i < n_poly_points-1; ++i)
        shifted_point_list[i+1] = point_list[i];

      for (unsigned int i = 0; i < n_poly_points; ++i)
        {
          const std::array<Point<2>,2> list = {{point_list[i], shifted_point_list[i]}};
          distances[i] = distance_to_line(list, point);
        }

      // Return the minimum of the distances of the point to all polygon segments
      return *std::min_element(distances.begin(),distances.end()) * sign;
    }



    double
    distance_to_line(const std::array<dealii::Point<2>,2> &point_list,
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



    Tensor<2,3>
    rotation_matrix_from_axis (const Tensor<1,3> &rotation_axis,
                               const double rotation_angle)
    {
      Tensor<2,3> rotation_matrix;
      rotation_matrix[0][0] = (1-std::cos(rotation_angle)) * rotation_axis[0]*rotation_axis[0] + std::cos(rotation_angle);
      rotation_matrix[0][1] = (1-std::cos(rotation_angle)) * rotation_axis[0]*rotation_axis[1] - rotation_axis[2] * std::sin(rotation_angle);
      rotation_matrix[0][2] = (1-std::cos(rotation_angle)) * rotation_axis[0]*rotation_axis[2] + rotation_axis[1] * std::sin(rotation_angle);
      rotation_matrix[1][0] = (1-std::cos(rotation_angle)) * rotation_axis[1]*rotation_axis[0] + rotation_axis[2] * std::sin(rotation_angle);
      rotation_matrix[1][1] = (1-std::cos(rotation_angle)) * rotation_axis[1]*rotation_axis[1] + std::cos(rotation_angle);
      rotation_matrix[1][2] = (1-std::cos(rotation_angle)) * rotation_axis[1]*rotation_axis[2] - rotation_axis[0] * std::sin(rotation_angle);
      rotation_matrix[2][0] = (1-std::cos(rotation_angle)) * rotation_axis[2]*rotation_axis[0] - rotation_axis[1] * std::sin(rotation_angle);
      rotation_matrix[2][1] = (1-std::cos(rotation_angle)) * rotation_axis[2]*rotation_axis[1] + rotation_axis[0] * std::sin(rotation_angle);
      rotation_matrix[2][2] = (1-std::cos(rotation_angle)) * rotation_axis[2]*rotation_axis[2] + std::cos(rotation_angle);
      return rotation_matrix;
    }



    Tensor<2,3>
    compute_rotation_matrix_for_slice (const Tensor<1,3> &point_one,
                                       const Tensor<1,3> &point_two)
    {
      AssertThrow(point_one.norm() > std::numeric_limits<double>::min()
                  && point_two.norm() > std::numeric_limits<double>::min(),
                  ExcMessage("The points that are used to define the slice that "
                             "should be rotated onto the x-y-plane can not lie "
                             "at the origin of the coordinate system."));

      // Set up the normal vector of an unrotated 2d spherical shell
      // that by default lies in the x-y plane.
      const Tensor<1,3> unrotated_normal_vector ({0.0,0.0,1.0});

      // Compute the normal vector of the plane that contains
      // the origin and the two points specified as the function arguments.
      Tensor<1,3> rotated_normal_vector = cross_product_3d(point_one, point_two);

      AssertThrow(rotated_normal_vector.norm() > std::numeric_limits<double>::min(),
                  ExcMessage("The points that are used to define the slice that "
                             "should be rotated onto the x-y-plane can not lie "
                             "along the line that also goes through the origin "
                             "of the coordinate system."));

      rotated_normal_vector /= rotated_normal_vector.norm();

      Tensor<2,3> rotation_matrix ({{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}});

      // Calculate the crossing line of the two normals,
      // which will be the rotation axis to transform the one
      // normal into the other
      Tensor<1,3> rotation_axis = cross_product_3d(unrotated_normal_vector, rotated_normal_vector);

      // If the normal vector of the slice already points in z-direction, we do not have to
      // apply the first rotation.
      if (rotation_axis.norm() > std::numeric_limits<double>::min())
        {
          rotation_axis /= rotation_axis.norm();

          // Calculate the rotation angle from the inner product rule
          const double rotation_angle = std::acos(rotated_normal_vector * unrotated_normal_vector);
          rotation_matrix = rotation_matrix_from_axis(rotation_axis, rotation_angle);
        }

      // Now apply the rotation that will project point_one onto the known point
      // (0,1,0).
      const Tensor<1,3> normalized_point_one = point_one / point_one.norm();
      const Tensor<1,3> rotated_point_one = transpose(rotation_matrix) * normalized_point_one;
      const Tensor<1,3> final_point_one ({0.0,1.0,0.0});

      const double second_rotation_angle = std::acos(rotated_point_one * final_point_one);
      Tensor<1,3> second_rotation_axis = cross_product_3d(final_point_one, rotated_point_one);

      // If point 1 is already located at (0,1,0) after the first rotation, we do not
      // have to apply the second rotation.
      if (second_rotation_axis.norm() > std::numeric_limits<double>::min())
        {
          second_rotation_axis /= second_rotation_axis.norm();
          const Tensor<2,3> second_rotation_matrix = rotation_matrix_from_axis(second_rotation_axis, second_rotation_angle);

          // The final rotation used for the model will be the combined
          // rotation of the two operations above. This is achieved by a
          // matrix multiplication of the rotation matrices.
          // This concatenation of rotations is the reason for using a
          // rotation matrix instead of a combined rotation_axis + angle.
          rotation_matrix = rotation_matrix * second_rotation_matrix;
        }
      return rotation_matrix;
    }



    std::pair<double,double> real_spherical_harmonic( const unsigned int l,
                                                      const unsigned int m,
                                                      const double theta,
                                                      const double phi)
    {
      Assert (l<65, ExcMessage("ASPECT uses the implementation of spherical harmonics provided by the "
                               "BOOST library, which is only accurate up a limited degree (see "
                               "https://www.boost.org/doc/libs/1_62_0/libs/math/doc/html/math_toolkit/sf_poly/sph_harm.html#math_toolkit.sf_poly.sph_harm.accuracy). "
                               "The function Utilities::real_spherical_harmonic was asked for a higher degree of "
                               + std::to_string(l) + ", which is not guaranteed to produce a correct value. "
                               "If you absolutely need to use higher degrees, uncomment this assertion on "
                               "your own risk and benchmark your results before using them."));

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
      std::ifstream ifile(filename);

      // return whether construction of the input file has succeeded;
      // success requires the file to exist and to be readable
      return static_cast<bool>(ifile);
    }



    bool
    fexists(const std::string &filename, const MPI_Comm comm)
    {
      bool file_exists = false;
      if (Utilities::MPI::this_mpi_process(comm) == 0)
        {
          std::ifstream ifile(filename);

          // return whether construction of the input file has succeeded;
          // success requires the file to exist and to be readable
          file_exists = static_cast<bool>(ifile);
        }
      return Utilities::MPI::broadcast(comm, file_exists);
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
                                     const MPI_Comm comm)
    {
      std::string data_string;

      if (Utilities::MPI::this_mpi_process(comm) == 0)
        {
          std::size_t filesize;

          // Check to see if the prm file will be reading data from disk or
          // from a provided URL
          if (filename_is_url(filename))
            {
#ifdef ASPECT_WITH_LIBDAP
              std::unique_ptr<libdap::Connect> url
                = std::make_unique<libdap::Connect>(filename);
              libdap::BaseTypeFactory factory;
              libdap::DataDDS dds(&factory);
              libdap::DAS das;

              url->request_data(dds, "");
              url->request_das(das);


              // Temporary vector that will hold the different arrays stored in urlArray
              std::vector<libdap::dods_float32> tmp;
              // Vector that will hold the arrays (columns) and the values within those arrays
              std::vector<std::vector<libdap::dods_float32>> columns;

              // Check dds values to make sure the arrays are of the same length and of type string
              for (libdap::DDS::Vars_iter i = dds.var_begin(); i != dds.var_end(); ++i)
                {
                  libdap::BaseType *btp = *i;
                  if ((*i)->type() == libdap::dods_array_c)
                    {
                      // Array to store the url data
                      libdap::Array *urlArray;
                      urlArray = static_cast <libdap::Array *>(btp);
                      if (urlArray->var() != nullptr && urlArray->var()->type() == libdap::dods_float32_c)
                        {
                          tmp.resize(urlArray->length());

                          // The url Array contains a separate array for each column of data.
                          // This will put each of these individual arrays into its own vector.
                          urlArray->value(&tmp[0]);
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
              for (libdap::AttrTable::Attr_iter i = das.var_begin(); i != das.var_end(); ++i)
                {
                  libdap::AttrTable *table = das.get_table(i);
                  if (table->get_attr("POINTS") != "")
                    points.push_back(table->get_attr("POINTS"));
                  if (table->get_attr("points") != "")
                    points.push_back(table->get_attr("points"));
                }

              std::stringstream urlString;

              // Append the gathered POINTS in the proper format:
              // "# POINTS: <val1> <val2> <val3>"
              urlString << "# POINTS:";
              for (unsigned int i = 0; i < points.size(); ++i)
                {
                  urlString << ' ' << points[i];
                }
              urlString << "\n";

              // Add the values from the arrays into the stringstream. The values are passed in
              // per row with a character return added at the end of each row.
              // TODO: Add a check to make sure that each column is the same size before writing
              //     to the stringstream
              for (unsigned int i = 0; i < tmp.size(); ++i)
                {
                  for (unsigned int j = 0; j < columns.size(); ++j)
                    {
                      urlString << columns[j][i];
                      urlString << ' ';
                    }
                  urlString << "\n";
                }

              data_string = urlString.str();
              filesize = data_string.size();

#else // ASPECT_WITH_LIBDAP

              // Broadcast failure state, then throw. We signal the failure by
              // setting the file size to an invalid size, then trigger an assert.
              {
                std::size_t invalid_filesize = numbers::invalid_size_type;
                const int ierr = MPI_Bcast(&invalid_filesize, 1, Utilities::internal::MPI::mpi_type_id(&filesize), 0, comm);
                AssertThrowMPI(ierr);
              }
              AssertThrow(false,
                          ExcMessage(std::string("Reading of file ") + filename + " failed. " +
                                     "Make sure you have the dependencies for reading a url " +
                                     "(run cmake with -DASPECT_WITH_LIBDAP=ON)"));

#endif // ASPECT_WITH_LIBDAP
            }
          else
            {
              std::ifstream filestream;
              const bool filename_ends_in_gz = std::regex_search(filename, std::regex("\\.gz$"));
              if (filename_ends_in_gz == true)
                filestream.open(filename, std::ios_base::in | std::ios_base::binary);
              else
                filestream.open(filename);

              if (!filestream)
                {
                  // broadcast failure state, then throw
                  std::size_t invalid_filesize = numbers::invalid_size_type;
                  const int ierr = MPI_Bcast(&invalid_filesize, 1, Utilities::internal::MPI::mpi_type_id(&filesize), 0, comm);
                  AssertThrowMPI(ierr);
                  AssertThrow (false,
                               ExcMessage (std::string("Could not open file <") + filename + ">."));
                }

              // Read data from disk
              std::stringstream datastream;

              try
                {
                  boost::iostreams::filtering_istreambuf in;
                  if (filename_ends_in_gz == true)
                    in.push(boost::iostreams::gzip_decompressor());

                  in.push(filestream);
                  boost::iostreams::copy(in, datastream);
                }
              catch (const std::ios::failure &)
                {
                  // broadcast failure state, then throw
                  std::size_t invalid_filesize = numbers::invalid_size_type;
                  const int ierr = MPI_Bcast(&invalid_filesize, 1, Utilities::internal::MPI::mpi_type_id(&filesize), 0, comm);
                  AssertThrowMPI(ierr);
                  AssertThrow (false,
                               ExcMessage (std::string("Could not read file content from <") + filename + ">."));
                }

              data_string = datastream.str();
              filesize = data_string.size();
            }

          // Distribute data_size and data across processes
          int ierr = MPI_Bcast(&filesize, 1, Utilities::internal::MPI::mpi_type_id(&filesize), 0, comm);
          AssertThrowMPI(ierr);

          big_mpi::broadcast(&data_string[0], filesize, 0, comm);
        }
      else
        {
          // Prepare for receiving data
          std::size_t filesize;
          int ierr = MPI_Bcast(&filesize, 1, Utilities::internal::MPI::mpi_type_id(&filesize), 0, comm);
          AssertThrowMPI(ierr);
          if (filesize == numbers::invalid_size_type)
            throw QuietException();

          data_string.resize(filesize);

          // Receive and store data
          big_mpi::broadcast(&data_string[0], filesize, 0, comm);
        }

      return data_string;
    }



    void
    collect_and_write_file_content(const std::string &filename,
                                   const std::string &file_content,
                                   const MPI_Comm comm)
    {
      const std::vector<std::string> collected_content = Utilities::MPI::gather(comm, file_content);

      if (Utilities::MPI::this_mpi_process(comm) == 0)
        {
          std::ofstream filestream(filename);

          AssertThrow (filestream.good(),
                       ExcMessage (std::string("Could not open file <") + filename + ">."));

          try
            {
              for (const auto &content : collected_content)
                filestream << content;

              bool success = filestream.good();
              const int ierr = MPI_Bcast(&success, 1, Utilities::internal::MPI::mpi_type_id(&success), 0, comm);
              AssertThrowMPI(ierr);
            }
          catch (const std::ios::failure &)
            {
              // broadcast failure state, then throw
              bool success = false;
              const int ierr = MPI_Bcast(&success, 1, Utilities::internal::MPI::mpi_type_id(&success), 0, comm);
              AssertThrowMPI(ierr);
              AssertThrow (false,
                           ExcMessage (std::string("Could not write content to file <") + filename + ">."));
            }

          filestream.close();
        }
      else
        {
          // Check that the file was written successfully
          bool success;
          int ierr = MPI_Bcast(&success, 1, Utilities::internal::MPI::mpi_type_id(&success), 0, comm);
          AssertThrowMPI(ierr);
          if (success == false)
            throw QuietException();
        }
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
          const std::string subdir = pathname.substr(0,++pos);
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
                          const MPI_Comm comm,
                          const bool silent)
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
          std::vector<std::vector<double>> m_upper;

          /**
           * diagonals below the diagonal
           */
          std::vector<std::vector<double>> m_lower;
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
        for (auto &x : m_upper)
          x.resize(dim);
        for (auto &x : m_lower)
          x.resize(dim);
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
        for (int i = 0; i < this->dim(); ++i)
          {
            assert(this->operator()(i,i) != 0.0);
            this->saved_diag(i) = 1.0/this->operator()(i,i);
            j_min = std::max(0,i-this->num_lower());
            j_max = std::min(this->dim()-1,i+this->num_upper());
            for (int j = j_min; j <= j_max; ++j)
              {
                this->operator()(i,j) *= this->saved_diag(i);
              }
            this->operator()(i,i) = 1.0;          // prevents rounding errors
          }

        // Gauss LR-Decomposition
        for (int k = 0; k < this->dim(); ++k)
          {
            i_max = std::min(this->dim()-1,k+this->num_lower());  // num_lower not a mistake!
            for (int i = k+1; i <= i_max; ++i)
              {
                assert(this->operator()(k,k) != 0.0);
                x = -this->operator()(i,k)/this->operator()(k,k);
                this->operator()(i,k) = -x;                         // assembly part of L
                j_max = std::min(this->dim()-1, k + this->num_upper());
                for (int j = k+1; j <= j_max; ++j)
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
        for (int i = 0; i < this->dim(); ++i)
          {
            sum = 0;
            j_start = std::max(0,i-this->num_lower());
            for (int j = j_start; j < i; ++j) sum += this->operator()(i,j)*x[j];
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
            for (int j = i+1; j <= j_stop; ++j) sum += this->operator()(i,j)*x[j];
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
        for (unsigned int i = 0; i < n-1; ++i)
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
                for (unsigned int i=0; i < n-1; ++i)
                  {
                    dxs[i] = x[i+1]-x[i];
                    dys[i] = y[i+1]-y[i];
                    ms[i] = dys[i]/dxs[i];
                  }

                // get m_a parameter
                m_c.resize(n);
                m_c[0] = 0;

                for (unsigned int i = 0; i < n-2; ++i)
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
                for (unsigned int i = 0; i < m_c.size()-1; ++i)
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
                for (unsigned int i = 1; i<n-1; ++i)
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
                for (unsigned int i = 0; i<n-1; ++i)
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
            for (unsigned int i = 0; i<n-1; ++i)
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
        const int idx = std::max( static_cast<int>(it-m_x.begin())-1, 0);

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
      // Check for environment variable override to ASPECT_SOURCE_DIR
      char const *ASPECT_SOURCE_DIR_env = getenv("ASPECT_SOURCE_DIR");
      if (ASPECT_SOURCE_DIR_env != nullptr)
        {
          return Utilities::replace_in_string(location,
                                              "$ASPECT_SOURCE_DIR",
                                              ASPECT_SOURCE_DIR_env);
        }

      // Otherwise, use the default define from config.h
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
    string_to_bool(const std::string &s)
    {
      return (s == "true" || s == "yes");
    }



    std::vector<bool>
    string_to_bool(const std::vector<std::string> &s)
    {
      std::vector<bool> result;
      result.reserve(s.size());

      for (auto &i : s)
        result.push_back(string_to_bool(i));

      return result;
    }



    unsigned int
    string_to_unsigned_int(const std::string &s)
    {
      const int value = dealii::Utilities::string_to_int(s);
      AssertThrow (value >= 0, ExcMessage("Negative number in string_to_unsigned_int() detected."));
      return static_cast<unsigned int>(value);
    }



    std::vector<unsigned int>
    string_to_unsigned_int(const std::vector<std::string> &s)
    {
      std::vector<unsigned int> result;
      result.reserve(s.size());

      for (auto &str : s)
        result.emplace_back(string_to_unsigned_int(str));

      return result;
    }



    bool
    has_unique_entries (const std::vector<std::string> &strings)
    {
      const std::set<std::string> set_of_strings(strings.begin(),strings.end());
      return (set_of_strings.size() == strings.size());
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
          return Utilities::fixed_power<-2>(averaged_parameter_derivative_part_1/sum_of_weights) * averaged_parameter_derivative_part_2/sum_of_weights;
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


      // The factor in the Newton matrix is going to be of the form
      //   2*eta I + (a \otimes b + b \otimes a)
      // where a=strain_rate and b=dviscosities_dstrain_rate.
      //
      // If a,b are parallel, this simplifies to
      //   [2*eta + 2 a:b] I =  2 [eta + a:b] I
      // and we need to make sure that
      //   [eta + alpha a:b] > (1-safety_factor)*eta
      // by choosing alpha appropriately.

      // So, first check: If
      //   [eta + a:b] > (1-safety_factor)*eta
      // is already satisfied, then we can choose alpha=1
      const double a_colon_b = strain_rate * dviscosities_dstrain_rate;
      if (eta + a_colon_b > eta * (1. - SPD_safety_factor))
        return 1.0;
      else
        {
          // Otherwise solve the equation above for alpha, which yields
          //   a:b = -safety_factor*eta / a:b
          // This can only ever happen if a:b < 0, so we get
          //   a:b = safety_factor * abs(eta / a:b)
          Assert (a_colon_b < 0, ExcInternalError());
          return SPD_safety_factor * std::abs(eta / a_colon_b);
        }
    }



    template <int dim>
    Point<dim> convert_array_to_point(const std::array<double,dim> &array)
    {
      Point<dim> point;
      for (unsigned int i = 0; i < dim; ++i)
        point[i] = array[i];

      return point;
    }



    template <int dim>
    std::array<double,dim> convert_point_to_array(const Point<dim> &point)
    {
      std::array<double,dim> array;
      for (unsigned int i = 0; i < dim; ++i)
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



    template <int dim>
    const std::array<double,dim> &NaturalCoordinate<dim>::get_coordinates() const
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
                       const std::vector<Point<dim>> &,
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
            for (unsigned int q=0; q<n_q_points; ++q)
              {
                const double JxW = fe_values.JxW(q);
                for (unsigned int i=0; i<dofs_per_cell; ++i)
                  {
                    if (dof_handler.get_fe().system_to_component_index(i).first == component_index)
                      cell_vector(i) +=
                        rhs_values[q] *
                        fe_values[extractor].value(i,q) *
                        JxW;

                    for (unsigned int j=0; j<dofs_per_cell; ++j)
                      if ((dof_handler.get_fe().system_to_component_index(i).first ==
                           component_index)
                          &&
                          (dof_handler.get_fe().system_to_component_index(j).first ==
                           component_index))
                        local_mass_matrix(j,i) += (fe_values[extractor].value(i,q) * fe_values[extractor].value(j,q) *
                                                   JxW);
                      else if (i == j)
                        local_mass_matrix(i,j) = 1.;
                  }
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



    void throw_linear_solver_failure_exception(const std::string &solver_name,
                                               const std::string &function_name,
                                               const std::vector<SolverControl> &solver_controls,
                                               const std::exception &exc,
                                               const MPI_Comm mpi_communicator,
                                               const std::string &output_filename)
    {
      if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
        {
          std::ostringstream exception_message;
          exception_message << std::scientific
                            << "The " + solver_name
                            << " in " + function_name
                            << " did not converge. " << std::endl << std::endl;

          if (solver_controls.front().last_step() != numbers::invalid_unsigned_int)
            exception_message << "The initial residual was: "
                              << solver_controls.front().initial_value() << std::endl;

          if (solver_controls.back().last_step() != numbers::invalid_unsigned_int)
            exception_message << "The final residual is: "
                              << solver_controls.back().last_value() << std::endl;

          exception_message << "The required residual for convergence is: "
                            << solver_controls.front().tolerance() << std::endl;

          if (output_filename != "")
            {
              // output solver history
              std::ofstream f((output_filename));

              for (const auto &solver_control: solver_controls)
                {
                  // Skip the output if no iterations were run for this solver
                  if (solver_control.last_step() == numbers::invalid_unsigned_int)
                    continue;

                  // Add an empty line between solvers
                  if (&solver_control != &(solver_controls.front()))
                    f << std::endl;

                  unsigned int j=0;
                  for (const auto &residual: solver_control.get_history_data())
                    f << j++ << ' ' << residual << std::endl;
                }

              f.close();

              exception_message << "See " << output_filename
                                << " for the full convergence history." << std::endl;
            }

          exception_message << std::endl
                            << "The solver reported the following error:"
                            << std::endl;
          exception_message << exc.what();

          AssertThrow (false,
                       ExcMessage (exception_message.str()));
        }
      else
        throw QuietException();
    }



    std::vector<Tensor<2,3>>
    rotation_matrices_random_draw_volume_weighting(const std::vector<double> &volume_fraction,
                                                   const std::vector<Tensor<2,3>> &rotation_matrices,
                                                   const unsigned int n_output_matrices,
                                                   std::mt19937 &random_number_generator)
    {
      const unsigned int n_grains = volume_fraction.size();

      // Get volume weighted euler angles, using random draws to convert odf
      // to a discrete number of orientations, weighted by volume
      // 1a. Sort the volume fractions and matrices based on the volume fractions size
      const auto p = compute_sorting_permutation(volume_fraction);

      const std::vector<double> fv_sorted = apply_permutation(volume_fraction, p);
      const std::vector<Tensor<2,3>> matrices_sorted = apply_permutation(rotation_matrices, p);

      // 2. Get cumulative weight for volume fraction
      std::vector<double> cum_weight(n_grains);
      std::partial_sum(fv_sorted.begin(),fv_sorted.end(),cum_weight.begin());

      // 3. Generate random indices
      std::uniform_real_distribution<> dist(0, cum_weight[n_grains-1]);
      std::vector<double> idxgrain(n_output_matrices);
      for (unsigned int grain_i = 0; grain_i < n_output_matrices; ++grain_i)
        {
          idxgrain[grain_i] = dist(random_number_generator);
        }

      // 4. Find the maximum cum_weight that is less than the random value.
      // the euler angle index is +1. For example, if the idxGrain(g) < cumWeight(1),
      // the index should be 1 not zero)
      std::vector<Tensor<2,3>> matrices_out(n_output_matrices);
      for (unsigned int grain_i = 0; grain_i < n_output_matrices; ++grain_i)
        {
          const std::vector<double>::iterator selected_matrix =
            std::lower_bound(cum_weight.begin(),
                             cum_weight.end(),
                             idxgrain[grain_i]);

          const unsigned int matrix_index =
            std::distance(cum_weight.begin(), selected_matrix);

          matrices_out[grain_i] = matrices_sorted[matrix_index];
        }
      return matrices_out;
    }



    double
    wrap_angle(const double angle)
    {
      return angle - 360.0*std::floor(angle/360.0);
    }



    std::array<double,3>
    zxz_euler_angles_from_rotation_matrix(const Tensor<2,3> &rotation_matrix)
    {
      // ZXZ Euler angles
      std::array<double,3> euler_angles;
      for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
          Assert(std::abs(rotation_matrix[i][j]) <= 1.0,
                 ExcMessage("rotation_matrix[" + std::to_string(i) + "][" + std::to_string(j) +
                            "] is larger than one: " + std::to_string(rotation_matrix[i][j]) + " (" + std::to_string(rotation_matrix[i][j]-1.0) + "). rotation_matrix = \n"
                            + std::to_string(rotation_matrix[0][0]) + " " + std::to_string(rotation_matrix[0][1]) + " " + std::to_string(rotation_matrix[0][2]) + "\n"
                            + std::to_string(rotation_matrix[1][0]) + " " + std::to_string(rotation_matrix[1][1]) + " " + std::to_string(rotation_matrix[1][2]) + "\n"
                            + std::to_string(rotation_matrix[2][0]) + " " + std::to_string(rotation_matrix[2][1]) + " " + std::to_string(rotation_matrix[2][2])));


      AssertThrow(rotation_matrix[2][2] <= 1.0, ExcMessage("rot_matrix[2][2] > 1.0"));

      const double theta = std::acos(rotation_matrix[2][2]);
      double phi1 = 0.0;
      double phi2 = 0.0;

      if (theta != 0.0 && theta != dealii::numbers::PI)
        {
          //
          phi1  = std::atan2(rotation_matrix[2][0]/-std::sin(theta),rotation_matrix[2][1]/-std::sin(theta));
          phi2  = std::atan2(rotation_matrix[0][2]/-std::sin(theta),rotation_matrix[1][2]/std::sin(theta));
        }

      // note that in the case theta is 0 or phi a dimension is lost
      // see: https://en.wikipedia.org/wiki/Gimbal_lock. We set phi1
      // to 0 and compute the corresponding phi2. The resulting direction
      // (cosine matrix) should be the same.
      else if (theta == 0.0)
        {
          phi2 = - phi1 - std::atan2(rotation_matrix[0][1],rotation_matrix[0][0]);
        }
      else
        {
          phi2 = phi1 + std::atan2(rotation_matrix[0][1],rotation_matrix[0][0]);
        }


      AssertThrow(!std::isnan(phi1), ExcMessage("phi1 is not a number. theta = " + std::to_string(theta) + ", rotation_matrix[2][2]= " + std::to_string(rotation_matrix[2][2])
                                                + ", acos(rotation_matrix[2][2]) = " + std::to_string(std::acos(rotation_matrix[2][2])) + ", acos(1.0) = " + std::to_string(std::acos(1.0))));
      AssertThrow(!std::isnan(theta), ExcMessage("theta is not a number."));
      AssertThrow(!std::isnan(phi2), ExcMessage("phi2 is not a number."));

      euler_angles[0] = wrap_angle(phi1 * constants::radians_to_degree);
      euler_angles[1] = wrap_angle(theta * constants::radians_to_degree);
      euler_angles[2] = wrap_angle(phi2 * constants::radians_to_degree);


      AssertThrow(!std::isnan(euler_angles[0]), ExcMessage(" euler_angles[0] is not a number."));
      AssertThrow(!std::isnan(euler_angles[1]), ExcMessage(" euler_angles[1] is not a number."));
      AssertThrow(!std::isnan(euler_angles[2]), ExcMessage(" euler_angles[2] is not a number."));

      return euler_angles;
    }



    Tensor<2,3>
    zxz_euler_angles_to_rotation_matrix(const double phi1_degrees, const double theta_degrees, const double phi2_degrees)
    {
      // ZXZ Euler angles
      const double phi1 = phi1_degrees * constants::degree_to_radians;
      const double theta = theta_degrees * constants::degree_to_radians;
      const double phi2 = phi2_degrees * constants::degree_to_radians;
      Tensor<2,3> rot_matrix;

      rot_matrix[0][0] = std::cos(phi2)*std::cos(phi1) - std::cos(theta)*std::sin(phi1)*std::sin(phi2);
      rot_matrix[0][1] = -std::cos(phi2)*std::sin(phi1) - std::cos(theta)*std::cos(phi1)*std::sin(phi2);
      rot_matrix[0][2] = -std::sin(phi2)*std::sin(theta);

      rot_matrix[1][0] = std::sin(phi2)*std::cos(phi1) + std::cos(theta)*std::sin(phi1)*std::cos(phi2);
      rot_matrix[1][1] = -std::sin(phi2)*std::sin(phi1) + std::cos(theta)*std::cos(phi1)*std::cos(phi2);
      rot_matrix[1][2] = std::cos(phi2)*std::sin(theta);

      rot_matrix[2][0] = -std::sin(theta)*std::sin(phi1);
      rot_matrix[2][1] = -std::sin(theta)*std::cos(phi1);
      rot_matrix[2][2] = std::cos(theta);
      AssertThrow(rot_matrix[2][2] <= 1.0, ExcMessage("rot_matrix[2][2] > 1.0"));

      return rot_matrix;
    }



    namespace Tensors
    {
      SymmetricTensor<4,3>
      rotate_full_stiffness_tensor(const Tensor<2,3> &rotation_tensor, const SymmetricTensor<4,3> &input_tensor)
      {
        SymmetricTensor<4,3> output;

        // Dealii symmetric tensor is C_{ijkl} == C_{jikl} == C_{ijlk}, but not C_{klij}.
        // So we make sure that those entries are not added twice in this loop by having
        // the second and 4th loop starting with the first and third index respectively.
        for (unsigned short int i1 = 0; i1 < 3; ++i1)
          {
            for (unsigned short int i2 = i1; i2 < 3; ++i2)
              {
                for (unsigned short int i3 = 0; i3 < 3; ++i3)
                  {
                    for (unsigned short int i4 = i3; i4 < 3; ++i4)
                      {
                        for (unsigned short int j1 = 0; j1 < 3; ++j1)
                          {
                            for (unsigned short int j2 = 0; j2 < 3; ++j2)
                              {
                                for (unsigned short int j3 = 0; j3 < 3; ++j3)
                                  {
                                    for (unsigned short int j4 = 0; j4 < 3; ++j4)
                                      {
                                        output[i1][i2][i3][i4] += rotation_tensor[i1][j1]*rotation_tensor[i2][j2]*rotation_tensor[i3][j3]*rotation_tensor[i4][j4]*input_tensor[j1][j2][j3][j4];
                                      }
                                  }
                              }
                          }
                      }
                  }
              }
          }

        return output;
      }



      SymmetricTensor<2,6>
      rotate_voigt_stiffness_matrix(const Tensor<2,3> &rotation_tensor, const SymmetricTensor<2,6> &input_tensor)
      {
        // we can represent the rotation of the 4th order tensor as a rotation in the Voigt
        // notation by computing $C'=MCM^{-1}$. Because M is orthogonal we can replace $M^{-1}$
        // with $M^T$ resulting in $C'=MCM^{T}$ (Carcione, J. M. (2007). Wave Fields in Real Media:
        // Wave Propagation in Anisotropic, Anelastic, Porous and Electromagnetic Media. Netherlands:
        // Elsevier Science. Pages 8-9).
        Tensor<2,6> rotation_matrix;
        // top left block
        rotation_matrix[0][0] = rotation_tensor[0][0] * rotation_tensor[0][0];
        rotation_matrix[1][0] = rotation_tensor[1][0] * rotation_tensor[1][0];
        rotation_matrix[2][0] = rotation_tensor[2][0] * rotation_tensor[2][0];
        rotation_matrix[0][1] = rotation_tensor[0][1] * rotation_tensor[0][1];
        rotation_matrix[1][1] = rotation_tensor[1][1] * rotation_tensor[1][1];
        rotation_matrix[2][1] = rotation_tensor[2][1] * rotation_tensor[2][1];
        rotation_matrix[0][2] = rotation_tensor[0][2] * rotation_tensor[0][2];
        rotation_matrix[1][2] = rotation_tensor[1][2] * rotation_tensor[1][2];
        rotation_matrix[2][2] = rotation_tensor[2][2] * rotation_tensor[2][2];

        // top right block
        rotation_matrix[0][3] = 2.0 * rotation_tensor[0][1] * rotation_tensor[0][2];
        rotation_matrix[1][3] = 2.0 * rotation_tensor[1][1] * rotation_tensor[1][2];
        rotation_matrix[2][3] = 2.0 * rotation_tensor[2][1] * rotation_tensor[2][2];
        rotation_matrix[0][4] = 2.0 * rotation_tensor[0][2] * rotation_tensor[0][0];
        rotation_matrix[1][4] = 2.0 * rotation_tensor[1][2] * rotation_tensor[1][0];
        rotation_matrix[2][4] = 2.0 * rotation_tensor[2][2] * rotation_tensor[2][0];
        rotation_matrix[0][5] = 2.0 * rotation_tensor[0][0] * rotation_tensor[0][1];
        rotation_matrix[1][5] = 2.0 * rotation_tensor[1][0] * rotation_tensor[1][1];
        rotation_matrix[2][5] = 2.0 * rotation_tensor[2][0] * rotation_tensor[2][1];

        // bottom left block
        rotation_matrix[3][0] = rotation_tensor[1][0] * rotation_tensor[2][0];
        rotation_matrix[4][0] = rotation_tensor[2][0] * rotation_tensor[0][0];
        rotation_matrix[5][0] = rotation_tensor[0][0] * rotation_tensor[1][0];
        rotation_matrix[3][1] = rotation_tensor[1][1] * rotation_tensor[2][1];
        rotation_matrix[4][1] = rotation_tensor[2][1] * rotation_tensor[0][1];
        rotation_matrix[5][1] = rotation_tensor[0][1] * rotation_tensor[1][1];
        rotation_matrix[3][2] = rotation_tensor[1][2] * rotation_tensor[2][2];
        rotation_matrix[4][2] = rotation_tensor[2][2] * rotation_tensor[0][2];
        rotation_matrix[5][2] = rotation_tensor[0][2] * rotation_tensor[1][2];

        // bottom right block
        rotation_matrix[3][3] = rotation_tensor[1][1] * rotation_tensor[2][2] + rotation_tensor[1][2] * rotation_tensor[2][1];
        rotation_matrix[4][3] = rotation_tensor[0][1] * rotation_tensor[2][2] + rotation_tensor[0][2] * rotation_tensor[2][1];
        rotation_matrix[5][3] = rotation_tensor[0][1] * rotation_tensor[1][2] + rotation_tensor[0][2] * rotation_tensor[1][1];
        rotation_matrix[3][4] = rotation_tensor[1][0] * rotation_tensor[2][2] + rotation_tensor[1][2] * rotation_tensor[2][0];
        rotation_matrix[4][4] = rotation_tensor[0][2] * rotation_tensor[2][0] + rotation_tensor[0][0] * rotation_tensor[2][2];
        rotation_matrix[5][4] = rotation_tensor[0][2] * rotation_tensor[1][0] + rotation_tensor[0][0] * rotation_tensor[1][2];
        rotation_matrix[3][5] = rotation_tensor[1][1] * rotation_tensor[2][0] + rotation_tensor[1][0] * rotation_tensor[2][1];
        rotation_matrix[4][5] = rotation_tensor[0][0] * rotation_tensor[2][1] + rotation_tensor[0][1] * rotation_tensor[2][0];
        rotation_matrix[5][5] = rotation_tensor[0][0] * rotation_tensor[1][1] + rotation_tensor[0][1] * rotation_tensor[1][0];

        const Tensor<2,6> rotation_matrix_transposed = transpose(rotation_matrix);

        return symmetrize((rotation_matrix*input_tensor)*rotation_matrix_transposed);
      }



      SymmetricTensor<2,6>
      to_voigt_stiffness_matrix(const SymmetricTensor<4,3> &input_tensor)
      {
        SymmetricTensor<2,6> output;

        for (unsigned short int i = 0; i < 3; i++)
          {
            output[i][i] = input_tensor[i][i][i][i];
          }

        for (unsigned short int i = 1; i < 3; i++)
          {
            output[0][i] = 0.5*(input_tensor[0][0][i][i] + input_tensor[i][i][0][0]);
            //symmetry: output[0][i] = output[i][0];
          }
        output[1][2]=0.5*(input_tensor[1][1][2][2]+input_tensor[2][2][1][1]);
        //symmetry: output[2][1]=output[1][2];

        for (unsigned short int i = 0; i < 3; i++)
          {
            output[i][3]=0.25*(input_tensor[i][i][1][2]+input_tensor[i][i][2][1]+ input_tensor[1][2][i][i]+input_tensor[2][1][i][i]);
            //symmetry: output[3][i]=output[i][3];
          }

        for (unsigned short int i = 0; i < 3; i++)
          {
            output[i][4]=0.25*(input_tensor[i][i][0][2]+input_tensor[i][i][2][0]+ input_tensor[0][2][i][i]+input_tensor[2][0][i][i]);
            //symmetry: output[4][i]=output[i][4];
          }

        for (unsigned short int i = 0; i < 3; i++)
          {
            output[i][5]=0.25*(input_tensor[i][i][0][1]+input_tensor[i][i][1][0]+input_tensor[0][1][i][i]+input_tensor[1][0][i][i]);
            //symmetry: output[5][i]=output[i][5];
          }

        output[3][3]=0.25*(input_tensor[1][2][1][2]+input_tensor[1][2][2][1]+input_tensor[2][1][1][2]+input_tensor[2][1][2][1]);
        output[4][4]=0.25*(input_tensor[0][2][0][2]+input_tensor[0][2][2][0]+input_tensor[2][0][0][2]+input_tensor[2][0][2][0]);
        output[5][5]=0.25*(input_tensor[1][0][1][0]+input_tensor[1][0][0][1]+input_tensor[0][1][1][0]+input_tensor[0][1][0][1]);

        output[3][4]=0.125*(input_tensor[1][2][0][2]+input_tensor[1][2][2][0]+input_tensor[2][1][0][2]+input_tensor[2][1][2][0]+input_tensor[0][2][1][2]+input_tensor[0][2][2][1]+input_tensor[2][0][1][2]+input_tensor[2][0][2][1]);
        //symmetry: output[4][3]=output[3][4];
        output[3][5]=0.125*(input_tensor[1][2][0][1]+input_tensor[1][2][1][0]+input_tensor[2][1][0][1]+input_tensor[2][1][1][0]+input_tensor[0][1][1][2]+input_tensor[0][1][2][1]+input_tensor[1][0][1][2]+input_tensor[1][0][2][1]);
        //symmetry: output[5][3]=output[3][5];
        output[4][5]=0.125*(input_tensor[0][2][0][1]+input_tensor[0][2][1][0]+input_tensor[2][0][0][1]+input_tensor[2][0][1][0]+input_tensor[0][1][0][2]+input_tensor[0][1][2][0]+input_tensor[1][0][0][2]+input_tensor[1][0][2][0]);
        //symmetry: output[5][4]=output[4][5];

        return output;
      }



      SymmetricTensor<4,3>
      to_full_stiffness_tensor(const SymmetricTensor<2,6> &input_tensor)
      {
        SymmetricTensor<4,3> output;

        for (unsigned short int i = 0; i < 3; i++)
          for (unsigned short int j = 0; j < 3; j++)
            for (unsigned short int k = 0; k < 3; k++)
              for (unsigned short int l = 0; l < 3; l++)
                {
                  // The first part of the inline if statement gets the diagonal.
                  // The second part is never higher than 5 (which is the limit of the tensor index)
                  // because to reach this part the variables need to be different, which results in
                  // at least a minus 1.
                  const unsigned short int p = (i == j ? i : 6 - i - j);
                  const unsigned short int q = (k == l ? k : 6 - k - l);
                  output[i][j][k][l] = input_tensor[p][q];
                }
        return output;
      }



      Tensor<1,21>
      to_voigt_stiffness_vector(const SymmetricTensor<2,6> &input)
      {
        return Tensor<1,21,double> (
        {
          input[0][0],           // 0  // 1
          input[1][1],           // 1  // 2
          input[2][2],           // 2  // 3
          numbers::SQRT2 *input[1][2],  // 3  // 4
          numbers::SQRT2 *input[0][2],  // 4  // 5
          numbers::SQRT2 *input[0][1],  // 5  // 6
          2*input[3][3],         // 6  // 7
          2*input[4][4],         // 7  // 8
          2*input[5][5],         // 8  // 9
          2*input[0][3],         // 9  // 10
          2*input[1][4],         // 10 // 11
          2*input[2][5],         // 11 // 12
          2*input[2][3],         // 12 // 13
          2*input[0][4],         // 13 // 14
          2*input[1][5],         // 14 // 15
          2*input[1][3],         // 15 // 16
          2*input[2][4],         // 16 // 17
          2*input[0][5],         // 17 // 18
          2*numbers::SQRT2 *input[4][5], // 18 // 19
          2*numbers::SQRT2 *input[3][5], // 19 // 20
          2*numbers::SQRT2 *input[3][4] // 20 // 21
        });

      }



      SymmetricTensor<2,6>
      to_voigt_stiffness_matrix(const Tensor<1,21> &input)
      {
        SymmetricTensor<2,6> result;

        const double sqrt_2_inv = 1/numbers::SQRT2;

        result[0][0] = input[0];
        result[1][1] = input[1];
        result[2][2] = input[2];
        result[1][2] = sqrt_2_inv * input[3];
        result[0][2] = sqrt_2_inv * input[4];
        result[0][1] = sqrt_2_inv * input[5];
        result[3][3] = 0.5 * input[6];
        result[4][4] = 0.5 * input[7];
        result[5][5] = 0.5 * input[8];
        result[0][3] = 0.5 * input[9];
        result[1][4] = 0.5 * input[10];
        result[2][5] = 0.5 * input[11];
        result[2][3] = 0.5 * input[12];
        result[0][4] = 0.5 * input[13];
        result[1][5] = 0.5 * input[14];
        result[1][3] = 0.5 * input[15];
        result[2][4] = 0.5 * input[16];
        result[0][5] = 0.5 * input[17];
        result[4][5] = 0.5 * sqrt_2_inv * input[18];
        result[3][5] = 0.5 * sqrt_2_inv * input[19];
        result[3][4] = 0.5 * sqrt_2_inv * input[20];

        return result;

      }



      Tensor<1,21>
      to_voigt_stiffness_vector(const SymmetricTensor<4,3> &input_tensor)
      {
        return Tensor<1,21,double> (
        {
          input_tensor[0][0][0][0],           // 0  // 1
          input_tensor[1][1][1][1],           // 1  // 2
          input_tensor[2][2][2][2],           // 2  // 3
          numbers::SQRT2*0.5*(input_tensor[1][1][2][2] + input_tensor[2][2][1][1]),   // 3  // 4
          numbers::SQRT2*0.5*(input_tensor[0][0][2][2] + input_tensor[2][2][0][0]),   // 4  // 5
          numbers::SQRT2*0.5*(input_tensor[0][0][1][1] + input_tensor[1][1][0][0]),   // 5  // 6
          0.5*(input_tensor[1][2][1][2]+input_tensor[1][2][2][1]+input_tensor[2][1][1][2]+input_tensor[2][1][2][1]),         // 6  // 7
          0.5*(input_tensor[0][2][0][2]+input_tensor[0][2][2][0]+input_tensor[2][0][0][2]+input_tensor[2][0][2][0]),         // 7  // 8
          0.5*(input_tensor[1][0][1][0]+input_tensor[1][0][0][1]+input_tensor[0][1][1][0]+input_tensor[0][1][0][1]),         // 8  // 9
          0.5*(input_tensor[0][0][1][2]+input_tensor[0][0][2][1]+input_tensor[1][2][0][0]+input_tensor[2][1][0][0]),         // 9  // 10
          0.5*(input_tensor[1][1][0][2]+input_tensor[1][1][2][0]+input_tensor[0][2][1][1]+input_tensor[2][0][1][1]),         // 10 // 11
          0.5*(input_tensor[2][2][0][1]+input_tensor[2][2][1][0]+input_tensor[0][1][2][2]+input_tensor[1][0][2][2]),         // 11 // 12
          0.5*(input_tensor[2][2][1][2]+input_tensor[2][2][2][1]+input_tensor[1][2][2][2]+input_tensor[2][1][2][2]),         // 12 // 13
          0.5*(input_tensor[0][0][0][2]+input_tensor[0][0][2][0]+input_tensor[0][2][0][0]+input_tensor[2][0][0][0]),         // 13 // 14
          0.5*(input_tensor[1][1][0][1]+input_tensor[1][1][1][0]+input_tensor[0][1][1][1]+input_tensor[1][0][1][1]),         // 14 // 15
          0.5*(input_tensor[1][1][1][2]+input_tensor[1][1][2][1]+input_tensor[1][2][1][1]+input_tensor[2][1][1][1]),         // 15 // 16
          0.5*(input_tensor[2][2][0][2]+input_tensor[2][2][2][0]+input_tensor[0][2][2][2]+input_tensor[2][0][2][2]),         // 16 // 17
          0.5*(input_tensor[0][0][0][1]+input_tensor[0][0][1][0]+input_tensor[0][1][0][0]+input_tensor[1][0][0][0]),         // 17 // 18
          numbers::SQRT2*0.25*(input_tensor[0][2][0][1]+input_tensor[0][2][1][0]+input_tensor[2][0][0][1]+input_tensor[2][0][1][0]+input_tensor[0][1][0][2]+input_tensor[0][1][2][0]+input_tensor[1][0][0][2]+input_tensor[1][0][2][0]), // 18 // 19
          numbers::SQRT2*0.25*(input_tensor[1][2][0][1]+input_tensor[1][2][1][0]+input_tensor[2][1][0][1]+input_tensor[2][1][1][0]+input_tensor[0][1][1][2]+input_tensor[0][1][2][1]+input_tensor[1][0][1][2]+input_tensor[1][0][2][1]), // 19 // 20
          numbers::SQRT2*0.25*(input_tensor[1][2][0][2]+input_tensor[1][2][2][0]+input_tensor[2][1][0][2]+input_tensor[2][1][2][0]+input_tensor[0][2][1][2]+input_tensor[0][2][2][1]+input_tensor[2][0][1][2]+input_tensor[2][0][2][1])  // 20 // 21
        });

      }


      template <>
      const Tensor<3,3> &levi_civita<3>()
      {
        static const Tensor<3,3> t =
          []()
        {
          Tensor<3,3> permutation_operator_3d;

          permutation_operator_3d[0][1][2]  = 1;
          permutation_operator_3d[1][2][0]  = 1;
          permutation_operator_3d[2][0][1]  = 1;
          permutation_operator_3d[0][2][1]  = -1;
          permutation_operator_3d[1][0][2]  = -1;
          permutation_operator_3d[2][1][0]  = -1;
          return permutation_operator_3d;
        }();

        return t;
      }
    }


// Explicit instantiations

#define INSTANTIATE(dim) \
  template \
  IndexSet extract_locally_active_dofs_with_component(const DoFHandler<dim> &, \
                                                      const ComponentMask &); \
  \
  template \
  std::vector<Point<dim>> \
  get_unit_support_points(const SimulatorAccess<dim> &simulator_access); \
  \
  template \
  std::vector<std::string> \
  expand_dimensional_variable_names<dim> (const std::vector<std::string> &var_declarations); \
  \
  template \
  class VectorFunctionFromVelocityFunctionObject<dim>; \
  \
  template \
  Point<dim> Coordinates::spherical_to_cartesian_coordinates<dim>(const std::array<double,dim> &scoord); \
  \
  template \
  std::array<double,dim> Coordinates::cartesian_to_spherical_coordinates<dim>(const Point<dim> &position); \
  \
  template \
  Tensor<1,dim> Coordinates::spherical_to_cartesian_vector<dim>(const Tensor<1,dim> &spherical_vector, \
                                                                const Point<dim> &position); \
  \
  template \
  bool polygon_contains_point<dim>(const std::vector<Point<2>> &pointList, \
                                   const dealii::Point<2> &point); \
  \
  template \
  bool point_is_in_triangulation<dim>(const Mapping<dim> &mapping, \
                                      const parallel::distributed::Triangulation<dim> &triangulation, \
                                      const Point<dim> &point, \
                                      const MPI_Comm mpi_communicator); \
  \
  template \
  double signed_distance_to_polygon<dim>(const std::vector<Point<2>> &pointList, \
                                         const dealii::Point<2> &point); \
  \
  template \
  std::array<Tensor<1,dim>,dim-1> orthogonal_vectors (const Tensor<1,dim> &v); \
  \
  template \
  dealii::SymmetricTensor<2, dim, double> \
  derivative_of_weighted_p_norm_average (const double averaged_parameter, \
                                         const std::vector<double> &weights, \
                                         const std::vector<double> &values, \
                                         const std::vector<dealii::SymmetricTensor<2, dim, double>> &derivatives, \
                                         const double p); \
  \
  template \
  double compute_spd_factor(const double eta, \
                            const SymmetricTensor<2,dim> &strain_rate, \
                            const SymmetricTensor<2,dim> &dviscosities_dstrain_rate, \
                            const double safety_factor); \
  \
  template \
  Point<dim> convert_array_to_point<dim>(const std::array<double,dim> &array); \
  \
  template \
  std::array<double,dim> convert_point_to_array<dim>(const Point<dim> &point); \
  \
  template \
  SymmetricTensor<2,dim> nth_basis_for_symmetric_tensors (const unsigned int k); \
  \
  template \
  class NaturalCoordinate<dim>; \
  \
  template \
  void \
  project_cellwise(const Mapping<dim> &mapping, \
                   const DoFHandler<dim> &dof_handler, \
                   const unsigned int component_index, \
                   const Quadrature<dim> &quadrature, \
                   const std::function<void( \
                                             const DoFHandler<dim>::active_cell_iterator &, \
                                             const std::vector<Point<dim>> &, \
                                             std::vector<double> &)> &function, \
                   dealii::LinearAlgebra::distributed::Vector<double> &vec_result); \
  \
  template \
  void \
  project_cellwise(const Mapping<dim> &mapping, \
                   const DoFHandler<dim> &dof_handler, \
                   const unsigned int component_index, \
                   const Quadrature<dim> &quadrature, \
                   const std::function<void( \
                                             const DoFHandler<dim>::active_cell_iterator &, \
                                             const std::vector<Point<dim>> &, \
                                             std::vector<double> &)> &function, \
                   LinearAlgebra::BlockVector &vec_result);

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE

    // only instantiate for dim=3:
    template                \
    std::array<double,3> Coordinates::WGS84_coordinates<3>(const Point<3> &position);


    template double
    derivative_of_weighted_p_norm_average (const double averaged_parameter,
                                           const std::vector<double> &weights,
                                           const std::vector<double> &values,
                                           const std::vector<double> &derivatives,
                                           const double p);

    template Table<2,double> parse_input_table(const std::string &input_string,
                                               const unsigned int n_rows,
                                               const unsigned int n_columns,
                                               const std::string &property_name);
  }
}
