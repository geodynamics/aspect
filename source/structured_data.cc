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
#include <aspect/structured_data.h>

#include <aspect/geometry_model/box.h>
#include <aspect/geometry_model/two_merged_boxes.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/geometry_model/sphere.h>
#include <aspect/geometry_model/chunk.h>
#include <aspect/geometry_model/initial_topography_model/ascii_data.h>

#include <boost/lexical_cast.hpp>

namespace aspect
{
  namespace Utilities
  {

    template <int dim>
    StructuredDataLookup<dim>::StructuredDataLookup(const unsigned int components,
                                                    const double scale_factor)
      :
      components(components),
      data(components),
      maximum_component_value(components),
      scale_factor(scale_factor)
    {}



    template <int dim>
    StructuredDataLookup<dim>::StructuredDataLookup(const double scale_factor)
      :
      components(numbers::invalid_unsigned_int),
      data(),
      maximum_component_value(),
      scale_factor(scale_factor)
    {}



    template <int dim>
    std::vector<std::string>
    StructuredDataLookup<dim>::get_column_names() const
    {
      return data_component_names;
    }



    template <int dim>
    unsigned int
    StructuredDataLookup<dim>::get_column_index_from_name(const std::string &column_name) const
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
    StructuredDataLookup<dim>::get_column_name_from_index(const unsigned int column_index) const
    {
      AssertThrow(data_component_names.size() > column_index,
                  ExcMessage("There is no data column number " + Utilities::to_string(column_index)
                             + " in the current data file. The data file only contains "
                             + Utilities::to_string(data_component_names.size()) + " named columns."));

      return data_component_names[column_index];
    }



    template <int dim>
    double
    StructuredDataLookup<dim>::get_maximum_component_value(const unsigned int component) const
    {
      return maximum_component_value[component];
    }



    template <int dim>
    void
    StructuredDataLookup<dim>::reinit(const std::vector<std::string> &column_names,
                                      const std::vector<std::vector<double>> &coordinate_values_,
                                      const std::vector<Table<dim,double> > &raw_data)
    {
      Assert(coordinate_values_.size()==dim, ExcMessage("Invalid size of coordinate_values."));

      std::array<std::vector<double>,dim> coordinate_values;
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

      // min and max of the coordinates in the data file.
      std::array<std::pair<double,double>,dim> grid_extent;


      bool coordinate_values_are_equidistant = true;
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
    StructuredDataLookup<dim>::load_file(const std::string &filename,
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
    StructuredDataLookup<dim>::get_data(const Point<dim> &position,
                                        const unsigned int component) const
    {
      Assert(component<components, ExcMessage("Invalid component index"));
      return data[component]->value(position);
    }



    template <int dim>
    Tensor<1,dim>
    StructuredDataLookup<dim>::get_gradients(const Point<dim> &position,
                                             const unsigned int component)
    {
      return data[component]->gradient(position,0);
    }



    template <int dim>
    TableIndices<dim>
    StructuredDataLookup<dim>::compute_table_indices(const TableIndices<dim> &sizes, const unsigned int i) const
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
                                        std_cxx14::make_unique<Utilities::StructuredDataLookup<dim-1>>
                                        (components,
                                         this->scale_factor)));

          old_lookups.insert(std::make_pair(boundary_id,
                                            std_cxx14::make_unique<Utilities::StructuredDataLookup<dim-1>>
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

          lookups.push_back(std_cxx14::make_unique<Utilities::StructuredDataLookup<dim-1>> (components,
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

      lookup = std_cxx14::make_unique<Utilities::StructuredDataLookup<dim>> (components,
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
      lookup = std_cxx14::make_unique<Utilities::StructuredDataLookup<1>> (this->scale_factor);

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


// Explicit instantiations

    template class StructuredDataLookup<1>;
    template class StructuredDataLookup<2>;
    template class StructuredDataLookup<3>;
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
  }
}
