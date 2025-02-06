/*
  Copyright (C) 2011 - 2023 by the authors of the ASPECT code.

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
#include <aspect/geometry_model/two_merged_chunks.h>
#include <aspect/geometry_model/ellipsoidal_chunk.h>
#include <aspect/geometry_model/initial_topography_model/ascii_data.h>
#include <aspect/geometry_model/two_merged_chunks.h>

#include <deal.II/base/exceptions.h>

#include <boost/lexical_cast.hpp>
#include <regex>

#ifdef ASPECT_WITH_NETCDF

#include <netcdf.h>
#define AssertThrowNetCDF(error_code) \
  AssertThrow(error_code == NC_NOERR, dealii::ExcMessage("A NetCDF Error with code " + std::to_string(error_code) + " occurred."))

#else

#define AssertThrowNetCDF(error_code) \
  AssertThrow(false, ExcInternalError())

#endif



namespace aspect
{
  namespace Utilities
  {

    template <int dim>
    StructuredDataLookup<dim>::StructuredDataLookup(const unsigned int n_components,
                                                    const double scale_factor)
      :
      n_components(n_components),
      data(n_components),
      maximum_component_value(n_components),
      scale_factor(scale_factor),
      coordinate_values_are_equidistant(false)
    {}



    template <int dim>
    StructuredDataLookup<dim>::StructuredDataLookup(const double scale_factor)
      :
      n_components(numbers::invalid_unsigned_int),
      data(),
      maximum_component_value(),
      scale_factor(scale_factor),
      coordinate_values_are_equidistant(false)
    {}



    template <int dim>
    std::vector<std::string>
    StructuredDataLookup<dim>::get_column_names() const
    {
      return data_component_names;
    }



    template <int dim>
    bool
    StructuredDataLookup<dim>::has_equidistant_coordinates() const
    {
      return coordinate_values_are_equidistant;
    }



    template <int dim>
    const std::vector<double> &
    StructuredDataLookup<dim>::get_interpolation_point_coordinates(const unsigned int dimension) const
    {
      AssertThrow(dimension < dim,
                  ExcMessage("There is no spatial dimension number " + std::to_string(dimension)
                             + " in the current data file."));

      return coordinate_values[dimension];
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
    unsigned
    StructuredDataLookup<dim>::get_number_of_coordinates(const unsigned int dimension) const
    {
      return table_points[dimension];
    }



    namespace
    {
      /**
       * Return whether a set of coordinate points in dim space dimensions
       * are all equidistantly spaced within some small tolerance.
       */
      template <int dim>
      bool data_is_equidistant (const std::array<std::vector<double>,dim> &coordinate_values)
      {
        bool coordinate_values_are_equidistant = true;

        for (unsigned int d=0; d<dim; ++d)
          {
            const double grid_spacing = coordinate_values[d][1] - coordinate_values[d][0];

            for (unsigned int n = 1; n < coordinate_values[d].size(); ++n)
              {
                const double current_grid_spacing = coordinate_values[d][n] - coordinate_values[d][n-1];

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

        return coordinate_values_are_equidistant;
      }
    }


    template <int dim>
    void
    StructuredDataLookup<dim>::reinit(const std::vector<std::string> &column_names,
                                      std::vector<std::vector<double>> &&coordinate_values_,
                                      std::vector<Table<dim,double>> &&data_table,
                                      const MPI_Comm mpi_communicator,
                                      const unsigned int root_process)
    {
      // If this is the root process, or if the user did not request
      // sharing, set up the various member variables we need to compute
      // from the input data:
      if ((root_process == numbers::invalid_unsigned_int)
          ||
          ((root_process != numbers::invalid_unsigned_int)
           &&
           (Utilities::MPI::this_mpi_process(mpi_communicator) == root_process)))
        {
          Assert(coordinate_values_.size()==dim, ExcMessage("Invalid size of coordinate_values."));
          for (unsigned int d=0; d<dim; ++d)
            {
              this->coordinate_values[d] = std::move(coordinate_values_[d]);
              AssertThrow(this->coordinate_values[d].size()>1,
                          ExcMessage("Error: At least 2 entries per coordinate direction are required."));
              table_points[d] = this->coordinate_values[d].size();
            }

          n_components = column_names.size();
          data_component_names = column_names;
          Assert(data_table.size() == n_components,
                 ExcMessage("Error: Incorrect number of columns specified."));
          for (unsigned int c=0; c<n_components; ++c)
            Assert(data_table[c].size() == table_points,
                   ExcMessage("Error: One of the data tables has an incorrect size."));

          // compute maximum_component_value for each component:
          maximum_component_value = std::vector<double>(n_components,-std::numeric_limits<double>::max());
          for (unsigned int c=0; c<n_components; ++c)
            {
              const std::size_t n_elements = data_table[c].n_elements();
              for (std::size_t idx=0; idx<n_elements; ++idx)
                maximum_component_value[c] = std::max(maximum_component_value[c], data_table[c](
                                                        compute_table_indices(table_points, idx)));
            }

          // In case the data is specified on a grid that is equidistant
          // in each coordinate direction, we only need to store
          // (besides the data) the number of intervals in each direction and
          // the begin- and end-points of the coordinates.
          // In case the grid is not equidistant, we need to keep
          // all the coordinates in each direction, which is more costly.
          coordinate_values_are_equidistant = data_is_equidistant<dim> (coordinate_values);
        }

      // If the caller of this function requested it, then we have
      // set up member variables on the root process, but not on any of
      // the other processes. Broadcast the data to the remaining
      // processes
      if (root_process != numbers::invalid_unsigned_int)
        {
          coordinate_values                 = Utilities::MPI::broadcast (mpi_communicator,
                                                                         coordinate_values,
                                                                         root_process);
          n_components                        = Utilities::MPI::broadcast (mpi_communicator,
                                                                           n_components,
                                                                           root_process);
          data_component_names              = Utilities::MPI::broadcast (mpi_communicator,
                                                                         data_component_names,
                                                                         root_process);
          maximum_component_value           = Utilities::MPI::broadcast (mpi_communicator,
                                                                         maximum_component_value,
                                                                         root_process);
          coordinate_values_are_equidistant = Utilities::MPI::broadcast (mpi_communicator,
                                                                         coordinate_values_are_equidistant,
                                                                         root_process);
          table_points                      = Utilities::MPI::broadcast (mpi_communicator,
                                                                         table_points,
                                                                         root_process);

          // We can then also prepare the data tables for sharing between
          // processes
          for (unsigned int c = 0; c < n_components; ++c)
            data_table[c].replicate_across_communicator (mpi_communicator,
                                                         root_process);
        }

      Assert(data_table.size() == n_components,
             ExcMessage("Error: Incorrect number of columns specified."));
      for (unsigned int c=0; c<n_components; ++c)
        Assert(data_table[c].size() == table_points,
               ExcMessage("Error: One of the data tables has an incorrect size."));


      // For each data component, set up a GridData,
      // its type depending on the read-in grid.
      data.resize(n_components);
      for (unsigned int c = 0; c < n_components; ++c)
        {
          if (coordinate_values_are_equidistant)
            {
              std::array<unsigned int,dim> table_intervals;
              for (unsigned int d=0; d<dim; ++d)
                table_intervals[d] = table_points[d]-1;

              // The min and max of the coordinates in the data file.
              std::array<std::pair<double,double>,dim> grid_extent;
              for (unsigned int d=0; d<dim; ++d)
                {
                  grid_extent[d].first = coordinate_values[d][0];
                  grid_extent[d].second = coordinate_values[d][table_points[d]-1];

                  Assert(table_intervals[d] >= 1,
                         ExcMessage("There needs to be at least one subinterval in each "
                                    "coordinate direction."));
                  Assert(grid_extent[d].first <
                         grid_extent[d].second,
                         ExcMessage("The interval in each coordinate direction needs "
                                    "to have positive size"));
                }

              data[c]
                = std::make_unique<Functions::InterpolatedUniformGridData<dim>> (std::move(grid_extent),
                                                                                  std::move(table_intervals),
                                                                                  std::move(data_table[c]));
            }
          else
            // Create the object and move the big objects. Due to an old design flaw,
            // the current class stores a copy of the 'coordinate_values' and some
            // plugins actually use it too, i.e., we can't move the data out of
            // this object. In another design flaw, the deal.II classes
            // do not make the coordinate values accessible, so we have to continue
            // storing a copy. In other words, we can't just move stuff --
            // we have to make a copy of 'coordinate_values' and then move
            // that copy.
            //
            // (The call to std::move on the first argument is unnecessary: We
            // create a temporary object, and that's an rvalue that the constructor
            // we call would bind to. But never a bad idea to be explicit.)
            data[c]
              = std::make_unique<Functions::InterpolatedTensorProductGridData<dim>>
                (std::move(std::array<std::vector<double>,dim>(this->coordinate_values)),
                 std::move(data_table[c]));
        }
    }



    template <int dim>
    void
    StructuredDataLookup<dim>::load_ascii(const std::string &filename,
                                          const MPI_Comm comm)
    {
      const unsigned int root_process = 0;

      std::vector<std::string> column_names;
      std::vector<Table<dim,double>> data_tables;
      std::vector<std::vector<double>> coordinate_values(dim);

      // If this is the root process, set up the various member variables we need to compute
      // from the input data
      if (Utilities::MPI::this_mpi_process(comm) == root_process)
        {
          // Grab the values already stored in this class (if they exist), this way we can
          // check if somebody changes the size of the table over time and error out (see below)
          TableIndices<dim> new_table_points = this->table_points;

          // We do not need to distribute the contents as we are using shared data
          // to place it later. Therefore, just pass MPI_COMM_SELF (i.e.,
          // a communicator with just a single MPI process) and no distribution
          // will happen.
          std::stringstream in(read_and_distribute_file_content(filename,
                                                                MPI_COMM_SELF));

          // Read header lines and table size
          while (in.peek() == '#')
            {
              std::string line;
              std::getline(in,line);
              std::stringstream linestream(line);
              std::string word;
              while (linestream >> word)
                if (word == "POINTS:")
                  for (unsigned int i = 0; i < dim; ++i)
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

          for (unsigned int i = 0; i < dim; ++i)
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
                  // and have read the first data field. Save number of n_components, and
                  // make sure there is no contradiction if the n_components were already given to
                  // the constructor of this class.
                  if (n_components == numbers::invalid_unsigned_int)
                    n_components = name_column_index - dim;
                  else if (name_column_index != 0)
                    AssertThrow (n_components+dim == name_column_index,
                                 ExcMessage("The number of expected data columns and the "
                                            "list of column names at the beginning of the data file "
                                            + filename + " do not match. The file should contain "
                                            + Utilities::int_to_string(name_column_index) + " column "
                                            "names (one for each dimension and one per data column), "
                                            "but it only has " + Utilities::int_to_string(n_components+dim) +
                                            " column names."));
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
          AssertThrow (n_components != numbers::invalid_unsigned_int,
                       ExcMessage("ERROR: number of n_components in " + filename + " could not be "
                                  "determined automatically. Either add a header with column "
                                  "names or pass the number of columns in the StructuredData "
                                  "constructor."));
          data_tables.resize(n_components, data_table);

          for (unsigned int d=0; d<dim; ++d)
            coordinate_values[d].resize(new_table_points[d]);

          if (column_names.size()==0)
            {
              // set default column names:
              for (unsigned int c=0; c<n_components; ++c)
                column_names.push_back("column " + Utilities::int_to_string(c,2));
            }

          // Make sure the data file actually has as many columns as we think it has
          // (either based on the header, or based on what was passed to the constructor).
          const std::streampos position = in.tellg();
          std::string first_data_row;
          std::getline(in, first_data_row);
          std::stringstream linestream(first_data_row);
          std::string column_entry;

          // We have already read in the first data entry above in the try/catch block,
          // so there's one more column in the file than in the line we just read in.
          unsigned int number_of_entries = 1;
          while (linestream >> column_entry)
            number_of_entries += 1;

          AssertThrow ((number_of_entries) == column_names.size()+dim,
                       ExcMessage("ERROR: The number of columns in the data file " + filename +
                                  " is incorrect. It needs to have " + Utilities::int_to_string(column_names.size()+dim) +
                                  " columns, but the first row has " + Utilities::int_to_string(number_of_entries) +
                                  " columns."));

          // Go back to the position in the file where we started the check for the column numbers.
          in.seekg (position);

          // Finally read data lines:
          std::size_t read_data_entries = 0;
          do
            {
              // what row and column of the file are we in?
              const std::size_t column_num = read_data_entries%(n_components+dim);
              const std::size_t row_num = read_data_entries/(n_components+dim);
              const TableIndices<dim> idx = compute_table_indices(new_table_points, row_num);

              if (column_num < dim)
                {
                  // This is a coordinate. Store (and check that they are consistent)
                  const double old_value = coordinate_values[column_num][idx[column_num]];

                  AssertThrow(old_value == 0. ||
                              (std::abs(old_value-temp_data) < 1e-8*std::abs(old_value)),
                              ExcMessage("Invalid coordinate in column "
                                         + Utilities::int_to_string(column_num) + " in row "
                                         + Utilities::int_to_string(row_num)
                                         + " in file " + filename +
                                         "\nThis class expects the coordinates to be structured, meaning "
                                         "the coordinate values in each coordinate direction repeat exactly "
                                         "each time. This also means each row in the data file has to have "
                                         "the same number of columns as the first row containing data."));

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

          const std::size_t n_expected_data_entries = (n_components + dim) * data_table.n_elements();
          AssertThrow(read_data_entries == n_expected_data_entries,
                      ExcMessage ("While reading the data file '" + filename + "' the ascii data "
                                  "plugin has reached the end of the file, but has not found the "
                                  "expected number of data values considering the spatial dimension, "
                                  "data columns, and number of lines prescribed by the POINTS header "
                                  "of the file. Please check the number of data "
                                  "lines against the POINTS header in the file."));
        }

      // deal.II supports sharing data (since 9.4), so we have to
      // set up member variables on the root process, but not on any of
      // the other processes. So broadcast the data to the remaining
      // processes -- the code above really only wrote into one
      // member variable ('n_components'), so that is the only one we
      // have to broadcast.
      //
      // The first three arguments to the call to reinit() below will
      // only be read on the root process, and so it is totally ok
      // that we are passing empty tables on all other processes. In
      // the case of 'data_table', we do have to make sure that it is
      // an array of the right size, though, even though the array
      // contains only empty tables.
      {
        n_components = Utilities::MPI::broadcast (comm,
                                                  n_components,
                                                  root_process);
        coordinate_values = Utilities::MPI::broadcast (comm,
                                                       coordinate_values,
                                                       root_process);
        column_names = Utilities::MPI::broadcast (comm,
                                                  column_names,
                                                  root_process);

        if (Utilities::MPI::this_mpi_process(comm) != root_process)
          data_tables.resize (n_components);
      }

      // Finally create the data. We want to call the move-version of reinit() so
      // that the data doesn't have to be copied, so use std::move on all big
      // objects.
      this->reinit(column_names,
                   std::move(coordinate_values),
                   std::move(data_tables),
                   comm,
                   root_process);
    }



    template <int dim>
    void
    StructuredDataLookup<dim>::load_netcdf(const std::string &filename, const std::vector<std::string> &data_column_names_)
    {
#ifndef ASPECT_WITH_NETCDF
      (void)filename;
      (void)data_column_names_;
      AssertThrow(false, ExcMessage("Loading NetCDF files is only supported if ASPECT is configured with the NetCDF library!"));
#else
      TableIndices<dim> new_table_points;
      std::vector<std::string> coordinate_column_names(dim);
      std::vector<Table<dim,double>> data_tables;
      std::vector<std::string> data_column_names = data_column_names_;
      std::vector<std::vector<double>> coordinate_values(dim);

      int ncid;
      int status;

      status = nc_open(filename.c_str(), NC_NOWRITE, &ncid);
      AssertThrowNetCDF(status);

      int ndims, nvars, ngatts, unlimdimid;
      status = nc_inq(ncid, &ndims, &nvars, &ngatts, &unlimdimid);
      AssertThrowNetCDF(status);

      // The number of dimensions ndims in the netcdf file can be
      // different than dim. In fact, each variable (data column in
      // our notation) has a subset of those dimensions associated
      // with it. That also means that variables can have a different
      // number of dimensions and/or a different subset. We can only
      // use variables with dim dimensions (our template argument) and
      // only those that all use the same dimids inside the netcdf
      // file.
      int dimids_to_use[dim] = {}; // dimids of the coordinate columns to use
      std::vector<int> varids_to_use; // all netCDF varids of the data columns

      if (data_column_names.empty())
        {
          // The user did not ask for a specific list of data
          // columns. Let's find all columns we can possible load. We
          // find the first data column with the correct number of
          // dims. This one will determine the dimids of the
          // coordinates to use. Following that, we pick all other
          // data columns with the same coordinates.

          for (int varid=0; varid<nvars; ++varid)
            {
              // Each netCDF dimension also has an associated variable that stores the
              // coordinate data. We are looking for data columns, so skip them:
              if (varid < ndims)
                continue;

              char  var_name[NC_MAX_NAME];
              nc_type xtype;
              int var_ndims;
              int var_dimids[NC_MAX_VAR_DIMS];
              int var_natts;

              status = nc_inq_var (ncid, varid, var_name, &xtype, &var_ndims, var_dimids,
                                   &var_natts);
              AssertThrowNetCDF(status);

              // only consider data that has dim variables:
              if (var_ndims == dim)
                {
                  bool use = true;
                  if (varids_to_use.size()>0)
                    {
                      // This is not the first data column, so we can only use it if
                      // it uses the same dim variables as the first data column we
                      // found.
                      for (int i=0; i<dim; ++i)
                        if (dimids_to_use[i]!=var_dimids[i])
                          {
                            use=false;
                            break;
                          }
                    }
                  else
                    {
                      // This is the first data column we found, so grab the ids of
                      // the dimensions to use and store in dimids_to_use:
                      for (int i=0; i<dim; ++i)
                        {
                          size_t length;
                          status = nc_inq_dim(ncid, var_dimids[i], nullptr, &length);
                          dimids_to_use[i] = var_dimids[i];
                          // dimensions are specified in reverse order in the nc file:
                          new_table_points[dim-1-i] = length;
                        }
                    }

                  if (use)
                    {
                      varids_to_use.push_back(varid);
                      data_column_names.push_back(var_name);
                    }
                }
            }

        }
      else
        {
          // The user wants a specific list of columns, so lets find them.

          for (const auto &cur_name: data_column_names)
            {
              bool found = false;

              for (int varid=0; varid<nvars; ++varid)
                {
                  // Each netCDF dimension also has an associated variable that stores the
                  // coordinate data. We are looking for data columns, so skip them:
                  if (varid < ndims)
                    continue;

                  char  var_name[NC_MAX_NAME];
                  nc_type xtype;
                  int var_ndims;
                  int var_dimids[NC_MAX_VAR_DIMS];
                  int var_natts;

                  status = nc_inq_var (ncid, varid, var_name, &xtype, &var_ndims, var_dimids,
                                       &var_natts);
                  AssertThrowNetCDF(status);
                  if (cur_name != var_name)
                    continue;

                  found = true;

                  if (var_ndims == dim)
                    {
                      bool use = true;
                      if (varids_to_use.size()>0)
                        {
                          for (int i=0; i<dim; ++i)
                            if (dimids_to_use[i]!=var_dimids[i])
                              {
                                use=false;
                                break;
                              }
                        }
                      else
                        {
                          for (int i=0; i<dim; ++i)
                            {
                              size_t length;
                              status = nc_inq_dim(ncid, var_dimids[i], nullptr, &length);
                              dimids_to_use[i] = var_dimids[i];
                              // dimensions are specified in reverse order in the nc file:
                              new_table_points[dim-1-i] = length;
                            }
                        }

                      AssertThrow(use, ExcMessage(
                                    "You asked to include column '" + cur_name + "', but it unfortunately has different dimensions than the first column you chose."
                                  ));


                      varids_to_use.push_back(varid);
                    }
                  else
                    AssertThrow(false, ExcMessage(
                                  "You asked to include column '" + cur_name + "', but it unfortunately has an incorrect number of dimensions."
                                ));

                }
              AssertThrow(found, ExcMessage(
                            "You asked to include column '" + cur_name + "', but it was not found!"
                          ));

            }
        }


      // Extract names of the coordinates
      for (int idx=0; idx<dim; ++idx)
        {
          char  name[NC_MAX_NAME];
          size_t length;
          status = nc_inq_dim(ncid, dimids_to_use[idx], name, &length);
          AssertThrowNetCDF(status);
          coordinate_column_names[idx] = name;
        }

      {
        // Now load coordinate data

        for (int d=0; d<dim; ++d)
          {
            // dimensions are specified in reverse order in the nc file:
            int varid = dimids_to_use[dim-1-d];

            nc_type xtype;
            int ndims;
            int dimids[NC_MAX_VAR_DIMS];
            char  name[NC_MAX_NAME];
            int natts;

            status = nc_inq_var (ncid, varid, name, &xtype, &ndims, dimids,
                                 &natts);
            AssertThrow(ndims == 1, ExcMessage("A variable of a dimension should have only one dimension."));

            AssertThrow(xtype == NC_DOUBLE || xtype == NC_FLOAT, ExcMessage("We only support float or double data."));

            coordinate_values[d].resize(new_table_points[d]);
            status = nc_get_var_double(ncid, varid, coordinate_values[d].data());
            AssertThrowNetCDF(status);
          }
      }

      {
        // Finally load the data for each column
        data_tables.resize(varids_to_use.size());
        std::vector<double> raw_data;

        for (unsigned int var = 0; var<varids_to_use.size(); ++var)
          {
            // Allocate space
            data_tables[var].TableBase<dim,double>::reinit(new_table_points);
            const std::size_t n_elements = data_tables[var].n_elements();
            raw_data.resize(n_elements);

            // Load the data
            status = nc_get_var_double(ncid, varids_to_use[var], raw_data.data());
            AssertThrowNetCDF(status);

            // .. and copy it over:
            for (std::size_t n = 0; n < n_elements; ++n)
              {
                TableIndices<dim> ind = compute_table_indices(new_table_points, n);
                TableIndices<dim> ind_to_use;
                for (int i=0; i<dim; ++i)
                  ind_to_use[i] = ind[dimids_to_use[i]];

                data_tables[var](ind) = scale_factor * raw_data[n];
              }
          }

      }


      status = nc_close(ncid);
      AssertThrowNetCDF(status);

      // ready to go:
      this->reinit(data_column_names, std::move(coordinate_values),std::move(data_tables));
#endif
    }

    template <int dim>
    void
    StructuredDataLookup<dim>::load_file(const std::string &filename,
                                         const MPI_Comm communicator)
    {
      const bool is_netcdf_filename = std::regex_search(filename, std::regex("\\.(nc|NC)$"));
      if (is_netcdf_filename)
        load_netcdf(filename);
      else
        load_ascii(filename, communicator);
    }

    template <int dim>
    double
    StructuredDataLookup<dim>::get_data(const Point<dim> &position,
                                        const unsigned int component,
                                        const bool crash_if_not_in_range) const
    {
      Assert(component<n_components, ExcMessage("Invalid component index"));

      if (crash_if_not_in_range)
        {
          const std::vector<double> &x_coordinates = get_interpolation_point_coordinates(0);

          AssertThrow (position[0] >= x_coordinates[0] && position[0] <= x_coordinates[x_coordinates.size()-1],
                       ExcMessage("The requested position "
                                  + std::to_string(position[0])
                                  + " is outside the range of the data (minimum value = "
                                  + std::to_string(x_coordinates[0])
                                  + " , maximum value = "
                                  + std::to_string(x_coordinates[x_coordinates.size()-1])
                                  + ")."
                                 ));

          const std::vector<double> &y_coordinates = get_interpolation_point_coordinates(1);

          AssertThrow (position[1] >= y_coordinates[0] && position[1] <= y_coordinates[y_coordinates.size()-1],
                       ExcMessage("The requested position "
                                  + std::to_string(position[1])
                                  + " is outside the range of the data (minimum value = "
                                  + std::to_string(y_coordinates[0])
                                  + " , maximum value = "
                                  + std::to_string(y_coordinates[y_coordinates.size()-1])
                                  + ")."
                                 ));
        }

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
    StructuredDataLookup<dim>::compute_table_indices(const TableIndices<dim> &sizes, const std::size_t idx) const
    {
      TableIndices<dim> result;
      result[0] = idx % sizes[0];
      if (dim >= 2)
        result[1] = (idx / sizes[0]) % sizes[1];
      if (dim == 3)
        result[2] = idx / (sizes[0] * sizes[1]);

      return result;
    }



    template <int dim>
    AsciiDataBase<dim>::AsciiDataBase ()
      = default;


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

    namespace
    {
      template <int dim>
      void
      check_supported_geometry_models(const GeometryModel::Interface<dim> &geometry_model)
      {
        AssertThrow ((Plugins::plugin_type_matches<const GeometryModel::SphericalShell<dim>> (geometry_model))
                     || (Plugins::plugin_type_matches<const GeometryModel::Chunk<dim>> (geometry_model))
                     || (Plugins::plugin_type_matches<const GeometryModel::TwoMergedChunks<dim>> (geometry_model))
                     || (Plugins::plugin_type_matches<const GeometryModel::EllipsoidalChunk<dim>> (geometry_model))
                     || (Plugins::plugin_type_matches<const GeometryModel::Sphere<dim>> (geometry_model))
                     || (Plugins::plugin_type_matches<const GeometryModel::Box<dim>> (geometry_model))
                     || (Plugins::plugin_type_matches<const GeometryModel::TwoMergedBoxes<dim>> (geometry_model)),
                     ExcMessage ("This ascii data plugin can only be used with supported "
                                 "geometry models."));
      }
    }



    template <int dim>
    AsciiDataBoundary<dim>::AsciiDataBoundary ()
      :
      current_file_number(numbers::invalid_unsigned_int),
      first_data_file_number(numbers::invalid_unsigned_int),
      decreasing_file_order(false),
      data_file_time_step(numbers::signaling_nan<double>()),
      time_weight(numbers::signaling_nan<double>()),
      time_dependent(false),
      lookups(),
      old_lookups()
    {}



    template <int dim>
    void
    AsciiDataBoundary<dim>::initialize(const std::set<types::boundary_id> &boundary_ids,
                                       const unsigned int n_components)
    {
      check_supported_geometry_models(this->get_geometry_model());

      for (const auto &boundary_id : boundary_ids)
        {
          lookups.insert(std::make_pair(boundary_id,
                                        std::make_unique<Utilities::StructuredDataLookup<dim-1>>
                                        (n_components,
                                         this->scale_factor)));

          // Set the first file number and load the first files
          current_file_number = first_data_file_number;

          const std::string filename (create_filename (current_file_number, boundary_id));

          this->get_pcout() << std::endl << "   Loading Ascii data boundary file "
                            << filename << '.' << std::endl << std::endl;


          AssertThrow(Utilities::fexists(filename, this->get_mpi_communicator()) || filename_is_url(filename),
                      ExcMessage (std::string("Ascii data file <")
                                  +
                                  filename
                                  +
                                  "> not found!"));
          lookups.find(boundary_id)->second->load_file(filename,this->get_mpi_communicator());

          if (time_dependent == true)
            {
              old_lookups.insert(std::make_pair(boundary_id,
                                                std::make_unique<Utilities::StructuredDataLookup<dim-1>>
                                                (n_components,
                                                 this->scale_factor)));

              const int next_file_number =
                (decreasing_file_order) ?
                current_file_number - 1
                :
                current_file_number + 1;

              const std::string filename (create_filename (next_file_number, boundary_id));
              if (Utilities::fexists(filename, this->get_mpi_communicator()))
                {
                  this->get_pcout() << std::endl << "   Also loading next Ascii data boundary file "
                                    << filename << '.' << std::endl << std::endl;
                  lookups.find(boundary_id)->second.swap(old_lookups.find(boundary_id)->second);
                  lookups.find(boundary_id)->second->load_file(filename, this->get_mpi_communicator());
                }
              else
                {
                  // next file not found, issue warning and end looking for new files
                  end_time_dependence ();
                }
            }
        }
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
       * by @p filenumber. Options: (1) specify both boundary_name and
       * filenumber placeholders (%s and %d), (2) do not specify any
       * boundary_name and filenumber placeholders, (3) only specify
       * boundary_name placeholder (%s). Do not only specify filenumber
       * placeholder (%d). Placeholders order is first boundary_name
       * (%s), then filenumber (%d).
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
      if (fexists(result, this->get_mpi_communicator()))
        return result;

      // Backwards compatibility check: people might still be using the old
      // names of the top/bottom boundary. If they do, print a warning but
      // accept those files.
      std::string compatible_result;
      if (boundary_name == "top")
        {
          compatible_result = replace_placeholders(templ, "surface", filenumber);
          if (!fexists(compatible_result, this->get_mpi_communicator()))
            compatible_result = replace_placeholders(templ, "outer", filenumber);
        }
      else if (boundary_name == "bottom")
        compatible_result = replace_placeholders(templ, "inner", filenumber);

      if (!fexists(result, this->get_mpi_communicator()) && fexists(compatible_result, this->get_mpi_communicator()))
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
      if (time_dependent == true)
        {
          // always initialize with the start time during model setup, even
          // when restarting (to have identical setup in both cases)
          double model_time = this->get_parameters().start_time;

          // if we are past initialization use the current time instead
          if (this->simulator_is_past_initialization())
            model_time = this->get_time();

          const double time_steps_since_start = model_time / data_file_time_step;
          // whether we need to update our data files. This looks so complicated
          // because we need to catch increasing and decreasing file orders and all
          // possible first_data_file_numbers.
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
                            << filename << '.' << std::endl << std::endl;
          if (Utilities::fexists(filename, this->get_mpi_communicator()))
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
                        << filename << '.' << std::endl << std::endl;
      if (Utilities::fexists(filename, this->get_mpi_communicator()))
        {
          lookups.find(boundary_id)->second.swap(old_lookups.find(boundary_id)->second);
          lookups.find(boundary_id)->second->load_file(filename,this->get_mpi_communicator());
        }

      // If next file does not exist, end time dependent part with current_time_step and issue warning.
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



    namespace
    {
      /**
       * Determines which of the dimensions of a position are aligned with
       * a certain @p boundary_id. E.g. the left boundary of a box
       * model extents in the y and z direction (position[1] and
       * position[2]), therefore the function would return [1,2] for dim==3
       * or [1] for dim==2. We are lucky that these indices are identical
       * for all existing geometries (if we use natural coordinates),
       * therefore we do not need to distinguish between them. If we
       * introduce a geometry model for which boundaries are not aligned
       * with natural coordinates this needs to become a function in
       * the interface of the geometry model and needs to be specialized
       * for distinct geometry models.
       */
      template <int dim>
      std::array<unsigned int,dim-1>
      get_boundary_dimensions (const types::boundary_id boundary_id)
      {
        std::array<unsigned int,dim-1> boundary_dimensions;
        boundary_dimensions.fill(numbers::invalid_unsigned_int);

        switch (dim)
          {
            case 2:
              if ((boundary_id == 2) || (boundary_id == 3) || (boundary_id == 4) || (boundary_id == 5))
                {
                  boundary_dimensions[0] = 0;
                }
              else if ((boundary_id == 0) || (boundary_id == 1))
                {
                  boundary_dimensions[0] = 1;
                }
              else
                {
                  AssertThrow(false,ExcNotImplemented());
                }
              break;

            case 3:
              if ((boundary_id == 4) || (boundary_id == 5) || (boundary_id == 8) || (boundary_id == 9))
                {
                  boundary_dimensions[0] = 0;
                  boundary_dimensions[1] = 1;
                }
              else if ((boundary_id == 0) || (boundary_id == 1))
                {
                  boundary_dimensions[0] = 1;
                  boundary_dimensions[1] = 2;
                }
              else if ((boundary_id == 2) || (boundary_id == 3) || (boundary_id == 6) || (boundary_id == 7))
                {
                  boundary_dimensions[0] = 0;
                  boundary_dimensions[1] = 2;
                }
              else
                {
                  AssertThrow(false,ExcNotImplemented());
                }

              break;

            default:
              AssertThrow(false,ExcNotImplemented());
          }
        return boundary_dimensions;
      }



      /**
       * Convert a point @p position in cartesian coordinates into the set of coordinates
       * used in the structured data input files, i.e. coordinates in the
       * natural coordinate system of the given @p geometry_model.
       */
      template <int dim>
      Point<dim>
      data_coordinates_from_position (const Point<dim> &position,
                                      const aspect::GeometryModel::Interface<dim> &geometry_model)
      {
        const std::array<double,dim> natural_position = geometry_model.cartesian_to_natural_coordinates(position);
        Point<dim> data_coordinates = Utilities::convert_array_to_point<dim>(natural_position);

        // The chunk model has latitude as natural coordinate. We need to convert this to colatitude
        // to allow for consistent data files between chunk geometries and spherical geometries.
        if (dim == 3 &&
            (Plugins::plugin_type_matches<const GeometryModel::Chunk<dim>> (geometry_model) ||
             Plugins::plugin_type_matches<const GeometryModel::TwoMergedChunks<dim>> (geometry_model) ||
             Plugins::plugin_type_matches<const GeometryModel::EllipsoidalChunk<dim>> (geometry_model)))
          {
            data_coordinates[2] = numbers::PI/2. - data_coordinates[2];
          }

        return data_coordinates;
      }



      /**
       * Convert a point @p data_coordinates using coordinates appropriate for
       * volumetric structured data input files (e.g. produced by the
       * data_coordinates_from_position() function above) into appropriate
       * coordinates on the boundary given by @p boundary_indicator.
       * Because all existing geometry models use boundaries that are aligned
       * with natural coordinate directions this means finding the correct n_components
       * of the input coordinates and returning a point with one dimension less
       * than the input point that is determined by these coordinates.
       */
      template <int dim>
      Point<dim-1>
      boundary_coordinates_from_data_coordinates (const Point<dim> &data_coordinates,
                                                  const types::boundary_id boundary_indicator)
      {
        const std::array<unsigned int,dim-1> boundary_dimensions =
          get_boundary_dimensions<dim>(boundary_indicator);

        Point<dim-1> boundary_coordinates;
        for (unsigned int i = 0; i < dim-1; ++i)
          boundary_coordinates[i] = data_coordinates[boundary_dimensions[i]];

        return boundary_coordinates;
      }
    }



    template <int dim>
    double
    AsciiDataBoundary<dim>::
    get_data_component (const types::boundary_id             boundary_indicator,
                        const Point<dim>                    &position,
                        const unsigned int                   component) const
    {
      const Point<dim> data_coordinates = data_coordinates_from_position(position, this->get_geometry_model());
      const Point<dim-1> boundary_coordinates = boundary_coordinates_from_data_coordinates(data_coordinates, boundary_indicator);

      Assert (lookups.find(boundary_indicator) != lookups.end(),
              ExcInternalError());
      const double data = lookups.find(boundary_indicator)->second->get_data(boundary_coordinates,component);

      if (!time_dependent)
        return data;

      const double old_data = old_lookups.find(boundary_indicator)->second->get_data(boundary_coordinates,component);

      return time_weight * data + (1 - time_weight) * old_data;
    }


    template <int dim>
    Tensor<1,dim-1>
    AsciiDataBoundary<dim>::vector_gradient (const types::boundary_id             boundary_indicator,
                                             const Point<dim>                    &position,
                                             const unsigned int                   component) const
    {
      const Point<dim> data_coordinates = data_coordinates_from_position(position, this->get_geometry_model());
      const Point<dim-1> boundary_coordinates = boundary_coordinates_from_data_coordinates(data_coordinates, boundary_indicator);

      const Tensor<1,dim-1>  gradients = lookups.find(boundary_indicator)->second->get_gradients(boundary_coordinates,component);

      if (!time_dependent)
        return gradients;

      const Tensor<1,dim-1> old_gradients = old_lookups.find(boundary_indicator)->second->get_gradients(boundary_coordinates,component);

      return time_weight * gradients + (1 - time_weight) * old_gradients;
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
                                                const std::string &subsection_name,
                                                const bool declare_time_dependent_parameters)
    {
      Utilities::AsciiDataBase<dim>::declare_parameters(prm,
                                                        default_directory,
                                                        default_filename,
                                                        subsection_name);

      prm.enter_subsection (subsection_name);
      {
        if (declare_time_dependent_parameters == true)
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
                               "The `First data file model time' parameter "
                               "has been deactivated and will be removed in a future release. "
                               "Do not use this parameter and instead provide data files "
                               "starting from the model start time.");
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
        else
          {
            prm.declare_entry ("Data file name",
                               default_filename,
                               Patterns::Anything (),
                               "The file name of the model data. Provide file in format: "
                               "(File name).\\%s, where \\%s is a string specifying "
                               "the boundary of the model according to the names of the boundary "
                               "indicators (of the chosen geometry model).");
          }
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    AsciiDataBoundary<dim>::parse_parameters (ParameterHandler &prm,
                                              const std::string &subsection_name,
                                              const bool parse_time_dependent_parameters)
    {
      Utilities::AsciiDataBase<dim>::parse_parameters(prm,
                                                      subsection_name);

      prm.enter_subsection(subsection_name);
      {
        if (parse_time_dependent_parameters == true)
          {
            data_file_time_step             = prm.get_double ("Data file time step");
            const double first_data_file_model_time      = prm.get_double ("First data file model time");

            AssertThrow (first_data_file_model_time == 0.0,
                         ExcMessage("The `First data file model time' parameter "
                                    "has been deactivated and will be removed in a future release. "
                                    "Do not use this parameter and instead provide data files "
                                    "starting from the model start time."));

            first_data_file_number          = prm.get_integer("First data file number");
            decreasing_file_order           = prm.get_bool   ("Decreasing file order");

            if (this->convert_output_to_years() == true)
              {
                data_file_time_step        *= year_in_seconds;
              }
          }
        else
          {
            AssertThrow (create_filename (0, 0) == create_filename (1, 0),
                         ExcMessage("A boundary data file name was used that contained a placeholder for "
                                    "the current file number, but this AsciiDataBoundary object does not support "
                                    "time dependent information. Please remove the file number placeholder."));
          }

        // if filename does not contain a placeholder for timestep, no time dependence
        // do not issue a warning, the parameter file is specifying exactly one file.
        if (create_filename (0, 0) == create_filename (1, 0))
          time_dependent = false;
        else
          time_dependent = true;
      }
      prm.leave_subsection();
    }



    template <int dim>
    AsciiDataLayered<dim>::AsciiDataLayered ()
      = default;



    template <int dim>
    void
    AsciiDataLayered<dim>::initialize(const unsigned int n_components)
    {
      check_supported_geometry_models(this->get_geometry_model());

      // Create the lookups for each file
      number_of_layer_boundaries = data_file_names.size();
      for (unsigned int i=0; i<number_of_layer_boundaries; ++i)
        {
          const std::string filename = data_directory + data_file_names[i];
          AssertThrow(Utilities::fexists(filename, this->get_mpi_communicator()) || filename_is_url(filename),
                      ExcMessage (std::string("Ascii data file <")
                                  +
                                  filename
                                  +
                                  "> not found!"));

          lookups.push_back(std::make_unique<Utilities::StructuredDataLookup<dim-1>> (n_components,
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
      const Point<dim> data_coordinates = data_coordinates_from_position(position, this->get_geometry_model());

      double vertical_position;
      Point<dim-1> horizontal_position;
      if (this->get_geometry_model().natural_coordinate_system() == Utilities::Coordinates::CoordinateSystem::cartesian)
        {
          // in cartesian coordinates, the vertical component comes last
          vertical_position = data_coordinates[dim-1];
          for (unsigned int i = 0; i < dim-1; ++i)
            horizontal_position[i] = data_coordinates[i];
        }
      else
        {
          // in spherical coordinates, the vertical component comes first
          vertical_position = data_coordinates[0];
          for (unsigned int i = 0; i < dim-1; ++i)
            horizontal_position[i] = data_coordinates[i+1];
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
      : slice_data(false),
        rotation_matrix()
    {}



    template <int dim>
    void
    AsciiDataInitial<dim>::initialize (const unsigned int n_components)
    {
      const std::string filename = this->data_directory + this->data_file_name;

      this->get_pcout() << std::endl << "   Loading Ascii data initial file "
                        << filename << '.' << std::endl << std::endl;


      AssertThrow(Utilities::fexists(filename, this->get_mpi_communicator()) || filename_is_url(filename),
                  ExcMessage (std::string("Ascii data file <")
                              +
                              filename
                              +
                              "> not found!"));

      if (slice_data == true)
        {
          slice_lookup = std::make_unique<Utilities::StructuredDataLookup<3>> (n_components,
                                                                                this->scale_factor);
          slice_lookup->load_file(filename, this->get_mpi_communicator());
        }
      else
        {
          lookup = std::make_unique<Utilities::StructuredDataLookup<dim>> (n_components,
                                                                            this->scale_factor);
          lookup->load_file(filename, this->get_mpi_communicator());
        }
    }



    template <int dim>
    double
    AsciiDataInitial<dim>::
    get_data_component (const Point<dim> &position,
                        const unsigned int component) const
    {
      // Handle the special case of slicing through data first
      if (slice_data == true)
        {
          // This implies a SphericalShell, a 2d model, and a 3d dataset as asserted
          // in the parse_parameters function below.

          // Compute the coordinates of a 3d point based on the 2d position.
          const Tensor<1,3> position_tensor({position[0], position[1], 0.0});
          const Point<3> rotated_position (rotation_matrix * position_tensor);

          const std::array<double,3> spherical_position =
            Utilities::Coordinates::cartesian_to_spherical_coordinates(rotated_position);

          return slice_lookup->get_data(Point<3>(Tensor<1,3>(ArrayView<const double>(spherical_position))),component);
        }

      // else slice_data == false: model and dataset have the same dimension
      const Point<dim> data_coordinates = data_coordinates_from_position(position, this->get_geometry_model());

      return lookup->get_data(data_coordinates,component);
    }



    template <int dim>
    void
    AsciiDataInitial<dim>::declare_parameters (ParameterHandler  &prm,
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
        prm.declare_entry("Slice dataset in 2D plane", "false",
                          Patterns::Bool (),
                          "Whether to use a 2d data slice of a 3d data file "
                          "or the entire data file. Slicing a 3d dataset is "
                          "only supported for 2d models.");
        prm.declare_entry ("First point on slice", "0.0,1.0,0.0",
                           Patterns::Anything (),
                           "Point that determines the plane in which the 2d slice lies in. "
                           "This variable is only used if 'Slice dataset in 2d plane' is true. "
                           "The slice will go through this point, the point defined by the "
                           "parameter 'Second point on slice', and the center of the model "
                           "domain. After the rotation, this first point will lie along the "
                           "(0,1,0) axis of the coordinate system. The coordinates of the "
                           "point have to be given in Cartesian coordinates.");
        prm.declare_entry ("Second point on slice", "1.0,0.0,0.0",
                           Patterns::Anything (),
                           "Second point that determines the plane in which the 2d slice lies in. "
                           "This variable is only used if 'Slice dataset in 2d plane' is true. "
                           "The slice will go through this point, the point defined by the "
                           "parameter 'First point on slice', and the center of the model "
                           "domain. The coordinates of the point have to be given in Cartesian "
                           "coordinates.");
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    AsciiDataInitial<dim>::parse_parameters (ParameterHandler &prm,
                                             const std::string &subsection_name)
    {
      Utilities::AsciiDataBase<dim>::parse_parameters(prm,
                                                      subsection_name);

      prm.enter_subsection(subsection_name);
      {
        slice_data = prm.get_bool ("Slice dataset in 2D plane");
        if (slice_data == true)
          {
            AssertThrow (Plugins::plugin_type_matches<const GeometryModel::SphericalShell<dim>> (this->get_geometry_model()),
                         ExcMessage ("Slicing an ascii dataset is only supported when using "
                                     "a spherical shell geometry."));

            AssertThrow(dim == 2,
                        ExcMessage("The AsciiDataInitial class can only slice data in 2d models."));

            std::vector<double> point_one = Utilities::string_to_double(Utilities::split_string_list(prm.get("First point on slice")));
            std::vector<double> point_two = Utilities::string_to_double(Utilities::split_string_list(prm.get("Second point on slice")));

            AssertThrow(point_one.size() == 3 && point_two.size() == 3,
                        ExcMessage("The points on the slice in the Ascii data model need "
                                   "to be given in three dimensions; in other words, as three "
                                   "numbers, separated by commas."));

            Point<3> first_point_on_slice;
            Point<3> second_point_on_slice;

            for (unsigned int d=0; d<3; ++d)
              first_point_on_slice[d] = point_one[d];

            for (unsigned int d=0; d<3; ++d)
              second_point_on_slice[d] = point_two[d];

            rotation_matrix = Utilities::compute_rotation_matrix_for_slice(first_point_on_slice, second_point_on_slice);
          }
      }
      prm.leave_subsection();
    }



    template <int dim>
    AsciiDataProfile<dim>::AsciiDataProfile ()
      = default;



    template <int dim>
    void
    AsciiDataProfile<dim>::initialize (const MPI_Comm communicator)
    {
      lookup = std::make_unique<Utilities::StructuredDataLookup<1>> (this->scale_factor);

      const std::string filename = this->data_directory + this->data_file_name;

      AssertThrow(Utilities::fexists(filename, communicator) || filename_is_url(filename),
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
    AsciiDataProfile<dim>::get_interpolation_point_coordinates() const
    {
      return lookup->get_interpolation_point_coordinates(0);
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
