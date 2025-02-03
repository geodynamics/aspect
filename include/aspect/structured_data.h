/*
  Copyright (C) 2014 - 2024 by the authors of the ASPECT code.

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


#ifndef _aspect_structured_data_h
#define _aspect_structured_data_h

#include <aspect/global.h>
#include <aspect/simulator_access.h>

#include <array>

namespace aspect
{
  namespace Utilities
  {
    /**
     * Because many places in ASPECT assume that all functions in the namespace
     * <code>dealii::Utilities</code> are available without qualification as
     * <code>Utilities::function</code>, just as all the function in the
     * namespace <code>aspect::Utilities</code>, we make sure all these functions
     * are available inside <code>aspect::Utilities</code>. This is maybe not
     * the cleanest solution, but it is most compatible with a lot of existing
     * code, and also allows to migrate ASPECT functions into deal.II when
     * useful without introducing incompatibilities.
     *
     * We need to do this in every header that introduces something into the
     * namespace <code>aspect::Utilities</code>, because it needs to happen
     * no matter which header files of ASPECT are included.
     */
    using namespace dealii::Utilities;

    /**
     * StructuredDataLookup (formerly AsciiDataLookup) represents structured
     * data that can be read from files including in ascii format.
     *
     * For ascii files the files need to be formatted as follows:
     * Note the required format of the input data: The first lines may contain
     * any number of comments if they begin with '#', but one of these lines
     * needs to contain the number of grid points in each dimension as for
     * example '# POINTS: 3 3'. The comments can optionally be followed by a
     * single line, which does not start with '#', containing the names of
     * the data columns.
     * The order of the following data columns has to be
     * 'coordinates data' with @p dim coordinate columns and @p components
     * data columns. Note that the data in the input files need to be sorted
     * in a specific order: the first coordinate needs to ascend first,
     * followed by the second and so on in order to assign the correct data to
     * the prescribed coordinates. The coordinates do not need to be
     * equidistant.
     */
    template <int dim>
    class StructuredDataLookup
    {
      public:
        /**
         * Constructor that explicitly prescribes the number of data columns
         * in the data file. If a list of data components is provided in the
         * data file it is checked that the length of this list is consistent
         * with this number of components. This constructor is mostly provided
         * for backwards compatibility. Not prescribing the number of components
         * and instead reading them from the input file allows for more
         * flexible files.
         */
        StructuredDataLookup(const unsigned int n_components,
                             const double scale_factor);

        /**
         * This constructor relies on the list of column names at the beginning
         * of the model file to determine the number of data components,
         * therefore when using this constructor it is necessary to provide
         * this list in the first uncommented line of the data file.
         */
        explicit StructuredDataLookup(const double scale_factor);

        /**
         * Replace the data stored in this class by the data given to this function.
         *
         * @note This is a manual way to fill the data. Consider load_file() if your data
         * is stored in a txt/csv file.
         *
         * The data consists of @p n_components (implicitly given by the size of @p column_names)
         * specified at @p coordinate_values[d] points in each of the dim coordinate directions @p d.
         *
         * The data in @p data_table consists of a Table for each of the @p n_components components.
         *
         * The `coordinate_values` and `data_table` arguments are rvalue references,
         * and the function will move the data so provided into another storage
         * location. In other words, after the call, the variables passed as the
         * second and third arguments may be empty or otherwise altered.
         *
         * This class
         * is able to share data between processes located on the same
         * machine if the last argument `root_process` is
         * set to anything other than `numbers::invalid_unsigned_int`. Then only the
         * indicated root process within the given MPI communicator needs to
         * pass in valid arguments whereas on all other processes the values
         * of `column_names`, `coordinate_values`, and `data_table` are
         * ignored. Instead, the indicated root process will make sure that
         * the other processes obtain valid data, for example by sharing
         * data tables in shared memory spaces. This reduces both the
         * computational effort in reading data from disk and parsing it on
         * all processes, and can vastly reduce the amount of memory required
         * on machines where many processor cores have access to shared memory
         * resources.
         *
         * If `root_process` equals `numbers::invalid_unsigned_int`, then
         * every process needs to pass data for all arguments and the
         * `mpi_communicator` argument is ignored.
         */
        void reinit(const std::vector<std::string> &column_names,
                    std::vector<std::vector<double>> &&coordinate_values,
                    std::vector<Table<dim,double>> &&data_table,
                    const MPI_Comm mpi_communicator = MPI_COMM_SELF,
                    const unsigned int root_process = numbers::invalid_unsigned_int);

        /**
         * Loads a data text file replacing the current data.
         *
         * Throws an exception if the file does not
         * exist, if the data file format is incorrect, or if the file data
         * like the number and names of columns change (if this class contains
         * existing data for example when used with a collection of ascii files
         * changing over the simulation time).
         *
         * This function supports normal text / ASCII files that can optionally
         * be compressed using .gz. If the given filename is a URL and libDAB
         * support is enabled, the data is downloaded and parsed from the ASCII
         * data received.
         *
         * This function uses the given MPI communicator to only load the data on
         * rank 0 and broadcasting the information to all other ranks.
         */
        void
        load_ascii(const std::string &filename,
                   const MPI_Comm communicator);

        /**
         * Fill the current object with data read from a NetCDF file
         * with filename @p filename. This call will fail if ASPECT is not
         * configured with NetCDF support.
         *
         * @p data_column_names specifies the list of data columns to load (in the specified order). If
         * an empty vector is passed, all columns will be loaded.
         */
        void
        load_netcdf(const std::string &filename, const std::vector<std::string> &data_column_names = {});


        /**
         * Loads data from a file replacing the current data.
         *
         * Throws an exception if the file does not
         * exist, if the data file format is incorrect, or if the file data
         * like the number and names of columns change (if this class contains
         * existing data for example when used with a collection of ascii files
         * changing over the simulation time).
         *
         * The following formats are currently supported:
         * - ASCII files (typically ending in .txt)
         * - gzip compressed ASCII files (ending in .gz)
         * - URLs starting with "http" (handled by libDAB)
         */
        void
        load_file(const std::string &filename,
                  const MPI_Comm communicator);

        /**
         * Returns the computed data (velocity, temperature, etc. - according
         * to the used plugin) in Cartesian coordinates.
         *
         * @param position The current position to compute the data (velocity,
         * temperature, etc.)
         * @param component The index (starting at 0) of the data column to be
         * returned. The index is therefore less than the number of data
         * columns in the data file (or specified in the constructor).
         */
        double
        get_data(const Point<dim> &position,
                 const unsigned int component) const;

        /**
         * Returns the gradient of the function based on the bilinear
         * interpolation of the data (velocity, temperature, etc. - according
         * to the used plugin) in Cartesian coordinates.
         *
         * @param position The current position to compute the data (velocity,
         * temperature, etc.)
         * @param component The index of the data column to be returned.
         */
        Tensor<1,dim>
        get_gradients(const Point<dim> &position,
                      const unsigned int component);

        /**
         * Returns a vector that contains the names of all data columns in the
         * order of their appearance in the data file (and their order in the
         * memory data table). Returns an empty vector if no names are provided
         * or the file is not read in yet.
         */
        std::vector<std::string>
        get_column_names() const;

        /**
         * Returns whether the stored coordinates are equidistant. If
         * coordinates are equidistant the lookup is more efficient. Returns
         * false if no coordinates are loaded at the moment.
         */
        bool
        has_equidistant_coordinates() const;

        /**
         * Returns the coordinates of the interpolation points at which data is
         * stored. This function can be used to determine the number of data
         * points in each of the coordinate directions, or to query
         * data only at exactly the positions at which they are available (avoiding
         * interpolation).
         *
         * @param dimension The spatial direction for which to return the data
         * coordinates, e.g. 0 for $x$-direction, 1 for $y$-direction, or equivalent
         * values if your data coordinates are other dimensions such as
         * temperature, pressure.
         */
        const std::vector<double> &
        get_interpolation_point_coordinates(const unsigned int dimension) const;

        /**
         * Returns the column index of a column with the given name
         * @p column_name. Throws an exception if no such
         * column exists or no names were provided in the file.
         */
        unsigned int
        get_column_index_from_name(const std::string &column_name) const;

        /**
         * Returns a string that contains the name of the column with index
         * @p column_index. Throws an exception if no such
         * column exists or no name was provided in the file.
         */
        std::string
        get_column_name_from_index(const unsigned int column_index) const;

        /**
         * Return the maximum value of the component values.
         */
        double get_maximum_component_value(const unsigned int component) const;

        /**
         * Retrieve the number of table points for a given dimension.
         * Equivalent to calling get_interpolation_point_coordinates().size().
         *
         * @param dimension The index of the dimension for which to get the number of table points.
         * @return The number of points along the specified dimension.
         */
        unsigned int get_number_of_coordinates(const unsigned int dimension) const;

      private:
        /**
         * The number of data components read in (=columns in the data file).
         */
        unsigned int n_components;

        /**
         * The names of the data components in the columns of the read file.
         * Does not contain any strings if none are provided in the first
         * uncommented line of the file.
         */
        std::vector<std::string> data_component_names;

        /**
         * Interpolation functions to access the data.
         * Either InterpolatedUniformGridData or InterpolatedTensorProductGridData;
         * the type is determined from the grid specified in the data file.
         */
        std::vector<std::unique_ptr<Function<dim>>> data;

        /**
         * The coordinate values in each direction as specified in the data file.
         */
        std::array<std::vector<double>,dim> coordinate_values;

        /**
         * The maximum value of each component
         */
        std::vector<double> maximum_component_value;

        /**
         * Number of points in the data grid as specified in the data file.
         */
        TableIndices<dim> table_points;

        /**
         * Scales the data boundary condition by a scalar factor. Can be used
         * to transform the unit of the data.
         */
        const double scale_factor;

        /**
         * Stores whether the coordinate values are equidistant or not,
         * this determines the type of data function stored.
         */
        bool coordinate_values_are_equidistant;

        /**
         * Computes the table indices given the size @p sizes of the
         * entry with index @p idx.
         */
        TableIndices<dim>
        compute_table_indices(const TableIndices<dim> &sizes, const std::size_t idx) const;

    };

    /**
     * AsciDataBase is a generic plugin used for declaring and reading the
     * parameters from the parameter file.
     */
    template <int dim>
    class AsciiDataBase
    {
      public:
        /**
         * Constructor
         */
        AsciiDataBase();

        /**
         * Declare the parameters all derived classes take from input files.
         */
        static
        void
        declare_parameters (ParameterHandler  &prm,
                            const std::string &default_directory,
                            const std::string &default_filename,
                            const std::string &subsection_name = "Ascii data model");

        /**
         * Read the parameters from the parameter file.
         */
        void
        parse_parameters (ParameterHandler &prm,
                          const std::string &subsection_name = "Ascii data model");

        /**
         * Directory in which the data files are present.
         */
        std::string data_directory;

        /**
         * Filename of data file. The file names can contain the specifiers %s
         * and/or %c (in this order), meaning the name of the boundary and the
         * number of the data file time step.
         */
        std::string data_file_name;

        /**
         * Scale the data by a scalar factor. Can be used to transform the
         * unit of the data (if they are not specified in SI units (m/s or
         * m/yr depending on the "Use years in output instead of seconds"
         * parameter).
         */
        double scale_factor;
    };

    /**
     * A base class that implements boundary conditions determined from a
     * AsciiData input file.
     */
    template <int dim>
    class AsciiDataBoundary : public Utilities::AsciiDataBase<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor
         */
        AsciiDataBoundary();

        /**
         * Initialization function. This function is called once at the
         * beginning of the program. Checks preconditions.
         */
        virtual
        void
        initialize (const std::set<types::boundary_id> &boundary_ids,
                    const unsigned int components);

        /**
         * A function that is called at the beginning of each time step. For
         * the current plugin, this function loads the next data files if
         * necessary and outputs a warning if the end of the set of data files
         * is reached.
         */
        void
        update();

        /**
         * Returns the data component at the given position.
         */
        double
        get_data_component (const types::boundary_id             boundary_indicator,
                            const Point<dim>                    &position,
                            const unsigned int                   component) const;

        /**
         * Returns the maximum value of the given data component.
         */
        double
        get_maximum_component_value (const types::boundary_id boundary_indicator,
                                     const unsigned int       component) const;

        /**
         * Return the gradients of the parameters from the parameter file.
         */
        Tensor<1,dim-1>
        vector_gradient(const types::boundary_id boundary_indicator,
                        const Point<dim>        &p,
                        const unsigned int       component) const;

        /**
         * Declare the parameters all derived classes take from input files.
         *
         * @param prm The parameter handler in which the parameters are declared.
         * @param default_directory The default value for the data directory parameter.
         * @param default_filename The default value for the filename parameter.
         * @param subsection_name The name of the parameter file subsection all
         * parameters will be declared in. The function will enter this subsection,
         * declare all parameters, then leave this subsection, to return @p prm in
         * the same subsection it was in before.
         * @param declare_time_dependent_parameters Whether to declare the parameter
         * that are only needed for time dependent AsciiDataBoundary objects. If
         * the caller already knows time dependence is not supported for the current
         * application, disabling this parameter avoids introducing these parameters
         * to the parameter handler.
         */
        static
        void
        declare_parameters (ParameterHandler  &prm,
                            const std::string &default_directory,
                            const std::string &default_filename,
                            const std::string &subsection_name = "Ascii data model",
                            const bool declare_time_dependent_parameters = true);

        /**
         * Read the parameters from the parameter file.
         *
         * @param prm The parameter handler from which the parameters are parsed.
         * @param subsection_name The name of the parameter file subsection all
         * parameters will be parsed from. The function will enter this subsection,
         * parse all parameters, then leave this subsection, to return @p prm in
         * the same subsection it was in before.
         * @param parse_time_dependent_parameters Whether to parse the parameter
         * that are only needed for time dependent AsciiDataBoundary objects. This
         * parameter always needs to be set to the same value that was handed over
         * to declare_parameters().
         */
        void
        parse_parameters (ParameterHandler &prm,
                          const std::string &subsection_name = "Ascii data model",
                          const bool parse_time_dependent_parameters = true);

      protected:
        /**
         * A variable that stores the currently used data file of a series. It
         * gets updated if necessary by update().
         */
        int current_file_number;

        /**
         * Number of the first data file to be loaded when the model time is
         * larger than 'First data file model time'.
         */
        int first_data_file_number;

        /**
         * In some cases the boundary files are not numbered in increasing but
         * in decreasing order (e.g. 'Ma BP'). If this flag is set to 'True'
         * the plugin will first load the file with the number 'First data
         * file number' and decrease the file number during the model run.
         */
        bool decreasing_file_order;

        /**
         * Time in model units (depends on other model inputs) between two
         * data files.
         */
        double data_file_time_step;

        /**
         * Weight between data file n and n+1 while the current time is
         * between the two values t(n) and t(n+1).
         */
        double time_weight;

        /**
         * State whether we have time_dependent boundary conditions. Switched
         * off after finding no more data files to suppress attempts to read
         * in new files.
         */
        bool time_dependent;

        /**
         * Map between the boundary id and an object that reads and processes
         * data we get from text files.
         */
        std::map<types::boundary_id,
            std::unique_ptr<aspect::Utilities::StructuredDataLookup<dim-1>>> lookups;

        /**
         * Map between the boundary id and the old data objects.
         */
        std::map<types::boundary_id,
            std::unique_ptr<aspect::Utilities::StructuredDataLookup<dim-1>>> old_lookups;

        /**
         * Handles the update of the data in lookup.
         */
        void
        update_data (const types::boundary_id boundary_id,
                     const bool reload_both_files);

        /**
         * Handles settings and user notification in case the time-dependent
         * part of the boundary condition is over.
         */
        void
        end_time_dependence ();

        /**
         * Create a filename out of the name template.
         */
        std::string
        create_filename (const int filenumber,
                         const types::boundary_id boundary_id) const;
    };



    /**
     * A base class that implements initial conditions determined from a
     * AsciiData input file.
     */
    template <int dim>
    class AsciiDataInitial : public Utilities::AsciiDataBase<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor
         */
        AsciiDataInitial();

        /**
         * Initialization function. This function is called once at the
         * beginning of the program. Checks preconditions.
         */
        virtual
        void
        initialize (const unsigned int components);


        /**
         * Returns the data component at the given position.
         */
        double
        get_data_component (const Point<dim> &position,
                            const unsigned int component) const;

        /**
         * Declare the parameters all derived classes take from input files.
         */
        static
        void
        declare_parameters (ParameterHandler  &prm,
                            const std::string &default_directory,
                            const std::string &default_filename,
                            const std::string &subsection_name = "Ascii data model");

        /**
         * Read the parameters from the parameter file.
         */
        void
        parse_parameters (ParameterHandler &prm,
                          const std::string &subsection_name = "Ascii data model");

      protected:
        /**
         * Pointer to an object that reads and processes data we get from text
         * files.
         */
        std::unique_ptr<aspect::Utilities::StructuredDataLookup<dim>> lookup;

        /**
         * Pointer to an object that reads and processes data we get from text
         * files if the current model is a slice of the input file (e.g. a 2D
         * model and a 3D data file).
         */
        std::unique_ptr<aspect::Utilities::StructuredDataLookup<3>> slice_lookup;

        /**
         * Whether to use a dataset that has the same spatial dimensions as
         * the model or not. If true only a 2D slice of a 3D dataset is used.
         */
        bool slice_data;

        /**
         * The matrix that describes the rotation by which a 2D model
         * needs to be transformed to a plane that contains the origin and
         * the two prescribed points given in the input.
         */
        Tensor<2,3> rotation_matrix;
    };


    /**
     * A base class that implements conditions determined from a
     * layered AsciiData input file.
     */
    template <int dim>
    class AsciiDataLayered : public Utilities::AsciiDataBase<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor
         */
        AsciiDataLayered();

        /**
         * Initialization function. This function is called once at the
         * beginning of the program. Checks preconditions.
         */
        virtual
        void
        initialize (const unsigned int components);


        /**
         * Returns the data component at the given position.
         */
        double
        get_data_component (const Point<dim> &position,
                            const unsigned int component) const;


        /**
         * Declare the parameters all derived classes take from input files.
         */
        static
        void
        declare_parameters (ParameterHandler  &prm,
                            const std::string &default_directory,
                            const std::string &default_filename,
                            const std::string &subsection_name = "Ascii data model");

        /**
         * Read the parameters from the parameter file.
         */
        void
        parse_parameters (ParameterHandler &prm,
                          const std::string &subsection_name = "Ascii data model");

      protected:
        /**
         * Pointer to an object that reads and processes data we get from text
         * files.
         */
        std::vector<std::unique_ptr<aspect::Utilities::StructuredDataLookup<dim-1>>> lookups;

      private:

        /**
         * Directory in which the data files are present.
         */
        std::string data_directory;

        /**
         * Filenames of data files.
         */
        std::vector<std::string> data_file_names;

        /**
         * Number of layer boundaries in the model.
         */
        unsigned int number_of_layer_boundaries;

        /**
         * Interpolation scheme for profile averaging.
         */
        std::string interpolation_scheme;


    };


    /**
     * A base class that reads in a data profile and provides its values.
     */
    template <int dim>
    class AsciiDataProfile : public Utilities::AsciiDataBase<dim>
    {
      public:
        /**
         * Constructor
         */
        AsciiDataProfile();

        /**
         * Initialization function. This function is called once at the
         * beginning of the program. Checks preconditions.
         */
        virtual
        void
        initialize (const MPI_Comm communicator);


        /**
         * Returns the data component at the given position.
         */
        double
        get_data_component (const Point<1>                      &position,
                            const unsigned int                   component) const;

        /**
         * Returns a vector that contains the names of all data columns in the
         * order of their appearance in the data file (and their order in the
         * memory data table). Returns an empty vector if no names are provided
         * or the file is not read in yet.
         */
        std::vector<std::string>
        get_column_names() const;

        /**
         * Returns the coordinates of the interpolation points at which data is
         * stored. This function can be used to determine the number of data
         * points in each of the coordinate directions, or to query
         * data only at exactly the positions at which they are available (avoiding
         * interpolation).
         *
         * Because this class represents a one-dimensional profile, the returned
         * values correspond to the values of the sole coordinate of the
         * interpolation points, in contrast to the
         * StructuredDataLookup::get_interpolation_point_coordinates()
         * function that takes an integer argument indicating which coordinate
         * is to be selected.
         */
        const std::vector<double> &
        get_interpolation_point_coordinates() const;

        /**
         * Returns the column index of a column with the given name
         * @p column_name. Throws an exception if no such
         * column exists or no names were provided in the file.
         */
        unsigned int
        get_column_index_from_name(const std::string &column_name) const;

        /**
         * Returns the column index of a column with the given name
         * @p column_name. Returns an invalid unsigned int if no such
         * column exists or no names were provided in the file.
         */
        unsigned int
        maybe_get_column_index_from_name(const std::string &column_name) const;

        /**
         * Returns a string that contains the name of the column with index
         * @p column_index. Returns an empty string if no such
         * column exists or no name was provided in the file.
         */
        std::string
        get_column_name_from_index(const unsigned int column_index) const;
      protected:
        /**
         * Pointer to an object that reads and processes data we get from text
         * files.
         */
        std::unique_ptr<aspect::Utilities::StructuredDataLookup<1>> lookup;
    };
  }
}

#endif
