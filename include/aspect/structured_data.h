/*
  Copyright (C) 2014 - 2020 by the authors of the ASPECT code.

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
    using namespace dealii;

    /**
     * StructuredDataLookup (formerly AsciiDataLookup) represents structured
     * data that can be read from files including in ascii format.
     *
     * For ascii files the files need to be formated as follows:
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
        StructuredDataLookup(const unsigned int components,
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
         * The data in @p raw_data consists of a Table for each of the @p n_components components.
         */
        void reinit(const std::vector<std::string> &column_names,
                    const std::vector<std::vector<double>> &coordinate_values,
                    const std::vector<Table<dim,double> > &raw_data
                   );

        /**
         * Loads a data text file. Throws an exception if the file does not
         * exist, if the data file format is incorrect or if the file grid
         * changes over model runtime.
         */
        void
        load_file(const std::string &filename,
                  const MPI_Comm &communicator);

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

      private:
        /**
         * The number of data components read in (=columns in the data file).
         */
        unsigned int components;

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
         * Computes the table indices given the size @p sizes of the
         * i-th entry.
         */
        TableIndices<dim>
        compute_table_indices(const TableIndices<dim> &sizes, const unsigned int i) const;
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
         * Determines which of the dimensions of the position is used to find
         * the data point in the data grid. E.g. the left boundary of a box
         * model extents in the y and z direction (position[1] and
         * position[2]), therefore the function would return [1,2] for dim==3
         * or [1] for dim==2. We are lucky that these indices are identical
         * for the box and the spherical shell (if we use spherical
         * coordinates for the spherical shell), therefore we do not need to
         * distinguish between them. For the initial condition this function
         * is trivial, because the position in the data grid is the same as
         * the actual position (the function returns [0,1,2] or [0,1]), but
         * for the boundary conditions it matters.
         */
        std::array<unsigned int,dim-1>
        get_boundary_dimensions (const types::boundary_id boundary_id) const;

        /**
         * A variable that stores the currently used data file of a series. It
         * gets updated if necessary by update().
         */
        int current_file_number;

        /**
         * Time from which on the data file with number 'First data file
         * number' is used as boundary condition. Previous to this time, 0 is
         * returned for every field. Depending on the setting of the global
         * 'Use years in output instead of seconds' flag in the input file,
         * this number is either interpreted as seconds or as years."
         */
        double first_data_file_model_time;

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
            std::unique_ptr<aspect::Utilities::StructuredDataLookup<dim-1> > > lookups;

        /**
         * Map between the boundary id and the old data objects.
         */
        std::map<types::boundary_id,
            std::unique_ptr<aspect::Utilities::StructuredDataLookup<dim-1> > > old_lookups;

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
        get_data_component (const Point<dim>                    &position,
                            const unsigned int                   component) const;

      protected:
        /**
         * Pointer to an object that reads and processes data we get from text
         * files.
         */
        std::unique_ptr<aspect::Utilities::StructuredDataLookup<dim> > lookup;
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
        std::vector<std::unique_ptr<aspect::Utilities::StructuredDataLookup<dim-1> >> lookups;

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
        initialize (const MPI_Comm &communicator);


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
        std::unique_ptr<aspect::Utilities::StructuredDataLookup<1> > lookup;
    };



    template<int dim>
    using AsciiDataLookup DEAL_II_DEPRECATED = StructuredDataLookup<dim>;
  }
}

#endif
