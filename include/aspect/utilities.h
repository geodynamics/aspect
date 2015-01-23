/*
  Copyright (C) 2014 by the authors of the ASPECT code.

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


#ifndef __aspect__utilities_h
#define __aspect__utilities_h

#include <aspect/global.h>

#include <deal.II/base/std_cxx1x/array.h>
#include <deal.II/base/point.h>
#include <deal.II/base/conditional_ostream.h>

#include <deal.II/base/table_indices.h>
#include <deal.II/base/function_lib.h>

#include <aspect/geometry_model/interface.h>
#include <aspect/simulator_access.h>



namespace aspect
{
  /**
   * A namespace for utility functions that might be used in many different
   * places to prevent code duplication.
   */
  namespace Utilities
  {
    using namespace dealii;

    /**
     * Returns spherical coordinates of a cartesian point. The returned array
     * is filled with radius, phi and theta (polar angle). If the dimension is
     * set to 2 theta is omitted. Phi is always normalized to [0,2*pi].
     *
     */
    template <int dim>
    std_cxx1x::array<double,dim>
    spherical_coordinates(const Point<dim> &position);

    /**
     * Return the cartesian point of a spherical position defined by radius,
     * phi and theta (polar angle). If the dimension is set to 2 theta is
     * omitted.
     */
    template <int dim>
    Point<dim>
    cartesian_coordinates(const std_cxx1x::array<double,dim> &scoord);


    template <int dim, int grid_dim>
    class AsciiDataLookup
    {
      public:
        AsciiDataLookup(const GeometryModel::Interface<dim> &geometry_model,
                        const unsigned int components,
                        const double scale_factor,
                        const types::boundary_id boundary_id = numbers::invalid_boundary_id);

        /**
         * Checks whether a file named filename exists.
         *
         * @param filename File to check existence
         */
        bool fexists(const std::string &filename);

        /**
         * Outputs the AsciiData module information at model start.
         */
        void screen_output(const ConditionalOStream &pcout) const;

        /**
         * Loads a data text file. Throws an exception if the file does not exist,
         * if the data file format is incorrect or if the file grid changes over model runtime.
         */
        void
        load_file(const std::string &filename);

        /**
         * Returns the computed data (velocity, temperature, etc. - according to the used plugin)
         * in cartesian coordinates.
         *
         * @param position The current position to compute the data (velocity, temperature, etc.)
         * @param time_weight A weighting between the two current timesteps n and n+1
         */
        double
        get_data(const Point<dim> &position,
                 const unsigned int component,
                 const double time_weight) const;

      private:
        /**
         * The number of data components read in (=columns in the data file).
         */
        const unsigned int components;

        /**
         * A reference to the geometry model. Is needed to convert the
         * position into spherical coordinates if necessary.
         */
        const GeometryModel::Interface<dim> &geometry_model;

        /**
         * Interpolation functions to access the data.
         */
        std::vector<Functions::InterpolatedUniformGridData<grid_dim> *> data;
        std::vector<Functions::InterpolatedUniformGridData<grid_dim> *> old_data;

        /**
         * Model size
         */
        std_cxx11::array<std::pair<double,double>,grid_dim> grid_extent;

        /**
         * Number of points in the data grid.
         */
        TableIndices<grid_dim> table_points;

        /**
         * Dimensions of the boundary plane
         */
        unsigned int boundary_dimensions[grid_dim];

        /**
         * Scales the data boundary condition by a scalar factor. Can be
         * used to transform the unit of the data.
         */
        const double scale_factor;


        /**
         * Gets the extents of the model in the relevant dimensions and returns
         * the according minimum and maximum value of the boundary in each dimension.
         * In case of a spherical shell the function returns the values in
         * spherical coordinates.
         */
        std_cxx11::array<std::pair<double,double>,dim>
        get_model_extent (const GeometryModel::Interface<dim> &geometry_model) const;


        /**
         * Determines which of the dimensions of the position is used to find
         * the data point in the data grid. E.g. the left boundary of a box model extents in
         * the y and z direction (position[1] and position[2]), therefore the function
         * would return [1,2] for dim==3 or [1] for dim==2. We are lucky that these indices are
         * identical for the box and the spherical shell (if we use spherical coordinates for the
         * spherical shell), therefore we do not need to distinguish between them. For the initial
         * condition this function is trivial, because the position in the data grid is the same as
         * the actual position (the function returns [0,1,2] or [0,1]), but for the boundary
         * conditions it matters.
         */
        std_cxx11::array<unsigned int,grid_dim>
        get_boundary_dimensions (const types::boundary_id boundary_id) const;


        /**
         * Computes the table indices of each entry in the input data file.
         * The index depends on dim, grid_dim and the number of components.
         */
        TableIndices<grid_dim>
        compute_table_indices(const unsigned int i) const;

    };

    template <int dim>
    class AsciiDataBase : public SimulatorAccess<dim>
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
                            const std::string &default_filename);

        /**
         * Read the parameters from the parameter file.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);

        /**
         * Directory in which the data files are present.
         */
        std::string data_directory;

        /**
         * Filename of data file. The file names can contain
         * the specifiers %s and/or %c (in this order), meaning the name of the
         * boundary and the number of the data file time step.
         */
        std::string data_file_name;

        /**
         * Scale the data by a scalar factor. Can be
         * used to transform the unit of the data (if they are not
         * specified in SI units (m/s or m/yr depending on the
         * "Use years in output instead of seconds" parameter).
         */
        double scale_factor;
    };

    template <int dim>
    class AsciiDataBoundary : public AsciiDataBase<dim>
    {
      public:
        /**
         * Constructor
         */
        AsciiDataBoundary();

      protected:

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
         * necessary and outputs a warning if the end of the set of data
         * files is reached.
         */
        virtual
        void
        update();

        double
        get_data_component (const types::boundary_id             boundary_indicator,
                            const Point<dim>                    &position,
                            const unsigned int                   component) const;

        /**
         * A variable that stores the currently used data file of a
         * series. It gets updated if necessary by update_data().
         */
        int  current_file_number;

        /**
         * Time from which on the data file with number 'First data
         * file number' is used as boundary condition. Previous to this
         * time, 0 is returned for every field. Depending on the setting
         * of the global 'Use years in output instead of seconds' flag
         * in the input file, this number is either interpreted as seconds or as years."
         */
        double first_data_file_model_time;

        /**
         * Number of the first data file to be loaded when the model time
         * is larger than 'First data file model time'.
         */
        int first_data_file_number;

        /**
         * In some cases the boundary files are not numbered in increasing
         * but in decreasing order (e.g. 'Ma BP'). If this flag is set to
         * 'True' the plugin will first load the file with the number
         * 'First data file number' and decrease the file number during
         * the model run.
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
         * off after finding no more data files to suppress attempts to
         * read in new files.
         */
        bool time_dependent;

        /**
         * Map between the boundary id and an object that reads and processes
         * data we get from text files.
         */
        std::map<types::boundary_id,
                  std_cxx11::shared_ptr<::aspect::Utilities::AsciiDataLookup<dim,dim-1> > > lookups;

        /**
         * Handles the update of the data in lookup.
         */
        void
        update_data (const types::boundary_id boundary_id);

        /**
         * Handles settings and user notification in case the time-dependent
         * part of the boundary condition is over.
         */
        void
        end_time_dependence (const int timestep,
                             const types::boundary_id boundary_id);

        /**
         * Create a filename out of the name template.
         */
        std::string
        create_filename (const int timestep,
                         const types::boundary_id boundary_id) const;


        /**
         * Declare the parameters all derived classes take from input files.
         */
        static
        void
        declare_parameters (ParameterHandler  &prm,
                            const std::string &default_directory,
                            const std::string &default_filename);

        /**
         * Read the parameters from the parameter file.
         */
        void
        parse_parameters (ParameterHandler &prm);


    };

  }
}

#endif
