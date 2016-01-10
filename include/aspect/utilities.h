/*
  Copyright (C) 2014, 2015 by the authors of the ASPECT code.

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

#include <deal.II/base/std_cxx11/array.h>
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
    using namespace dealii::Utilities;

    /**
     * Returns spherical coordinates of a cartesian point. The returned array
     * is filled with radius, phi and theta (polar angle). If the dimension is
     * set to 2 theta is omitted. Phi is always normalized to [0,2*pi].
     *
     */
    template <int dim>
    std_cxx11::array<double,dim>
    spherical_coordinates(const Point<dim> &position);

    /**
     * Return the cartesian point of a spherical position defined by radius,
     * phi and theta (polar angle). If the dimension is set to 2 theta is
     * omitted.
     */
    template <int dim>
    Point<dim>
    cartesian_coordinates(const std_cxx11::array<double,dim> &scoord);

    /**
     * Checks whether a file named filename exists.
     *
     * @param filename File to check existence
     */
    bool fexists(const std::string &filename);

    /**
     * A namespace defining the cubic spline interpolation that can be used
     * between different spherical layers in the mantle.
     */
    namespace tk
    {
      // band matrix solver
      class band_matrix
      {
        private:
          std::vector< std::vector<double> > m_upper;  // upper band
          std::vector< std::vector<double> > m_lower;  // lower band
        public:
          band_matrix() {};                             // constructor
          band_matrix(int dim, int n_u, int n_l);       // constructor
          ~band_matrix() {};                            // destructor
          void resize(int dim, int n_u, int n_l);       // init with dim,n_u,n_l
          int dim() const;                              // matrix dimension
          int num_upper() const
          {
            return m_upper.size()-1;
          }
          int num_lower() const
          {
            return m_lower.size()-1;
          }
          // access operator
          double &operator () (int i, int j);             // write
          double   operator () (int i, int j) const;      // read
          // we can store an additional diogonal (in m_lower)
          double &saved_diag(int i);
          double  saved_diag(int i) const;
          void lu_decompose();
          std::vector<double> r_solve(const std::vector<double> &b) const;
          std::vector<double> l_solve(const std::vector<double> &b) const;
          std::vector<double> lu_solve(const std::vector<double> &b,
                                       bool is_lu_decomposed=false);
      };
      // spline interpolation
      class spline
      {
        private:
          std::vector<double> m_x,m_y;           // x,y coordinates of points
          // interpolation parameters
          // f(x) = a*(x-x_i)^3 + b*(x-x_i)^2 + c*(x-x_i) + y_i
          std::vector<double> m_a,m_b,m_c,m_d;
        public:
          void set_points(const std::vector<double> &x,
                          const std::vector<double> &y, bool cubic_spline=true);
          double operator() (double x) const;
      };
    }

    /**
     * Extract the compositional values at a single quadrature point with
     * index @p q from @p composition_values, which is indexed by
     * compositional index and quadrature point, and write them into @p
     * composition_values_at_q_point. In other words,
     * this extracts @p composition_values[i][q] for all @p i.
     */
    inline
    void
    extract_composition_values_at_q_point (const std::vector<std::vector<double> > &composition_values,
                                           const unsigned int q,
                                           std::vector<double> &composition_values_at_q_point)
    {
      Assert(q<composition_values.size(), ExcInternalError());
      Assert(composition_values_at_q_point.size() > 0,
             ExcInternalError());

      for (unsigned int k=0; k < composition_values_at_q_point.size(); ++k)
        {
          Assert(composition_values[k].size() == composition_values_at_q_point.size(),
                 ExcInternalError());
          composition_values_at_q_point[k] = composition_values[k][q];
        }
    }

    /**
     * Provide an object of type T filled with a signaling NaN that will cause an exception
     * when used in a computation. This basically serves the purpose of creating an object
     * that is not initialized.
     **/
    template <class T>
    T
    signaling_nan();

    template <>
    inline
    double
    signaling_nan<double>()
    {
      return std::numeric_limits<double>::signaling_NaN();
    }

    template <>
    inline
    SymmetricTensor<2,2>
    signaling_nan<SymmetricTensor<2,2> >()
    {
      const unsigned int dim = 2;
      SymmetricTensor<2,dim> nan_tensor;
      for (unsigned int i=0; i<dim; ++i)
        for (unsigned int j=0; j<dim; ++j)
          nan_tensor[i][j] = std::numeric_limits<double>::signaling_NaN();
      return nan_tensor;
    }
    template <>
    inline
    SymmetricTensor<2,3>
    signaling_nan<SymmetricTensor<2,3> >()
    {
      const unsigned int dim = 3;
      SymmetricTensor<2,dim> nan_tensor;
      for (unsigned int i=0; i<dim; ++i)
        for (unsigned int j=0; j<dim; ++j)
          nan_tensor[i][j] = std::numeric_limits<double>::signaling_NaN();
      return nan_tensor;
    }
    template <>
    inline
    Tensor<2,2>
    signaling_nan<Tensor<2,2> >()
    {
      const unsigned int dim = 2;
      Tensor<2,dim> nan_tensor;
      for (unsigned int i=0; i<dim; ++i)
        for (unsigned int j=0; j<dim; ++j)
          nan_tensor[i][j] = std::numeric_limits<double>::signaling_NaN();
      return nan_tensor;
    }
    template <>
    inline
    Tensor<2,3>
    signaling_nan<Tensor<2,3> >()
    {
      const unsigned int dim = 3;
      Tensor<2,dim> nan_tensor;
      for (unsigned int i=0; i<dim; ++i)
        for (unsigned int j=0; j<dim; ++j)
          nan_tensor[i][j] = std::numeric_limits<double>::signaling_NaN();
      return nan_tensor;
    }
    template <>
    inline
    Tensor<1,2>
    signaling_nan<Tensor<1,2> >()
    {
      const unsigned int dim = 2;
      Tensor<1,dim> nan_tensor;
      for (unsigned int i=0; i<dim; ++i)
        nan_tensor[i] = std::numeric_limits<double>::signaling_NaN();
      return nan_tensor;
    }
    template <>
    inline
    Tensor<1,3>
    signaling_nan<Tensor<1,3> >()
    {
      const unsigned int dim = 3;
      Tensor<1,dim> nan_tensor;
      for (unsigned int i=0; i<dim; ++i)
        nan_tensor[i] = std::numeric_limits<double>::signaling_NaN();
      return nan_tensor;
    }


    /**
     * AsciiDataLookup reads in files containing input data in ascii format.
     * Note the required format of the input data: The first lines may contain
     * any number of comments if they begin with '#', but one of these lines
     * needs to contain the number of grid points in each dimension as for
     * example '# POINTS: 3 3'. The order of the columns has to be
     * 'coordinates data' with @p dim coordinate columns and @p components
     * data columns. Note that the data in the input files need to be sorted
     * in a specific order: the first coordinate needs to ascend first,
     * followed by the second and so on in order to assign the correct data to
     * the prescribed coordinates.
     */
    template <int dim>
    class AsciiDataLookup
    {
      public:
        AsciiDataLookup(const unsigned int components,
                        const double scale_factor);

        /**
         * Loads a data text file. Throws an exception if the file does not
         * exist, if the data file format is incorrect or if the file grid
         * changes over model runtime.
         */
        void
        load_file(const std::string &filename);

        /**
         * Returns the computed data (velocity, temperature, etc. - according
         * to the used plugin) in cartesian coordinates.
         *
         * @param position The current position to compute the data (velocity,
         * temperature, etc.)
         * @param component The index of the data column to be returned.
         */
        double
        get_data(const Point<dim> &position,
                 const unsigned int component) const;

      private:
        /**
         * The number of data components read in (=columns in the data file).
         */
        const unsigned int components;

        /**
         * Interpolation functions to access the data.
         */
        std::vector<Functions::InterpolatedUniformGridData<dim> *> data;

        /**
         * Min and Max coordinates in data file
         */
        std_cxx11::array<std::pair<double,double>,dim> grid_extent;

        /**
         * Number of points in the data grid.
         */
        TableIndices<dim> table_points;

        /**
         * Scales the data boundary condition by a scalar factor. Can be used
         * to transform the unit of the data.
         */
        const double scale_factor;

        /**
         * Computes the table indices of each entry in the input data file.
         * The index depends on dim, grid_dim and the number of components.
         */
        TableIndices<dim>
        compute_table_indices(const unsigned int i) const;

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
                            const std::string &default_filename);

        /**
         * Read the parameters from the parameter file.
         */
        void
        parse_parameters (ParameterHandler &prm);

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
        std_cxx11::array<unsigned int,dim-1>
        get_boundary_dimensions (const types::boundary_id boundary_id) const;

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
         * A variable that stores the currently used data file of a series. It
         * gets updated if necessary by update().
         */
        int  current_file_number;

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
            std_cxx11::shared_ptr<aspect::Utilities::AsciiDataLookup<dim-1> > > lookups;

        /**
         * Map between the boundary id and the old data objects.
         */
        std::map<types::boundary_id,
            std_cxx11::shared_ptr<aspect::Utilities::AsciiDataLookup<dim-1> > > old_lookups;

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
        std_cxx11::shared_ptr<aspect::Utilities::AsciiDataLookup<dim> > lookup;
    };
  }
}

#endif
