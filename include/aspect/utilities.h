/*
  Copyright (C) 2014, 2015, 2016 by the authors of the ASPECT code.

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
     * Given an array @p values, consider three cases:
     * - If it has size @p N, return the original array.
     * - If it has size one, return an array of size @p N where all
     *   elements are equal to the one element of @p value.
     * - If it has any other size, throw an exception that uses
     *   @p id_text as an identifying string.
     *
     * This function is typically used for parameter lists that can either
     * contain different values for each of a set of objects (e.g., for
     * each compositional field), or contain a single value that is then
     * used for each object.
     */
    template <typename T>
    std::vector<T>
    possibly_extend_from_1_to_N (const std::vector<T> &values,
                                 const unsigned int N,
                                 const std::string &id_text);


    /**
     * Split the set of DoFs (typically locally owned or relevant) in @p whole_set into blocks
     * given by the @p dofs_per_block structure.
     *
     * The numbers of dofs per block need to add up to the size of the index space described
     * by @p whole_set.
     */
    void split_by_block (const std::vector<types::global_dof_index> &dofs_per_block,
                         const IndexSet &whole_set,
                         std::vector<IndexSet> &partitioned);

    namespace Coordinates
    {

      /**
       * Returns distance from the Earth's center, latitude and longitude from a
       * given ECEF Cartesian coordinates that account for ellipsoidal shape of
       * the Earth with WGS84 parameters.
       */
      template <int dim>
      std_cxx11::array<double,dim>
      WGS84_coordinates(const Point<dim> &position);

      /**
       * Returns spherical coordinates of a Cartesian point. The returned array
       * is filled with radius, phi and theta (polar angle). If the dimension is
       * set to 2 theta is omitted. Phi is always normalized to [0,2*pi].
       *
       */
      template <int dim>
      std_cxx11::array<double,dim>
      cartesian_to_spherical_coordinates(const Point<dim> &position);

      /**
       * Return the Cartesian point of a spherical position defined by radius,
       * phi and theta (polar angle). If the dimension is set to 2 theta is
       * omitted.
       */
      template <int dim>
      Point<dim>
      spherical_to_cartesian_coordinates(const std_cxx11::array<double,dim> &scoord);

      /**
       * Returns ellispoidal coordinates of a Cartesian point. The returned array
       * is filled with phi, theta and radius.
       *
       */
      template <int dim>
      std_cxx11::array<double,3>
      cartesian_to_ellipsoidal_coordinates(const Point<3> &position,
                                           const double semi_major_axis_a,
                                           const double eccentricity);

      /**
       * Return the Cartesian point of a ellispoidal position defined by phi,
       * phi and radius.
       */
      template <int dim>
      Point<3>
      ellipsoidal_to_cartesian_coordinates(const std_cxx11::array<double,3> &phi_theta_d,
                                           const double semi_major_axis_a,
                                           const double eccentricity);
    }

    /**
     * Given a vector @p v in @p dim dimensional space, return a set
     * of (dim-1) vectors that are orthogonal to @p v and to each
     * other. The lengths of these vectors equals that of the original
     * vector @p v to ensure a well-conditioned basis.
     */
    template <int dim>
    std_cxx11::array<Tensor<1,dim>,dim-1>
    orthogonal_vectors (const Tensor<1,dim> &v);

    /**
     * A function for evaluating real spherical harmonics. It takes the degree (l)
     * and the order (m) of the spherical harmonic, where l >= 0 and 0 <= m <=l.
     * It also takes the colatitude (theta) and longitude (phi), which are in
     * radians.
     *
     * There are an unfortunate number of normalization conventions in existence
     * for spherical harmonics. Here we use fully normalized spherical harmonics
     * including the Condon-Shortley phase. This corresponds to the definitions
     * given in equations B.72 and B.99-B.102 in Dahlen and Tromp (1998, ISBN: 9780691001241).
     * The functional form of the real spherical harmonic is given by
     *
     * \f[
     *    Y_{lm}(\theta, \phi) = \sqrt{2} X_{l \left| m \right| }(\theta) \cos m \phi \qquad \mathrm{if} \qquad -l \le m < 0
     * \f]
     * \f[
     *    Y_{lm}(\theta, \phi) = X_{l 0 }(\theta) \qquad \mathrm{if} \qquad m = 0
     * \f]
     * \f[
     *    Y_{lm}(\theta, \phi) = \sqrt{2} X_{lm}(\theta) \sin m \phi \qquad \mathrm{if}  \qquad 0< m \le m
     * \f]
     * where \f$X_{lm}( \theta )\f$ is an associated Legendre function.
     *
     * In practice it is often convenient to compute the sine (\f$-l \le m < 0\f$) and cosine (\f$0 < m \le l\f$)
     * variants of the real spherical harmonic at the same time. That is the approach taken
     * here, where we return a pair of numbers, the first corresponding the cosine part and the
     * second corresponding to the sine part. Given this, it is no longer necessary to distinguish
     * between postitive and negative \f$ m \f$, so this function only accepts \f$ m \ge 0 \f$.
     * For \f$ m = 0 \f$, there is only one part, which is stored in the first entry of the pair.
     *
     * @note This function uses the Boost spherical harmonics implementation internally,
     * which is not designed for very high order (> 100) spherical harmonics computation.
     * If you use spherical harmonics of a high order be sure to confirm the accuracy first.
     * For more information, see:
     * http://www.boost.org/doc/libs/1_49_0/libs/math/doc/sf_and_dist/html/math_toolkit/special/sf_poly/sph_harm.html
     */
    std::pair<double,double> real_spherical_harmonic( unsigned int l, //degree
                                                      unsigned int m, //order
                                                      double theta,   //colatitude (radians)
                                                      double phi );   //longitude (radians)

    /**
     * Checks whether a file named filename exists.
     *
     * @param filename File to check existence
     */
    bool fexists(const std::string &filename);

    /**
     * Reads the content of the ascii file @p filename on process 0 and
     * distributes the content by MPI_Bcast to all processes. The function
     * returns the content of the file on all processes.
     *
     * @param [in] filename The name of the ascii file to load.
     * @param [in] comm The MPI communicator in which the content is
     * distributed.
     * @return A string which contains the data in @p filename.
     */
    std::string
    read_and_distribute_file_content(const std::string &filename,
                                     const MPI_Comm &comm);

    /**
     * Creates a path as if created by the shell command "mkdir -p", therefore
     * generating directories from the highest to the lowest level if they are
     * not already existing.
     *
     * @param pathname String that contains the path to create. '/' is used as
     * directory separator.
     * @param mode Permissions (mode bits) of the created directories. See the
     * documentation of the chmod() command for more information.
     * @return The function returns the error value of the last mkdir call
     * inside. It returns zero on success. See the man page of mkdir() for
     * more information.
     */
    int
    mkdirp(std::string pathname, const mode_t mode = 0755);

    /**
     * Create directory @p pathname, optionally printing a message.
     *
     * @param pathname String that contains path to create. '/' is used as
     * directory separator.
     * @param comm MPI communicator, used to limit creation of directory to
     * processor 0.
     * @param silent Print a nicely formatted message on processor 0 if set
     * to true.
     */
    void create_directory(const std::string &pathname,
                          const MPI_Comm &comm,
                          bool silent);

    /**
     * A namespace defining the cubic spline interpolation that can be used
     * between different spherical layers in the mantle.
     */
    namespace tk
    {
      /**
       * Class for cubic spline interpolation
       */
      class spline
      {
        public:
          /**
           * Initialize the spline.
           *
           * @param x X coordinates of interpolation points.
           * @param y Values in the interpolation points.
           * @param cubic_spline Whether to construct a cubic spline or just do linear interpolation
           */
          void set_points(const std::vector<double> &x,
                          const std::vector<double> &y,
                          bool cubic_spline = true);
          /**
           * Evaluate spline at point @p x.
           */
          double operator() (double x) const;

        private:
          /**
           * x coordinates of points
           */
          std::vector<double> m_x;

          /**
           * interpolation parameters
           * \[
           * f(x) = a*(x-x_i)^3 + b*(x-x_i)^2 + c*(x-x_i) + y_i
           * \]
           */
          std::vector<double> m_a, m_b, m_c, m_y;
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

    template <typename T>
    inline
    std::vector<T>
    possibly_extend_from_1_to_N (const std::vector<T> &values,
                                 const unsigned int N,
                                 const std::string &id_text)
    {
      if (values.size() == 1)
        {
          return std::vector<T> (N, values[0]);
        }
      else if (values.size() == N)
        {
          return values;
        }
      else
        {
          // Non-specified behavior
          AssertThrow(false,
                      ExcMessage("Length of " + id_text + " list must be " +
                                 "either one or " + Utilities::to_string(N)));
        }

      // This should never happen, but return an empty vector so the compiler
      // will be happy
      return std::vector<T> ();
    }

    /**
     * Add standard call for replacing $ASPECT_SOURCE_DIR
     */
    inline
    std::string
    expand_ASPECT_SOURCE_DIR (std::string location)
    {
      return Utilities::replace_in_string(location,
                                          "$ASPECT_SOURCE_DIR",
                                          ASPECT_SOURCE_DIR);
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
        load_file(const std::string &filename,
                  const MPI_Comm &communicator);

        /**
         * Returns the computed data (velocity, temperature, etc. - according
         * to the used plugin) in Cartesian coordinates.
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
        std_cxx11::array<unsigned int,dim-1>
        get_boundary_dimensions (const types::boundary_id boundary_id) const;

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
