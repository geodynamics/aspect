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


#ifndef _aspect_utilities_h
#define _aspect_utilities_h

#include <aspect/global.h>

#include <array>
#include <deal.II/base/point.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/table_indices.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/component_mask.h>

#include <aspect/coordinate_systems.h>



namespace aspect
{
  template <int dim> class SimulatorAccess;

  namespace GeometryModel
  {
    template <int dim> class Interface;
  }

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
     * This function takes a string argument that is interpreted as a map
     * of the form "key1 : value1, key2 : value2, etc", and then parses
     * it to return a vector of these values where the values are ordered
     * in the same order as a given set of keys. If @p allow_multiple_values_per_key
     * is set to 'true' it also allows entries of the form
     * "key1: value1|value2|value3, etc", in which case the returned
     * vector will have more entries than the provided
     * @p list_of_keys.
     *
     * This function also considers a number of special cases:
     * - If the input string consists of only a comma separated
     *   set of values "value1, value2, value3, ..." (i.e., without
     *   the "keyx :" part), then the input string is interpreted
     *   as if it had had the form "key1 : value1, key2 : value2, ..."
     *   where "key1", "key2", ... are exactly the keys provided by the
     *   @p list_of_keys in the same order as provided. In this situation,
     *   if a background field is required, the background value is
     *   assigned to the first element of the output vector. This form only
     *   allows a single value per key.
     * - Whether or not a background field is required depends on
     *   @p expects_background_field. Requiring a background field
     *   inserts "background" into the @p list_of_keys as the first entry
     *   in the list. Some calling functions (like some material models)
     *   require values for background fields, while others may not
     *   need background values.
     * - Two special keys are recognized:
     *      all --> Assign the associated value to all fields.
     *              Only one key is allowed in this case.
     *      background --> Assign associated value to the background.
     *
     * @param[in] key_value_map The string representation of the map
     *   to be parsed.
     * @param[in] list_of_keys A list of N valid key names that are allowed
     *   to appear in the map. The order of these keys determines the order
     *   of values that are returned by this function.
     * @param[in] expects_background_field If true, expect N+1 values and allow
     *   setting of the background using the key "background".
     * @param[in] property_name A name that identifies the type of property
     *   that is being parsed by this function and that is used in generating
     *   error messages if the map does not conform to the expected format.
     * @param [in] allow_multiple_values_per_key If true allow having multiple values
     *   for each key. If false only allow a single value per key. In either
     *   case each key is only allowed to appear once.
     * @param [in,out] n_values_per_key A pointer to a vector of unsigned
     *   integers. If no pointer is handed over nothing happens. If a pointer
     *   to an empty vector is handed over, the vector
     *   is resized with as many components as
     *   keys (+1 if there is a background field). Each value then stores
     *   how many values were found for this key. The sum over all
     *   entries is the length of the return value of this function.
     *   If a pointer to an existing vector with one or more entries is
     *   handed over (e.g. a n_values_per_key vector created by a
     *   previous call to this function) then this vector is used as
     *   expected structure of the input parameter, and it is checked that
     *   @p key_value_map fulfills this structure. This can be used to assert
     *   that several input parameters prescribe the same number of values
     *   to each key.
     * @param [in] allow_missing_keys Whether to allow that some keys are
     *   not set to any values, i.e. they do not appear at all in the
     *   @p key_value_map. This also allows a completely empty map.
     *
     * @return A vector of values that are parsed from the map, provided
     *   in the order in which the keys appear in the @p list_of_keys argument.
     *   If multiple values per key are allowed, the vector contains first all
     *   values for key 1, then all values for key 2 and so forth. Using the
     *   @p n_values_per_key vector allows the caller to associate entries in the
     *   returned vector with specific keys.
     */
    std::vector<double>
    parse_map_to_double_array (const std::string &key_value_map,
                               const std::vector<std::string> &list_of_keys,
                               const bool expects_background_field,
                               const std::string &property_name,
                               const bool allow_multiple_values_per_key = false,
                               const std::shared_ptr<std::vector<unsigned int> > &n_values_per_key = nullptr,
                               const bool allow_missing_keys = false);

    /**
     * This function takes a string argument that is assumed to represent
     * an input table, in which each row is separated by
     * semicolons, and each column separated by commas. The function
     * returns the parsed entries as a table. In addition this function
     * utilizes the possibly_extend_from_1_to_N() function to accept
     * inputs with only a single column/row, which will be extended
     * to @p n_rows or @p n_columns respectively, to allow abbreviating
     * the input string (e.g. you can provide a single value instead of
     * n_rows by n_columns identical values). This function can for example
     * be used by material models to read densities for different
     * compositions and different phases for each composition.
     */
    template <typename T>
    Table<2,T>
    parse_input_table (const std::string &input_string,
                       const unsigned int n_rows,
                       const unsigned int n_columns,
                       const std::string &property_name);

    /**
     * Given a vector @p var_declarations expand any entries of the form
     * vector(str) or tensor(str) to sublists with component names of the form
     * str_x, str_y, str_z or str_xx, str_xy... for the correct dimension
     * value.
     *
     * This function is to be used for expanding lists of variable names where
     * one or more such variable is actually intended to be a list of
     * components.
     *
     * Returns the generated list of variable names
     */
    template <int dim>
    std::vector<std::string>
    expand_dimensional_variable_names (const std::vector<std::string> &var_declarations);


#if !DEAL_II_VERSION_GTE(9,2,0)
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
#endif

    /**
     * Returns an IndexSet that contains all locally active DoFs that belong to
     * the given component_mask.
     *
     * This function should be moved into deal.II at some point.
     */
    template <int dim>
    IndexSet extract_locally_active_dofs_with_component(const DoFHandler<dim> &dof_handler,
                                                        const ComponentMask &component_mask);



    namespace Coordinates
    {

      /**
       * Returns distance from the Earth's center, latitude and longitude from a
       * given ECEF Cartesian coordinates that account for ellipsoidal shape of
       * the Earth with WGS84 parameters.
       */
      template <int dim>
      std::array<double,dim>
      WGS84_coordinates(const dealii::Point<dim> &position);

      /**
       * Returns spherical coordinates of a Cartesian point. The returned array
       * is filled with radius, phi and theta (polar angle). If the dimension is
       * set to 2 theta is omitted. Phi is always normalized to [0,2*pi].
       *
       */
      template <int dim>
      std::array<double,dim>
      cartesian_to_spherical_coordinates(const dealii::Point<dim> &position);

      /**
       * Return the Cartesian point of a spherical position defined by radius,
       * phi and theta (polar angle). If the dimension is set to 2 theta is
       * omitted.
       */
      template <int dim>
      dealii::Point<dim>
      spherical_to_cartesian_coordinates(const std::array<double,dim> &scoord);

      /**
       * Given a vector defined in the radius, phi and theta directions, return
       * a vector defined in Cartesian coordinates. If the dimension is set to 2
       * theta is omitted. Position is given as a Point in Cartesian coordinates.
       */
      template <int dim>
      Tensor<1,dim>
      spherical_to_cartesian_vector(const Tensor<1,dim> &spherical_vector,
                                    const dealii::Point<dim> &position);


      /**
       * Returns ellipsoidal coordinates of a Cartesian point. The returned array
       * is filled with phi, theta and radius.
       *
       */
      template <int dim>
      std::array<double,3>
      cartesian_to_ellipsoidal_coordinates(const dealii::Point<3> &position,
                                           const double semi_major_axis_a,
                                           const double eccentricity);

      /**
       * Return the Cartesian point of a ellipsoidal position defined by phi,
       * theta and radius.
       */
      template <int dim>
      dealii::Point<3>
      ellipsoidal_to_cartesian_coordinates(const std::array<double,3> &phi_theta_d,
                                           const double semi_major_axis_a,
                                           const double eccentricity);


      /**
       * A function that takes a string representation of the name of a
       * coordinate system (as represented by the CoordinateSystem enum)
       * and returns the corresponding value.
       */
      CoordinateSystem
      string_to_coordinate_system (const std::string &);
    }


    /**
     * Given a 2d point and a list of points which form a polygon, computes if the point
     * falls within the polygon.
     */
    template <int dim>
    bool
    polygon_contains_point(const std::vector<Point<2> > &point_list,
                           const dealii::Point<2> &point);

    /**
     * Given a 2d point and a list of points which form a polygon, compute the smallest
     * distance of the point to the polygon. The sign is negative for points outside of
     * the polygon and positive for points inside the polygon.
     */
    template <int dim>
    double
    signed_distance_to_polygon(const std::vector<Point<2> > &point_list,
                               const dealii::Point<2> &point);


    /**
     * Given a 2d point and a list of two points that define a line, compute the smallest
     * distance of the point to the line segment. When the point's perpendicular
     * base does not lie on the line segment, the smallest distance to the segment's end
     * points is calculated.
     */
    double
    distance_to_line(const std::array<dealii::Point<2>,2 > &point_list,
                     const dealii::Point<2> &point);

    /**
     * Given a vector @p v in @p dim dimensional space, return a set
     * of (dim-1) vectors that are orthogonal to @p v and to each
     * other. The length of each of these vectors equals that of the original
     * vector @p v to ensure that the resulting set of vectors
     * represents a well-conditioned basis.
     */
    template <int dim>
    std::array<Tensor<1,dim>,dim-1>
    orthogonal_vectors (const Tensor<1,dim> &v);

    /**
     * A function for evaluating real spherical harmonics. It takes the degree (l)
     * and the order (m) of the spherical harmonic, where $l \geq 0$ and $0 \leq m \leq l$.
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
     * where $X_{lm}( \theta )$ is an associated Legendre function.
     *
     * In practice it is often convenient to compute the sine ($-l \le m < 0$) and cosine ($0 < m \le l$)
     * variants of the real spherical harmonic at the same time. That is the approach taken
     * here, where we return a pair of numbers, the first corresponding the cosine part and the
     * second corresponding to the sine part. Given this, it is no longer necessary to distinguish
     * between positive and negative $m$, so this function only accepts $ m \ge 0$.
     * For $m = 0$, there is only one part, which is stored in the first entry of the pair.
     *
     * @note This function uses the Boost spherical harmonics implementation internally,
     * which is not designed for very high order (> 100) spherical harmonics computation.
     * If you use spherical harmonics of a high order be sure to confirm the accuracy first.
     * For more information, see:
     * http://www.boost.org/doc/libs/1_49_0/libs/math/doc/sf_and_dist/html/math_toolkit/special/sf_poly/sph_harm.html
     */
    std::pair<double,double> real_spherical_harmonic( unsigned int l, // degree
                                                      unsigned int m, // order
                                                      double theta,   // colatitude (radians)
                                                      double phi );   // longitude (radians)

    /**
     * A struct to enable numerical output with a comma as thousands separator
     */
    struct ThousandSep : std::numpunct<char>
    {
      protected:
        char do_thousands_sep() const override
        {
          return ',';
        }

        std::string do_grouping() const override
        {
          return "\003";  // groups of 3 digits (this string is in octal format)
        }

    };

    /**
     * Checks whether a file named @p filename exists and is readable.
     *
     * @param filename File to check existence
     */
    bool fexists(const std::string &filename);

    /**
     * Checks to see if the user is trying to use data from a url.
     *
     * @param filename File to check
     */
    bool filename_is_url(const std::string &filename);

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
           * @param monotone_spline Whether the cubic spline should be a monotone cubic spline.
           * Requires cubic_spline to be set to true.
           */
          void set_points(const std::vector<double> &x,
                          const std::vector<double> &y,
                          const bool cubic_spline = true,
                          const bool monotone_spline = false);
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
                                 "either one or " + Utilities::to_string(N) +
                                 ". Currently it is " + Utilities::to_string(values.size()) + "."));
        }

      // This should never happen, but return an empty vector so the compiler
      // will be happy
      return std::vector<T> ();
    }

    /**
     * Replace the string <tt>\$ASPECT_SOURCE_DIR</tt> in @p location by the current
     * source directory of ASPECT and return the resulting string.
     */
    std::string
    expand_ASPECT_SOURCE_DIR (const std::string &location);

    /**
     * Given a string @p s, return it in the form ' ("s")' if nonempty.
     * Otherwise just return the empty string itself.
     */
    std::string parenthesize_if_nonempty (const std::string &s);

    /**
     * Returns if a vector of strings @p strings only contains unique
     * entries.
     */
    bool has_unique_entries (const std::vector<std::string> &strings);

    /**
     * AsciiDataLookup reads in files containing input data in ascii format.
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
    class AsciiDataLookup
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
        AsciiDataLookup(const unsigned int components,
                        const double scale_factor);

        /**
         * This constructor relies on the list of column names at the beginning
         * of the model file to determine the number of data components,
         * therefore when using this constructor it is necessary to provide
         * this list in the first uncommented line of the data file.
         */
        explicit AsciiDataLookup(const double scale_factor);

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
         * Returns whether the stored coordinates are equidistant. If
         * coordinates are equidistant the lookup is more efficient. Returns
         * false if no coordinates are loaded at the moment.
         */
        bool
        has_equidistant_coordinates() const;

        /**
         * Returns the coordinates at which data is stored. This function
         * can be used to determine the number of data points, or to query
         * data only at exactly the positions at which it is available (avoiding
         * interpolation).
         *
         * @param dimension The spatial direction for which to return the data
         * coordinates, e.g. 0 for x-direction, 1 for y-direction, or equivalent
         * values if your data coordinates are other dimensions such as
         * temperature, pressure.
         */
        const std::vector<double> &
        get_coordinates(const unsigned int dimension) const;

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
         * The coordinate values in each direction as specified in the data file.
         */
        std::array<std::vector<double>,dim> coordinate_values;

        /**
         * The maximum value of each component
         */
        std::vector<double> maximum_component_value;

        /**
         * The min and max of the coordinates in the data file.
         */
        std::array<std::pair<double,double>,dim> grid_extent;

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
            std::unique_ptr<aspect::Utilities::AsciiDataLookup<dim-1> > > lookups;

        /**
         * Map between the boundary id and the old data objects.
         */
        std::map<types::boundary_id,
            std::unique_ptr<aspect::Utilities::AsciiDataLookup<dim-1> > > old_lookups;

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
        std::unique_ptr<aspect::Utilities::AsciiDataLookup<dim> > lookup;
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
        std::vector<std::unique_ptr<aspect::Utilities::AsciiDataLookup<dim-1> >> lookups;

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
        * Returns the coordinates at which profile data is stored. This function
        * can be used to determine the number of data points, or to query
        * data only at exactly the positions at which it is available (avoiding
        * interpolation).
        */
        const std::vector<double> &
        get_coordinates() const;

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
        std::unique_ptr<aspect::Utilities::AsciiDataLookup<1> > lookup;
    };


    /**
     * This function computes the weighted average $\bar y$ of $y_i$  for a weighted p norm. This
     * leads for a general p to:
     * $\bar y = \left(\frac{\sum_{i=1}^k w_i y_i^p}{\sum_{i=1}^k w_i}\right)^{\frac{1}{p}}$.
     * When p = 0 we take the geometric average:
     * $\bar y = \exp\left(\frac{\sum_{i=1}^k w_i \log\left(y_i\right)}{\sum_{i=1}^k w_i}\right)$,
     * and when $p \le -1000$ or $p \ge 1000$ we take the minimum and maximum norm respectively.
     * This means that the smallest and largest value is respectively taken taken.
     *
     * This function has been set up to be very tolerant to strange values, such as negative weights.
     * The only things we require in for the general p is that the sum of the weights and the sum of
     * the weights times the values to the power p may not be smaller or equal to zero. Furthermore,
     * when a value is zero, the exponent is smaller then one and the correspondent weight is non-zero,
     * this corresponds to no resistance in a parallel system. This means that this 'path' will be followed,
     * and we return zero.
     *
     * The implemented special cases (which are minimum (p <= -1000), harmonic average (p = -1), geometric
     * average (p = 0), arithmetic average (p = 1), quadratic average (RMS) (p = 2), cubic average (p = 3)
     * and maximum (p >= 1000) ) is, except for the harmonic and quadratic averages even more tolerant of
     * negative values, because they only require the sum of weights to be non-zero.
     */
    double weighted_p_norm_average (const std::vector<double> &weights,
                                    const std::vector<double> &values,
                                    const double p);


    /**
     * This function computes the derivative ($\frac{\partial\bar y}{\partial x}$) of an average
     * of the values $y_i(x)$ with regard to $x$, using $\frac{\partial y_i(x)}{\partial x}$.
     * This leads for a general p to:
     * $\frac{\partial\bar y}{\partial x} =
     * \frac{1}{p}\left(\frac{\sum_{i=1}^k w_i y_i^p}{\sum_{i=1}^k w_i}\right)^{\frac{1}{p}-1}
     * \frac{\sum_{i=1}^k w_i p y_i^{p-1} y'_i}{\sum_{i=1}^k w_i}$.
     * When p = 0 we take the geometric average as a reference, which results in:
     * $\frac{\partial\bar y}{\partial x} =
     * \exp\left(\frac{\sum_{i=1}^k w_i \log\left(y_i\right)}{\sum_{i=1}^k w_i}\right)
     * \frac{\sum_{i=1}^k\frac{w_i y'_i}{y_i}}{\sum_{i=1}^k w_i}$
     * and when $p \le -1000$ or $p \ge 1000$ we take the min and max norm respectively.
     * This means that the derivative is taken which has the min/max value.
     *
     * This function has, like the function weighted_p_norm_average been set up to be very tolerant to
     * strange values, such as negative weights. The only things we require in for the general p is that
     * the sum of the weights and the sum of the weights times the values to the power p may not be smaller
     * or equal to zero. Furthermore, when a value is zero, the exponent is smaller then one and the
     * correspondent weight is non-zero, this corresponds to no resistance in a parallel system. This means
     * that this 'path' will be followed, and we return the corresponding derivative.
     *
     * The implemented special cases (which are minimum (p <= -1000), harmonic average (p = -1), geometric
     * average (p = 0), arithmetic average (p = 1), and maximum (p >= 1000) ) is, except for the harmonic
     * average even more tolerant of negative values, because they only require the sum of weights to be non-zero.
     */
    template <typename T>
    T derivative_of_weighted_p_norm_average (const double averaged_parameter,
                                             const std::vector<double> &weights,
                                             const std::vector<double> &values,
                                             const std::vector<T> &derivatives,
                                             const double p);
    /**
     * This function computes a factor which can be used to make sure that the
     * Jacobian remains positive definite.
     *
     * The goal of this function is to find a factor $\alpha$ so that
     * $2\eta(\varepsilon(\mathbf u)) I \otimes I +  \alpha\left[a \otimes b + b \otimes a\right]$ remains a
     * positive definite matrix. Here, $a=\varepsilon(\mathbf u)$ is the @p strain_rate
     * and $b=\frac{\partial\eta(\varepsilon(\mathbf u),p)}{\partial \varepsilon}$ is the derivative of the viscosity
     * with respect to the strain rate and is given by @p dviscosities_dstrain_rate. Since the viscosity $\eta$
     * must be positive, there is always a value of $\alpha$ (possibly small) so that the result is a positive
     * definite matrix. In the best case, we want to choose $\alpha=1$ because that corresponds to the full Newton step,
     * and so the function never returns anything larger than one.
     *
     * The factor is defined by:
     * $\frac{2\eta(\varepsilon(\mathbf u))}{\left[1-\frac{b:a}{\|a\| \|b\|} \right]^2\|a\|\|b\|}$. Alpha is
     * reset to a maximum of one, and if it is smaller then one, a safety_factor scales the alpha to make
     * sure that the 1-alpha won't get to close to zero.
     */
    template <int dim>
    double compute_spd_factor(const double eta,
                              const SymmetricTensor<2,dim> &strain_rate,
                              const SymmetricTensor<2,dim> &dviscosities_dstrain_rate,
                              const double SPD_safety_factor);

    /**
     * Converts an array of size dim to a Point of size dim.
     */
    template <int dim>
    Point<dim> convert_array_to_point(const std::array<double,dim> &array);

    /**
     * Converts a Point of size dim to an array of size dim.
     */
    template <int dim>
    std::array<double,dim> convert_point_to_array(const Point<dim> &point);

    /**
     * A class that represents a binary operator between two doubles. The type of
     * operation is specified on construction time, and can be checked later
     * by using the operator ==. The operator () executes the operation on two
     * double parameters and returns the result. This class is helpful for
     * user specified operations that are not known at compile time.
     */
    class Operator
    {
      public:
        /**
         * An enum of supported operations.
         */
        enum operation
        {
          uninitialized,
          add,
          subtract,
          minimum,
          maximum,
          replace_if_valid
        };

        /**
         * The default constructor creates an invalid operation that will fail
         * if ever executed.
         */
        Operator();

        /**
         * Construct the selected operator.
         */
        Operator(const operation op);

        /**
         * Execute the selected operation with the given parameters and
         * return the result.
         */
        double operator() (const double x, const double y) const;

        /**
         * Return the comparison result between the current operation and
         * the one provided as argument.
         */
        bool operator== (const operation op) const;

      private:
        /**
         * The selected operation of this object.
         */
        operation op;
    };

    /**
     * Create a vector of operator objects out of a list of strings. Each
     * entry in the list must match one of the allowed operations.
     */
    std::vector<Operator> create_model_operator_list(const std::vector<std::string> &operator_names);

    /**
     * Create a string of model operators for use in declare_parameters
     */
    const std::string get_model_operator_options();

    /**
     * A function that returns a SymmetricTensor, whose entries are zero, except for
     * the k'th component, which is set to one. If k is not on the main diagonal the
     * resulting tensor is symmetrized.
     */
    template <int dim>
    SymmetricTensor<2,dim> nth_basis_for_symmetric_tensors (const unsigned int k);

    /**
     * A class that represents a point in a chosen coordinate system.
     */
    template <int dim>
    class NaturalCoordinate
    {
      public:
        /**
         * Constructor based on providing the geometry model as a pointer.
         */
        NaturalCoordinate(Point<dim> &position,
                          const GeometryModel::Interface<dim> &geometry_model);

        /**
         * Constructor based on providing the coordinates and associated
         * coordinate system.
         */
        NaturalCoordinate(const std::array<double, dim> &coord,
                          const Utilities::Coordinates::CoordinateSystem &coord_system);

        /**
         * Returns the coordinates in the given coordinate system, which may
         * not be Cartesian.
         */
        std::array<double,dim> &get_coordinates();

        /**
         * The coordinate that represents the 'surface' directions in the
         * chosen coordinate system.
         */
        std::array<double,dim-1> get_surface_coordinates() const;

        /**
         * The coordinate that represents the 'depth' direction in the chosen
         * coordinate system.
         */
        double get_depth_coordinate() const;

      private:
        /**
         * An enum which stores the the coordinate system of this natural
         * point
         */
        Utilities::Coordinates::CoordinateSystem coordinate_system;

        /**
         * An array which stores the coordinates in the coordinates system
         */
        std::array<double,dim> coordinates;
    };


    /**
     * Compute the cellwise projection of component @p component_index of @p function to the finite element space
     * described by @p dof_handler.
     *
     * @param[in] mapping The mapping object to use.
     * @param[in] dof_handler The DoFHandler the describes the finite element space to
     * project into and that corresponds to @p vec_result.
     * @param[in] component_index The component index of the @p dof_handler for which
     * the projection is being performed. This component should be described by
     * a DG finite element.
     * @param[in] quadrature  The quadrature formula to use to evaluate @p function on each cell.
     * @param[in] function The function to project into the finite element space.
     * This function should store the value of the function at the points described
     * by the 2nd argument into the 3rd argument as an std::vector<double>
     * of size quadrature.size().
     * @param[out] vec_result The output vector where the projected function will be
     * stored in.
     */
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
                     VectorType                                                &vec_result);

    /**
    * Conversion object where one can provide a function that returns
    * a tensor for the velocity at a given point and it returns something
    * that matches the dealii::Function interface with a number of output
    * components equal to the number of components of the finite element
    * in use.
    */
    template <int dim>
    class VectorFunctionFromVelocityFunctionObject : public Function<dim>
    {
      public:
        /**
         * Given a function object that takes a Point and returns a Tensor<1,dim>,
         * convert this into an object that matches the Function@<dim@>
         * interface.
         *
         * @param n_components total number of components of the finite element system.
         * @param function_object The function that will form one component
         *     of the resulting Function object.
         */
        VectorFunctionFromVelocityFunctionObject (const unsigned int n_components,
                                                  const std::function<Tensor<1,dim> (const Point<dim> &)> &function_object);

        /**
         * Return the value of the
         * function at the given
         * point. Returns the value the
         * function given to the constructor
         * produces for this point.
         */
        double value (const Point<dim>   &p,
                      const unsigned int  component = 0) const override;

        /**
         * Return all components of a
         * vector-valued function at a
         * given point.
         *
         * <tt>values</tt> shall have the right
         * size beforehand,
         * i.e. #n_components.
         */
        void vector_value (const Point<dim>   &p,
                           Vector<double>     &values) const override;

      private:
        /**
         * The function object which we call when this class's value() or
         * value_list() functions are called.
         **/
        const std::function<Tensor<1,dim> (const Point<dim> &)> function_object;
    };
  }
}

#endif
