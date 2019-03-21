/*
  Copyright (C) 2011 - 2018 by the authors of the ASPECT code.

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


#ifndef _aspect_geometry_model_chunk_h
#define _aspect_geometry_model_chunk_h

#include <aspect/geometry_model/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/compat.h>

#include <deal.II/grid/manifold.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/grid/grid_out.h>

namespace aspect
{
  namespace GeometryModel
  {
    using namespace dealii;

    /**
     * A geometry model class that describes a chunk of a spherical shell.
     * In 2D, this class simulates a sector with inner and outer radius, and
     * minimum and maximum longitude. Longitude increases anticlockwise
     * from the positive x-axis, as per the mathematical convention of phi.
     * In 3D, the class simulates a chunk of a sphere, bounded by arbitrary
     * lines of longitude, latitude and radius. Boundary indicator names are
     * west, east, south, north, inner and outer.
     *
     * The parameters that describe this geometry and that are read from the
     * input file are the inner and outer radii of the shell, the minimum
     * and maximum longitude, minimum and maximum longitude, and the
     * number of cells initialised in each dimension.
     */
    template <int dim>
    class Chunk : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:

        void initialize ();

        /**
         * Generate a coarse mesh for the geometry described by this class.
         */
        virtual
        void create_coarse_mesh (parallel::distributed::Triangulation<dim> &coarse_grid) const;

#if !DEAL_II_VERSION_GTE(9,0,0)
        /**
         * Connect the set/clear manifold id functions to the post/pre_compute_no_normal_flux signals.
         * Or for dealii 9 and above, to set the pointer to the initial topography model again.
         */
        void
        connect_to_signal(SimulatorSignals<dim> &signals);
#endif


        /**
         * Return the set of boundary indicators that are used by this model.
         * This information is used to determine what boundary indicators can
         * be used in the input file.
         *
         * The chunk model uses boundary indicators zero through 2*dim-1, with
         * the first two being the faces perpendicular to the radius of the shell,
         * the next two along lines of longitude and, in 3d, the next two
         * along lines of latitude.
         */
        virtual
        std::set<types::boundary_id>
        get_used_boundary_indicators () const;

        /**
         * Return a mapping from symbolic names of each part of the boundary
         * to the corresponding boundary indicator. This allows users to
         * specify *names*, not just *numbers* in their input files when
         * describing which parts of the boundary have to satisfy which
         * boundary conditions.
         *
         * This geometry returns the map <code>{{"bottom"->0}, {"top"->1},
         * {"west"->2}, {"east"->3}}</code> in 2d, and <code>{{"bottom"->0},
         * {"top"->1}, {"west"->2}, {"east"->3}, {"south"->4},
         * {"north"->5}}</code> in 3d.
         */
        virtual
        std::map<std::string,types::boundary_id>
        get_symbolic_boundary_names_map () const;


        /**
         * Return the typical length scale one would expect of features in
         * this geometry, assuming realistic parameters.
         *
         * As described in the first ASPECT paper, a length scale of
         * 10km = 1e4m works well for the pressure scaling for earth
         * sized spherical shells. use a length scale that
         * yields this value for the R0,R1 corresponding to earth
         * but otherwise scales like (R1-R0)
         */
        virtual
        double length_scale () const;

        /**
         * Return the depth that corresponds to the given
         * position. The documentation of the base class (see
         * GeometryModel::Interface::depth()) describes in detail how
         * "depth" is interpreted in general.
         *
         * Computing a depth requires a geometry model to define a
         * "vertical" direction. The current class considers the
         * radial vector away from the origin as vertical and
         * considers the "top" boundary as the surface. In almost
         * all cases one will use a gravity model that also matches
         * these definitions.
         */
        virtual
        double depth(const Point<dim> &position) const;

        /**
         * Return the height of the given position relative to the outer
         * radius.
         */
        virtual
        double height_above_reference_surface(const Point<dim> &position) const;

       /**
        * Whereas the depth function returns the depth with respect
        * to the unperturbed surface, this function
        * returns the depth with respect to the surface 
        * including the initial topography. For models without
        * initial topography, the result will be the same. 
        *
        * Note that the perturbed surface only considers the 
        * initially prescribed topography, not any perturbations
        * due to a displacement of the free surface.
        */
        virtual
        double depth_wrt_topo(const Point<dim> &position) const;


        virtual
        Point<dim> representative_point(const double depth) const;

        /**
         * Return the longitude at the western edge of the chunk measured in
         * radians.
         */
        virtual
        double west_longitude() const;

        /**
         * Return the longitude at the eastern edge of the chunk measured in
         * radians.
         */
        virtual
        double east_longitude() const;

        /**
         * Return the longitude range of the chunk measured in radians.
         */
        virtual
        double longitude_range() const;

        /**
         * Return the latitude at the southern edge of the chunk measured in
         * radians from the equator.
         */
        virtual
        double south_latitude() const;

        /**
         * Return the latitude at the northern edge of the chunk measured in
         * radians from the equator.
         */
        virtual
        double north_latitude() const;

        /**
         * Return the latitude range of the chunk Measured in radians
         */
        virtual
        double latitude_range() const;

        /**
         * Return the maximum depth from the surface of the model measured in
         * meters.
         */
        virtual
        double maximal_depth() const;

        /**
         * Return the inner radius of the chunk measured in meters.
         */
        virtual
        double inner_radius() const;

        /**
         * Return the outer radius of the chunk measured in meters.
         */
        virtual
        double outer_radius() const;


        /**
         * @copydoc Interface::has_curved_elements()
         *
         * A chunk has curved boundaries and cells, so return true.
         */
        virtual
        bool
        has_curved_elements() const;

        /**
         * Return whether the given point lies within the domain specified
         * by the geometry. This function does not take into account
         * initial or dynamic surface topography.
         */
        virtual
        bool
        point_is_in_domain(const Point<dim> &p) const;

        /*
         * Returns what the natural coordinate system for this geometry model is,
         * which for a chunk is Spherical.
         */
        virtual
        aspect::Utilities::Coordinates::CoordinateSystem natural_coordinate_system() const;

        /**
         * Takes the Cartesian points (x,z or x,y,z) and returns standardized
         * coordinates which are most 'natural' to the geometry model. For a chunk
         * this is (radius, longitude) in 2d and (radius, longitude, latitude) in 3d.
         */
        virtual
        std::array<double,dim> cartesian_to_natural_coordinates(const Point<dim> &position) const;

        /**
         * Undoes the action of cartesian_to_natural_coordinates, and turns the
         * coordinate system which is most 'natural' to the geometry model into
         * Cartesian coordinates.
         */
        virtual
        Point<dim> natural_to_cartesian_coordinates(const std::array<double,dim> &position) const;

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);

      private:
        /**
         * Minimum depth, longitude-depth or
         * longitude-latitude-depth point
         */
        Point<dim> point1;

        /**
         * Maximum depth, longitude-depth or
         * longitude-latitude-depth point
         */
        Point<dim> point2;

        /**
         * The number of cells in each coordinate direction
         */
        unsigned int repetitions[dim];

        /**
         * ChunkGeometry is a class that implements the interface of
         * ChartManifold. The function push_forward takes a point
         * in the reference (radius,lon,lat) domain and transforms
         * it into real space (cartesian). The inverse function
         * pull_back reverses this operation.
         * The push_forward_gradient provides derivatives of the
         * reference coordinates to the real space coordinates,
         * which are used in computing normal vectors.
         * In set_min_longitude the minimum longitude is set,
         * which is used to test the quadrant of returned longitudes
         * in the pull_back function.
         */

        class ChunkGeometry : public ChartManifold<dim,dim>
        {
          public:
            /**
             * Constructor
             */
            ChunkGeometry();

            /**
             * Copy constructor
             */
            ChunkGeometry(const ChunkGeometry &other);

            /*
             * An initialization function necessary to make sure that the
             * manifold has access to the topography plugins.
             */
            void
            initialize(const InitialTopographyModel::Interface<dim> *topography);
            void
            set_topography_pointer(const InitialTopographyModel::Interface<dim> *topography);

            /**
             * This function receives a point in cartesian coordinates x, y and z,
             * including initial prescribed topography and returns
             * radius, longitude, latitude without topography.
             */
            virtual
            Point<dim>
            pull_back(const Point<dim> &space_point) const;

            /**
             * This function receives a point in spherical coordinates
             * radius, longitude, latitude and returns cartesian
             * coordinates x, y and z, including any initially prescribed
             * topography.
             */
            virtual
            Point<dim>
            push_forward(const Point<dim> &chart_point) const;

            /**
             * This function provides the derivatives of the push_forward
             * function to the spherical coordinates, which are needed
             * in the computation of vectors tangential to the domain boundaries.
             */
            virtual
            DerivativeForm<1, dim, dim>
            push_forward_gradient(const Point<dim> &chart_point) const;

            /**
             * This function receives a point in cartesian coordinates x, y and z,
             * and returns radius, longitude, latitude.
             */
            virtual
            Point<dim>
            pull_back_sphere(const Point<dim> &space_point) const;

            /**
             * This function receives a point in spherical coordinates
             * radius, longitude, latitude and returns cartesian
             * coordinates x, y and z.
             */
            virtual
            Point<dim>
            push_forward_sphere(const Point<dim> &chart_point) const;

            /**
             * This function computes the outer radius of the domain
             * at the longitude (and latitude) of the given point
             * (given in cartesian coordinates), i.e. the unperturbed
             * outer radius + the topography.
             */
            virtual
            double
            get_radius(const Point<dim> &space_point) const;

            virtual
            void
            set_min_longitude(const double p1_lon);

#if DEAL_II_VERSION_GTE(9,0,0)
            /**
             * Return a copy of this manifold.
             */
            virtual
            std::unique_ptr<Manifold<dim,dim> >
            clone() const;
#endif

            virtual
            void
            set_min_radius(const double p1_rad);

            virtual
            void
            set_max_depth(const double p2_rad_minus_p1_rad);


          private:
            // The minimum longitude of the domain
            double point1_lon;
            double inner_radius;

            double max_depth;

            /**
             * This function removes the initial topography from a
             * given point in spherical coordinates R+topo, lon, lat.
             * I.e. it returns R, lon, lat.
             */
            virtual
            Point<dim>
            pull_back_topo(const Point<dim> &space_point) const;

            /**
             * This function adds the initial topography to a
             * given point in spherical coordinates R, lon, lat.
             * I.e. it returns R+topo, lon, lat.
             */
            virtual
            Point<dim>
            push_forward_topo(const Point<dim> &chart_point) const;

            const InitialTopographyModel::Interface<dim> *topo;
        };

        /**
         * An object that describes the geometry.
         */
        ChunkGeometry manifold;

#if !DEAL_II_VERSION_GTE(9,0,0)
        void set_manifold_ids (typename parallel::distributed::Triangulation<dim> &triangulation);
        void clear_manifold_ids (typename parallel::distributed::Triangulation<dim> &triangulation);
#endif
    };
  }
}


#endif
