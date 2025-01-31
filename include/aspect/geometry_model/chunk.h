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
    namespace internal
    {
      /**
       * ChunkGeometry is a class that implements the interface of
       * ChartManifold. The function push_forward takes a point
       * in the reference (radius,lon,lat) domain and transforms
       * it into real space (cartesian). The inverse function
       * pull_back reverses this operation.
       * The push_forward_gradient provides derivatives of the
       * reference coordinates to the real space coordinates,
       * which are used in computing normal vectors.
       * The transformations can include topography added
       * to the initially radially symmetric mesh.
       */
      template <int dim>
      class ChunkGeometry : public ChartManifold<dim,dim>
      {
        public:
          /**
           * Constructor
           */
          ChunkGeometry(const InitialTopographyModel::Interface<dim> &topography,
                        const double min_longitude,
                        const double min_radius,
                        const double max_depth);

          /**
           * Copy constructor
           */
          ChunkGeometry(const ChunkGeometry &other) = default;

          /**
           * This function receives a point in cartesian coordinates x, y and z,
           * including initial prescribed topography and returns
           * radius, longitude, latitude without topography.
           */
          Point<dim>
          pull_back(const Point<dim> &space_point) const override;

          /**
           * This function receives a point in spherical coordinates
           * radius, longitude, latitude and returns cartesian
           * coordinates x, y and z, including any initially prescribed
           * topography.
           */
          Point<dim>
          push_forward(const Point<dim> &chart_point) const override;

          /**
           * This function provides the derivatives of the push_forward
           * function to the spherical coordinates, which are needed
           * in the computation of vectors tangential to the domain boundaries.
           */
          DerivativeForm<1, dim, dim>
          push_forward_gradient(const Point<dim> &chart_point) const override;

          /**
           * This function receives a point in cartesian coordinates x, y and z,
           * and returns radius, longitude, latitude.
           */
          Point<dim>
          pull_back_sphere(const Point<dim> &space_point) const;

          /**
           * This function receives a point in spherical coordinates
           * radius, longitude, latitude and returns cartesian
           * coordinates x, y and z.
           */
          Point<dim>
          push_forward_sphere(const Point<dim> &chart_point) const;

          /**
           * Return the (normalized) normal vector at the point @p p.
           */
          virtual Tensor<1, dim>
          normal_vector(
            const typename Triangulation<dim>::face_iterator &face,
            const Point<dim> &p) const override;

          /**
           * This function computes the outer radius of the domain
           * at the longitude (and latitude) of the given point
           * (given in cartesian coordinates), i.e. the unperturbed
           * outer radius + the topography.
           */
          double
          get_radius(const Point<dim> &space_point) const;

          /**
           * Return a copy of this manifold.
           */
          std::unique_ptr<Manifold<dim,dim>>
          clone() const override;

        private:
          /**
           * A pointer to the topography model.
           */
          const InitialTopographyModel::Interface<dim> *topo;

          /**
           * The minimum longitude of the domain.
           */
          double point1_lon;

          /**
           * The inner radius of the domain.
           */
          double inner_radius;

          /**
           * The maximum depth not taking into account
           * topography (outer radius minus inner radius).
           */
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
      };
    }

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
     * number of cells initialized in each dimension.
     *
     * Initial topography can be added through a radial displacement of the
     * mesh nodes.
     */
    template <int dim>
    class Chunk : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:

        /**
         * Initialization function. This function is called once at the
         * beginning of the program after parse_parameters is run and after
         * the SimulatorAccess (if applicable) is initialized.
         * This function calls the initialize function of the manifold
         * with a pointer to the initial topography model obtained
         * from SimulatorAccess.
         */
        void initialize () override;

        /**
         * Generate a coarse mesh for the geometry described by this class.
         */
        void create_coarse_mesh (parallel::distributed::Triangulation<dim> &coarse_grid) const override;


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
        std::set<types::boundary_id>
        get_used_boundary_indicators () const override;

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
        std::map<std::string,types::boundary_id>
        get_symbolic_boundary_names_map () const override;


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
        double length_scale () const override;

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
        double depth(const Point<dim> &position) const override;

        /**
         * Return the height of the given position relative to the outer
         * radius.
         */
        double height_above_reference_surface(const Point<dim> &position) const override;

        /**
         * @copydoc Interface::representative_point()
         */
        Point<dim> representative_point(const double depth) const override;

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
         * Return the latitude range of the chunk measured in radians
         */
        virtual
        double latitude_range() const;

        /**
         * Return the maximum depth from the surface of the model measured in
         * meters.
         */
        double maximal_depth() const override;

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
        bool
        has_curved_elements() const override;

        /**
         * Return whether the given point lies within the domain specified
         * by the geometry. This function does not take into account
         * initial or dynamic surface topography.
         */
        bool
        point_is_in_domain(const Point<dim> &point) const override;

        /**
         * Returns what the natural coordinate system for this geometry model is,
         * which for a chunk is Spherical.
         */
        aspect::Utilities::Coordinates::CoordinateSystem natural_coordinate_system() const override;

        /**
         * Takes the Cartesian points (x,z or x,y,z) and returns standardized
         * coordinates which are most 'natural' to the geometry model. For a chunk
         * this is (radius, longitude) in 2d and (radius, longitude, latitude) in 3d.
         */
        std::array<double,dim> cartesian_to_natural_coordinates(const Point<dim> &position) const override;

        /**
         * Undoes the action of cartesian_to_natural_coordinates, and turns the
         * coordinate system which is most 'natural' to the geometry model into
         * Cartesian coordinates.
         */
        Point<dim> natural_to_cartesian_coordinates(const std::array<double,dim> &position) const override;

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        void
        parse_parameters (ParameterHandler &prm) override;

      private:
        /**
         * Minimum longitude-depth or
         * longitude-latitude-depth point
         */
        Point<dim> point1;

        /**
         * Maximum longitude-depth or
         * longitude-latitude-depth point
         */
        Point<dim> point2;

        /**
         * The number of cells in each coordinate direction
         */
        std::array<unsigned int, dim> repetitions;

        /**
         * An object that describes the geometry. This pointer is
         * initialized in the initialize() function, and serves as the manifold
         * object that the triangulation is later given in create_coarse_mesh()
         * where the triangulation clones it.
         *
         * The object is marked as 'const' to make it clear that it should not
         * be modified once created. That is because the triangulation copies it,
         * and modifying the current object will not have any impact on the
         * manifold used by the triangulation.
         */
        std::unique_ptr<const internal::ChunkGeometry<dim>> manifold;

        /**
         * Give a symbolic name to the manifold id to be used by this class.
         */
        static constexpr types::manifold_id my_manifold_id = 15;
    };
  }
}


#endif
