/*
  Copyright (C) 2011 - 2021 by the authors of the ASPECT code.

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


#ifndef _aspect_geometry_model_two_merged_chunks_h
#define _aspect_geometry_model_two_merged_chunks_h

#include <aspect/geometry_model/chunk.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/compat.h>

#include <deal.II/grid/manifold.h>
#include <deal.II/base/function_lib.h>

namespace aspect
{
  namespace GeometryModel
  {
    using namespace dealii;


    /**
     * A geometry model class that describes a chunk of a spherical shell,
     * but with two boundary indicators per side boundary. This allows
     * for specifying two different types of boundary conditions on the
     * side boundaries, for example prescribed plate motions along
     * the edge of the lithosphere and an open boundary underneath.
     * In 2D, this class simulates a sector with inner and outer radius, and
     * minimum and maximum longitude. Longitude increases anticlockwise
     * from the positive x-axis, as per the mathematical convention of phi.
     * In 3D, the class simulates a chunk of a sphere, bounded by arbitrary
     * lines of longitude, latitude and radius. Boundary indicator names are
     * lowerwest, upperwest, lowereast, uppereast, lowersouth, uppersouth,
     * lowernorth, uppernorth, bottom and top.
     *
     * The parameters that describe this geometry and that are read from the
     * input file are the inner, outer and middle boundary radius of the shell,
     * the minimum and maximum longitude, minimum and maximum longitude, and the
     * number of cells initialised in each dimension. In the radial direction,
     * this number should be specified for both the lower and the upper part of
     * the chunk (with their boundary lying at the middle boundary radius).
     *
     * Initial topography can be added through a radial displacement of the
     * mesh nodes.
     */
    template <int dim>
    class TwoMergedChunks : public Interface<dim>, public SimulatorAccess<dim>
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
         * This function calls the initialize function of the manifold
         * with the given pointer to the initial topography model.
         */
        void set_topography_model (const InitialTopographyModel::Interface<dim> *topo_pointer);

        /**
         * Generate a coarse mesh for the geometry described by this class.
         */
        void create_coarse_mesh (parallel::distributed::Triangulation<dim> &coarse_grid) const override;

        /**
         * Return the set of boundary indicators that are used by this model.
         * This information is used to determine what boundary indicators can
         * be used in the input file.
         *
         * The chunk model uses boundary indicators zero through 2*dim+2*(dim-1)-1, with
         * the first two being the faces perpendicular to the radius of the shell,
         * the next four along lines of longitude and, in 3d, the next four
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
         * Whereas the depth function returns the depth with respect
         * to the unperturbed surface, this function
         * returns the depth with respect to the surface
         * including the initial topography. For models without
         * initial topography, the result will be the same.
         *
         * Note that the perturbed surface only considers the
         * initially prescribed topography, not any perturbations
         * due to a displacement of the free surface. Therefore,
         * be careful with using this function if the surface changes
         * over time.
         */
        double depth_wrt_topo(const Point<dim> &position) const;

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
         * Minimum longitude-depth (2D) or
         * longitude-latitude-depth (3D) point
         * of the entire merged chunk.
         */
        Point<dim> point1;

        /**
         * Maximum longitude-depth (2D) or
         * longitude-latitude-depth (3D) point
         * of the entire merged chunk.
         */
        Point<dim> point2;

        /**
         * Minimum longitude-depth (2D) or
         * longitude-latitude-depth (3D) point
         * for the upper chunk.
         */
        Point<dim> point3;

        /**
         * Maximum longitude-depth (2D) or
         * longitude-latitude-depth (3D) point
         * for the lower chunk.
         */
        Point<dim> point4;

        /**
         * The number of cells in each coordinate direction
         * for the lower and upper chunk.
         */
        std::vector<unsigned int> lower_repetitions;
        std::vector<unsigned int> upper_repetitions;

        /**
         * An object that describes the geometry.
         */
        internal::ChunkGeometry<dim> manifold;

        /**
         * Bind boundary indicators to child cells after each mesh refinement round.
         */
        virtual void set_boundary_indicators (parallel::distributed::Triangulation<dim> &triangulation) const;
    };
  }
}


#endif
