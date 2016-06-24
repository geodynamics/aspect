/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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


#ifndef __aspect__geometry_model_chunk_h
#define __aspect__geometry_model_chunk_h

#include <aspect/geometry_model/interface.h>
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
    class Chunk : public Interface<dim>
    {
      public:

        /**
         * Generate a coarse mesh for the geometry described by this class.
         */
        virtual
        void create_coarse_mesh (parallel::distributed::Triangulation<dim> &coarse_grid) const;

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
         * This geometry returns the map <code>{{"inner"->0}, {"outer"->1},
         * {"west"->2}, {"east"->3}}</code> in 2d, and <code>{{"inner"->0},
         * {"outer"->1}, {"west"->2}, {"east"->3}, {"south"->4},
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
         * considers the "outer" boundary as the "surface". In almost
         * all cases one will use a gravity model that also matches
         * these definitions.
         */
        virtual
        double depth(const Point<dim> &position) const;

        virtual
        Point<dim> representative_point(const double depth) const;

        /**
         * Return the longitude at the western edge of the chunk
         * Measured in radians
         */
        virtual
        double west_longitude() const;

        /**
         * Return the longitude at the eastern edge of the chunk
         * Measured in radians
         */
        virtual
        double east_longitude() const;

        /**
         * Returns the longitude range of the chunk
         * Measured in radians
         */
        virtual
        double longitude_range() const;

        /**
         * Return the latitude at the southern edge of the chunk
         * Measured in radians from the equator
         */
        virtual
        double south_latitude() const;

        /**
         * Return the latitude at the northern edge of the chunk
         * Measured in radians from the equator
         */
        virtual
        double north_latitude() const;

        /**
         * Return the latitude range of the chunk
         * Measured in radians
         */
        virtual
        double latitude_range() const;

        /**
         * Return the maximum depth from the surface of the model
         * Measured in meters
         */
        virtual
        double maximal_depth() const;

        /**
         * Return the inner radius of the chunk
         * Measured in meters
         */
        virtual
        double inner_radius() const;

        /**
         * Return the outer radius of the chunk
         * Measured in meters
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
         * in the reference (lat,long,radius) domain and transforms
         * it into real space (cartesian). The inverse function
         * pull_back reverses this operation.
         */

        class ChunkGeometry : public ChartManifold<dim,dim>
        {
          public:
            virtual
            Point<dim>
            pull_back(const Point<dim> &space_point) const;

            virtual
            Point<dim>
            push_forward(const Point<dim> &chart_point) const;
        };

        Point<dim> pull_back(const Point<dim>) const;
        Point<dim> push_forward(const Point<dim>) const;

        /**
         * An object that describes the geometry.
         */
        ChunkGeometry manifold;

        static void set_manifold_ids (Triangulation<dim> &triangulation);
        static void clear_manifold_ids (Triangulation<dim> &triangulation);

    };
  }
}


#endif
