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
         * The box model uses boundary indicators zero through 2*dim-1, with
         * the first two being the faces perpendicular to the x-axis, the next
         * two perpendicular to the y-axis, etc.
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
         * This geometry returns the map <code>{{"left"->0}, {"right"->1},
         * {"bottom"->2}, {"top"->3}}</code> in 2d, and <code>{{"left"->0},
         * {"right"->1}, {"front"->2}, {"back"->3}, {"bottom"->4},
         * {"top"->5}}</code> in 3d.
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


        virtual
        double depth(const Point<dim> &position) const;

        virtual
        Point<dim> representative_point(const double depth) const;

        virtual
        double start_longitude() const;

        virtual
        double end_longitude() const;

        virtual
        double longitude_range() const;

        virtual
        double maximal_depth() const;

        virtual
        double inner_radius() const;

        virtual
        double outer_radius() const;


        /**
         * @copydoc Interface::has_curved_elements()
         *
         * A box has only straight boundaries and cells, so return false.
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
