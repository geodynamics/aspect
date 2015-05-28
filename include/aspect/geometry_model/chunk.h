/*
  Copyright (C) 2011, 2012 by the authors of the ASPECT code.

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

namespace
{
  using namespace dealii;

  template <int dim>
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

}

namespace aspect
{
  namespace GeometryModel
  {
    using namespace dealii;

    /**
     * A class that describes a piece of a circle/sphere delimited by
     * minimum and maximum longitude and radius (2D), or
     * minimum and maximum longitude, latitude and radius (3D)
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
         * We return 1/100th of the diameter of the box.
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
         * Flag whether the box is periodic in longitude. TODO MAKE SURE THIS WORKS!
         */
        bool periodic;

        /**
         * The number of cells in each coordinate direction
         */
        unsigned int repetitions[dim];

        /**
         * An object that describes the geometry.
         */
        ChunkGeometry<dim> manifold;

        static void set_manifold_ids (Triangulation<dim> &triangulation);
        static void clear_manifold_ids (Triangulation<dim> &triangulation);

    };
  }
}


#endif
