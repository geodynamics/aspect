/*
  Copyright (C) 2011 - 2016 by the authors of the ASPECT code.

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


#ifndef _aspect_geometry_model_box_h
#define _aspect_geometry_model_box_h

#include <aspect/geometry_model/interface.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace GeometryModel
  {
    using namespace dealii;

    /**
     * A class that describes a box geometry of certain width, height, and
     * depth (in 3d), and, possibly, topography.
     */
    template <int dim>
    class Box : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Initialization function. This function is called once at the
         * beginning of the program after parse_parameters is run and after
         * the SimulatorAccess (if applicable) is initialized.
         */
        void initialize ();

        /**
         * Add initial topography to the mesh.
         */
        void topography (typename parallel::distributed::Triangulation<dim> &grid) const;

        /**
         * Relocate the vertical coordinate of the given point based on
         * the topography at the surface specified by the initial topography
         * model.
         */
        Point<dim> add_topography (const Point<dim> &x_y_z) const;

        /**
         * Generate a coarse mesh for the geometry described by this class.
         */
        virtual
        void create_coarse_mesh (parallel::distributed::Triangulation<dim> &coarse_grid) const;

        /**
         * Return a point that denotes the size of the box in each dimension
         * of the domain.
         */
        virtual
        Point<dim> get_extents () const;

        /**
         * Return a point that denotes the lower left corner of the box
         * domain.
         */
        virtual
        Point<dim> get_origin () const;

        /**
         * Return the typical length scale one would expect of features in
         * this geometry, assuming realistic parameters.
         *
         * We return 1/100th of the diameter of the box.
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
         * $(0,1)^T$ vector in 2d (and the $(0,0,1)^T$ vector in 3d)
         * as vertical and considers the "top" boundary as the
         * surface. In almost all cases one will use a gravity model
         * that also matches these definitions.
         *
         * Note that the depth is calculated with respect to the
         * surface without initial topography.
         */
        virtual
        double depth(const Point<dim> &position) const;

        virtual
        Point<dim> representative_point(const double depth) const;

        virtual
        double maximal_depth() const;

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
         * Return the set of periodic boundaries as described in the input
         * file.
         */
        virtual
        std::set< std::pair< std::pair<types::boundary_id, types::boundary_id>, unsigned int> >
        get_periodic_boundary_pairs () const;

        /**
         * @copydoc Interface::has_curved_elements()
         *
         * A box has only straight boundaries and cells, so return false.
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
         * Extent of the box in x-, y-, and z-direction (in 3d).
         */
        Point<dim> extents;

        /**
         * Origin of the box in x, y, and z (in 3d) coordinates.
         */
        Point<dim> box_origin;

        /**
         * Flag whether the box is periodic in the x-, y-, and z-direction.
         */
        bool periodic[dim];

        /**
         * The number of cells in each coordinate direction.
         */
        unsigned int repetitions[dim];

        /**
         * A pointer to the initial topography model.
         */
        InitialTopographyModel::Interface<dim> *topo_model;

    };
  }
}


#endif
