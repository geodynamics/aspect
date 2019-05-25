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


#ifndef _aspect_geometry_model_spherical_shell_h
#define _aspect_geometry_model_spherical_shell_h

#include <aspect/geometry_model/interface.h>
#include <aspect/simulator_access.h>

#include <deal.II/grid/manifold_lib.h>

namespace aspect
{
  namespace GeometryModel
  {
    using namespace dealii;

    /**
     * A class that describes the geometry as a spherical shell. To be more
     * precise, at least in 2d this class can also simulate just a sector of
     * the spherical shell geometry, in particular a half ring and a quarter
     * ring.
     *
     * The parameters that describe this geometry and that are read from the
     * input file are the inner and outer radii of the shell and the opening
     * angle of the section of the shell we want to build.
     */
    template <int dim>
    class CustomSphericalShell : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */
        CustomSphericalShell();

        /**
         * Generate a coarse mesh for the geometry described by this class.
         */
        virtual
        void create_coarse_mesh (parallel::distributed::Triangulation<dim> &coarse_grid) const;

        /**
         * Return the typical length scale one would expect of features in
         * this geometry, assuming realistic parameters.
         *
         * As discussed in the step-32 tutorial program, an appropriate length
         * scale for this geometry is 10km, so we return $10^4$ here.
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

        /**
         * Return the height of the given position relative to
         * the outer radius of the shell.
         */
        virtual
        double height_above_reference_surface(const Point<dim> &position) const;

        virtual
        Point<dim> representative_point(const double depth) const;

        virtual
        double maximal_depth() const;

        /**
         * Return the set of boundary indicators that are used by this model.
         * This information is used to determine what boundary indicators can
         * be used in the input file.
         *
         * The spherical shell model uses boundary indicators zero and one,
         * with zero corresponding to the inner surface and one corresponding
         * to the outer surface. In 2d, if the geometry is only a slice of the
         * shell, boundary indicators 2 and 3 indicate the left and right
         * radial bounding lines.
         */
        virtual
        std::set<types::boundary_id>
        get_used_boundary_indicators () const;

        /**
         * Return symbolic names for all boundary components. Their names are
         * described in the documentation of this plugin, at the bottom of the
         * .cc file.
         */
        virtual
        std::map<std::string,types::boundary_id>
        get_symbolic_boundary_names_map () const;

        /**
         * @copydoc Interface::has_curved_elements()
         *
         * Return true because we have a curved boundary.
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
        point_is_in_domain(const Point<dim> &point) const;

        /*
         * Returns what the natural coordinate system for this geometry model is,
         * which for a spherical shell is Spherical.
         */
        virtual
        aspect::Utilities::Coordinates::CoordinateSystem natural_coordinate_system() const;

        /**
         * Declare the parameters this class takes through input files. The
         * default implementation of this function does not describe any
         * parameters. Consequently, derived classes do not have to overload
         * this function if they do not take any runtime parameters.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         * The default implementation of this function does not read any
         * parameters. Consequently, derived classes do not have to overload
         * this function if they do not take any runtime parameters.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);

        /**
         * Return the inner radius of the shell.
         */
        double
        inner_radius () const;

        /**
         * Return the outer radius of the shell.
         */
        double
        outer_radius () const;

      private:
        /**
         * Inner and outer radii of the spherical shell.
         */
        double R0, R1;

        /**
         * Number of slices for extrusion in case of tailored spherical shell.
         */
        unsigned int n_slices;

        /**
         * Number of tangential mesh cells in the initial, coarse mesh.
         */
        int n_cells_along_circumference;

        /**
         * Height slices for extrusion in case of tailored spherical shell.
         */
        double height;

        /**
         * The manifold that describes the geometry.
         */
        const SphericalManifold<dim> spherical_manifold;

        /**
         * Set the manifold ids on all cells (also boundaries) before
         * refinement to generate well shaped cells.
         */
        void set_manifold_ids (parallel::distributed::Triangulation<dim> &triangulation) const;
    };
  }
}


#endif
