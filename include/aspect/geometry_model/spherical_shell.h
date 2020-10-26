/*
  Copyright (C) 2011 - 2019 by the authors of the ASPECT code.

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
    class SphericalShell : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */
        SphericalShell();

        /**
         * Generate a coarse mesh for the geometry described by this class.
         */
        void create_coarse_mesh (parallel::distributed::Triangulation<dim> &coarse_grid) const override;

        /**
         * Return the typical length scale one would expect of features in
         * this geometry, assuming realistic parameters.
         *
         * As discussed in the step-32 tutorial program, an appropriate length
         * scale for this geometry is 10km, so we return $10^4$ here.
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
         * considers the "outer" boundary as the "surface". In almost
         * all cases one will use a gravity model that also matches
         * these definitions.
         */
        double depth(const Point<dim> &position) const override;

        /**
         * Return the height of the given position relative to
         * the outer radius of the shell.
         */
        double height_above_reference_surface(const Point<dim> &position) const override;

        /**
         * @copydoc Interface<dim>::representative_point()
         */
        Point<dim> representative_point(const double depth) const override;

        /**
         * @copydoc Interface<dim>::maximal_depth()
         */
        double maximal_depth() const override;

        /**
         * Return the set of boundary indicators that are used by this model.
         * This information is used to determine what boundary indicators can
         * be used in the input file.
         *
         * The spherical shell may be generated as per in the original code
         * (with respect to the inner and outer radius, and an initial number
         * of cells along circumference) or following a custom mesh scheme:
         * list of radial values or number of slices. A surface mesh is first
         * generated and refined as desired, before it is extruded radially.
         * A list of radial values subdivides the spherical shell at specified
         * radii. The number of slices subdivides the spherical shell into N
         * slices of equal thickness. The custom spherical shell only works
         * with an opening angle of 360 degrees.
         *
         * The spherical shell model uses boundary indicators zero and one,
         * with zero corresponding to the inner surface and one corresponding
         * to the outer surface. In 2d, if the geometry is only a slice of the
         * shell, boundary indicators 2 and 3 indicate the left and right
         * radial bounding lines.
         */
        std::set<types::boundary_id>
        get_used_boundary_indicators () const override;

        /**
         * Return symbolic names for all boundary components. Their names are
         * described in the documentation of this plugin, at the bottom of the
         * .cc file.
         */
        std::map<std::string,types::boundary_id>
        get_symbolic_boundary_names_map () const override;

        /**
         * Return the set of periodic boundaries as described in the input
         * file.
         */
        std::set< std::pair< std::pair<types::boundary_id, types::boundary_id>, unsigned int> >
        get_periodic_boundary_pairs () const override;

        /**
         * @copydoc Interface::has_curved_elements()
         *
         * Return true because we have a curved boundary.
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
         * which for a spherical shell is Spherical.
         */
        aspect::Utilities::Coordinates::CoordinateSystem natural_coordinate_system() const override;

        /**
         * Takes the Cartesian points (x,z or x,y,z) and returns standardized
         * coordinates which are most 'natural' to the geometry model. For a spherical shell
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
        void
        parse_parameters (ParameterHandler &prm) override;

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

        /**
         * Return the opening angle of the shell sector.
         */
        double
        opening_angle () const;

        /**
         * Collects periodic boundaries constraints for the given geometry,
         * which will be added to the existing @p constraints.
         */
        void
        make_periodicity_constraints(const DoFHandler<dim> &dof_handler,
                                     AffineConstraints<double> &constraints) const override;

      private:
        /**
         * Specify the radial subdivision of the spherical shell
         * mesh.
         */
        enum CustomMeshRadialSubdivision
        {
          none,
          list,
          slices
        } custom_mesh;

        /**
         * Initial surface refinement for the custom mesh cases.
         */
        int initial_lateral_refinement;

        /**
         * Initial surface refinement for the custom mesh cases.
         */
        unsigned int n_slices;

        /**
         * List of radial values for the list custom mesh.
         */
        std::vector<double> R_values_list;

        /**
         * Inner and outer radii of the spherical shell.
         */
        double R0, R1;

        /**
         * Opening angle of the section of the shell that we simulate.
         */
        double phi;

        /**
         * Number of tangential mesh cells in the initial, coarse mesh.
         */
        int n_cells_along_circumference;

        /**
         * The manifold that describes the geometry.
         */
        const SphericalManifold<dim> spherical_manifold;

        /**
         * Set the manifold ids on all cells (also boundaries) before
         * refinement to generate well shaped cells.
         */
        void set_manifold_ids (parallel::distributed::Triangulation<dim> &triangulation) const;

        /**
         * Flag whether the 2D quarter shell is periodic in phi.
         */
        bool periodic;
    };
  }
}


#endif
