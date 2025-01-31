/*
  Copyright (C) 2014 - 2023 by the authors of the ASPECT code.

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
#ifndef _aspect_geometry_model_sphere_h
#define _aspect_geometry_model_sphere_h

#include <aspect/geometry_model/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace GeometryModel
  {
    template <int dim>
    class Sphere : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */
        Sphere() = default;

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
         * the radius of the sphere.
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
         * The sphere model only has one boundary, with indicator zero.
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
         * @copydoc Interface::has_curved_elements()
         *
         * Return true because we have a curved boundary.
         */
        bool
        has_curved_elements() const override;

        /**
         * Returns what the natural coordinate system for this geometry model is,
         * which for a sphere is Spherical.
         */
        aspect::Utilities::Coordinates::CoordinateSystem natural_coordinate_system() const override;

        /**
         * Return whether the given point lies within the domain specified
         * by the geometry. This function does not take into account
         * initial or dynamic surface topography.
         */
        bool
        point_is_in_domain(const Point<dim> &point) const override;

        /**
         * Takes the Cartesian points (x,z or x,y,z) and returns standardized
         * coordinates which are most 'natural' to the geometry model. For a sphere
         * this is (radius, longitude) in 2d and (radius, longitude, latitude) in 3d.
         */
        std::array<double,dim> cartesian_to_natural_coordinates(const Point<dim> &position) const override;

        /**
         * Undoes the action of cartesian_to_natural_coordinates, and turns the
         * coordinate system which is most 'natural' to the geometry model into
         * Cartesian coordinates.
         */
        Point<dim> natural_to_cartesian_coordinates(const std::array<double,dim> &position) const override;


        static
        void
        declare_parameters (ParameterHandler &prm);

        void
        parse_parameters (ParameterHandler &prm) override;

        /**
         * Return the radius of the sphere.
         */
        double
        radius () const;

      private:
        /**
         * Radius of the sphere
         */
        double R;
    };
  }
}


#endif
