/*
  Copyright (C) 2014 by the authors of the ASPECT code.

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
#ifndef __aspect__geometry_model_sphere_h
#define __aspect__geometry_model_sphere_h

#include <aspect/geometry_model/interface.h>


namespace aspect
{
  namespace GeometryModel
  {
    using namespace dealii;

    template <int dim>
    class Sphere : public Interface<dim>
    {
      public:
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
         * The sphere model only has one boundary, with indicator zero.
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

        static
        void
        declare_parameters (ParameterHandler &prm);

        virtual
        void
        parse_parameters (ParameterHandler &prm);

        /**
         * Return the radius of the sphere.
         */
        double
        radius () const;

      public:
        /**
         * Radius of the sphere
         */
        double R;

    };
  }
}


#endif
