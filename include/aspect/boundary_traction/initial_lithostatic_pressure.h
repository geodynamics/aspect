/*
  Copyright (C) 2011 - 2023 by the authors of the ASPECT code.

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


#ifndef _aspect_boundary_traction_initial_lithospheric_pressure_h
#define _aspect_boundary_traction_initial_lithospheric_pressure_h

#include <aspect/boundary_traction/interface.h>
#include <aspect/simulator_access.h>
#include <map>


namespace aspect
{
  namespace BoundaryTraction
  {
    /**
     * A class that implements traction boundary conditions by prescribing
     * the lithostatic pressure as the normal traction component.
     *
     * @ingroup BoundaryTractions
     */
    template <int dim>
    class InitialLithostaticPressure : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:

        /**
         * Initialization function. Because this function is called after
         * initializing the SimulatorAccess, all of the necessary information
         * is available to calculate the pressure profile based on the initial
         * temperature and pressure conditions.
         */
        void initialize () override;


        /**
         * Return the boundary traction as a function of position. The
         * (outward) normal vector to the domain is also provided as
         * a second argument.
         */
        Tensor<1,dim>
        boundary_traction (const types::boundary_id boundary_indicator,
                           const Point<dim> &position,
                           const Tensor<1,dim> &normal_vector) const override;


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
         * The number of integration points per profile.
         */
        unsigned int n_points;

        /**
         * Default representative point used for boundaries not listed in
         * representative_point_per_boundary_name. The depth coordinate is ignored.
         */
        Point<dim> default_representative_point;

        /**
         * Per-boundary representative points keyed by symbolic boundary name.
         * Overrides default_representative_point for the listed boundaries.
         */
        std::map<std::string, Point<dim>> representative_point_per_boundary_name;

        /**
         * Computed lithostatic pressure profiles, one entry per boundary indicator.
         */
        std::map<types::boundary_id, std::vector<double>> pressure_profiles;

        /**
         * Depth spacing between integration points, one entry per boundary indicator.
         */
        std::map<types::boundary_id, double> delta_z_per_boundary;

        /**
         * Compute the lithostatic pressure profile along the column defined by
         * the lateral coordinates of representative_point (depth coordinate ignored).
         * Returns {pressure_vector, delta_z}.
         */
        std::pair<std::vector<double>, double>
        compute_pressure_profile (Point<dim> representative_point) const;

        /**
         * Interpolate the stored pressure profile for boundary bid at position p.
         */
        double interpolate_pressure (const Point<dim> &p,
                                     types::boundary_id bid) const;

        /**
         * The id of the bottom boundary.
         */
        types::boundary_id bottom_boundary_id;

        /**
         * Whether or not to prescribe the
         * largest pressure in the lithostatic pressure
         * profile at the bottom boundary independent of
         * actual depth.
         */
        bool prescribe_constant_pressure_at_bottom_boundary;
    };
  }
}


#endif
