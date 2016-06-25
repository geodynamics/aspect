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


#ifndef __aspect__traction_boundary_conditions_lithospheric_pressure_h
#define __aspect__traction_boundary_conditions_lithospheric_pressure_h

#include <aspect/traction_boundary_conditions/interface.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace TractionBoundaryConditions
  {
    using namespace dealii;

    /**
     * A class that implements traction boundary conditions based on a
     * functional description provided in the input file.
     *
     * @ingroup TractionBoundaryConditionsModels
     */
    template <int dim>
    class LP : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:

        /**
         * Initialization function. Because this function is called after
         * initializing the SimulatorAccess, all of the necessary information
         * is available to calculate the pressure profile.
         */
        virtual void initialize ();


        /**
         * Return the boundary traction as a function of position. The
         * (outward) normal vector to the domain is also provided as
         * a second argument.
         */
        virtual
        Tensor<1,dim>
        traction (const Point<dim> &position,
                  const Tensor<1,dim> &normal_vector) const;


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

        std::set<types::boundary_id> traction_bi;

        /**
         * The number of integration points.
         */
        unsigned int n_points;

        /**
         * Interval spacing between each two data points in the tables above
         * with regard to the depth coordinate.
         */
        double delta_z;

        /*
         * The user-specified point where to calculate the pressure profile
         */
        Point<dim> representative_point;

        /**
         * The lithostatic pressure profile.
         */
        std::vector<double> pressure;

        /**
         * Return the lithostatic pressure at a given point of the domain.
         */
        double get_pressure (const Point<dim> &p) const;

        /**
         * Convert spherical coordinates to cartesian coordinates
         * TODO: use Utilities spherical_coordinates function?
         */
        Point<dim> spherical_to_cart(const Point<dim >spherical_coord) const;


    };
  }
}


#endif
