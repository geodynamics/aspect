/*
  Copyright (C) 2011 - 2021 by the authors of the ASPECT code.

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


#ifndef _aspect_geometry_model_initial_topography_model_function_h
#define _aspect_geometry_model_initial_topography_model_function_h

#include <aspect/geometry_model/initial_topography_model/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>

#include <deal.II/base/parsed_function.h>

namespace aspect
{
  namespace InitialTopographyModel
  {
    /**
     * A class that implements initial topography based
     * on a user-defined function..
     *
     * @ingroup InitialTopographyModels
     */
    template <int dim>
    class Function : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */
        Function ();

        /**
         * Return the value of the initial topography as a function
         * of surface position.
         */
        double
        value (const Point<dim-1> &p) const override;

        /**
         * Return the maximum value of the elevation.
         */
        double max_topography () const override;

        static
        void
        declare_parameters (ParameterHandler &prm);

        void
        parse_parameters (ParameterHandler &prm) override;

      private:

        /**
         * The maximum value the topography can take.
         */
        double max_topo;

        /**
         * A function object representing the topography.
         */
        Functions::ParsedFunction<dim> initial_topography_function;

        /**
         * The coordinate representation to evaluate the function. Possible
         * choices are cartesian and spherical.
         */
        Utilities::Coordinates::CoordinateSystem coordinate_system;
    };
  }
}


#endif
