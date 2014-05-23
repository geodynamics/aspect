/*
  Copyright (C) 2012 by the authors of the ASPECT code.

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
/*  $Id: two_plates.h - Nov 19, 2013 - Katrina M Arredondo - Modified from boussinesq_adiab.h $  */


#ifndef __aspect__initial_conditions_two_plates_h
#define __aspect__initial_conditions_two_plates_h

#include <aspect/initial_conditions/interface.h>
#include <aspect/simulator.h>

#include <deal.II/base/parsed_function.h>

namespace aspect
{
  namespace InitialConditions
  {
    using namespace dealii;

    /**
     * A class that implements adiabatic initial conditions
     * for the temperature field and, optional, upper and
     * lower thermal boundary layers calculated
     * using the half-space cooling model. The age of the
     * boundary layers are input parameters.
     *
     * @ingroup InitialConditionsModels
     */
    template <int dim>
    class Two_plates : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Return the initial temperature as a function of position.
         */
        virtual
        double initial_temperature (const Point<dim> &position) const;

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter
         * file.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);

      private:
        /**
         * Age of the thermal boundary layer at the
         * surface of the model. If set to zero,
         * no boundary layer will be present in the model.
         */
	double ovplate_age_ma;
        double plate_age_ma;
	double trench_distance;
	double T_surface;
        double T1;
        double K_0;
	double wedge_angle;
	double wedge_depth;	
    };
  }
}


#endif
