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
/*  $Id: function.h 1509 2012-12-12 05:50:56Z bangerth $  */


#ifndef __aspect__initial_conditions_solidus_h
#define __aspect__initial_conditions_solidus_h

#include <aspect/initial_conditions/interface.h>
#include <aspect/simulator.h>

namespace aspect
{
  namespace InitialConditions
  {

    /**
     * A class that implements temperature initial conditions based on a
     * functional description provided in the input file.
     *
     * @ingroup InitialConditionsModels
     */
    template <int dim>
    class Solidus : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */
        Solidus ();

        /**
         * Return the initial temperature as a function of position.
         */
        virtual
        double initial_temperature (const Point<dim> &position) const;

        /**
         * Declare the parameters this class takes through input files.
         * The default implementation of this function does not describe
         * any parameters. Consequently, derived classes do not have to
         * overload this function if they do not take any runtime parameters.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter
         * file. The default implementation of this function does not read
         * any parameters. Consequently, derived classes do not have to
         * overload this function if they do not take any runtime parameters.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);

      private:
        /**
         * Returns spherical coordinates of a cartesian position.
         */
        const Tensor<1,dim>
        spherical_surface_coordinates(const Tensor<1,dim> &position) const;

		double       litho_thick;
		double       Magnitude_T;
		double       Magnitude_lith;
		double       deltaT;
		int          lateral_wave_number_1;
		int          lateral_wave_number_2;
		std::string  solidus_filename;
    };
  }
  namespace melting
  {
		class Melting_curve
		{
  		public:
   			Melting_curve(const std::string &filename);
			void read(const std::string &filename);
    		double T(const double p, const double radius) const;
			bool is_radius;
			unsigned int Num_points;
  		private:
    		//unsigned int Num_points;
    		std::vector<double> T_array;
    		std::vector<double> P_array;
    		//bool is_radius;
		};
  }
}


#endif
