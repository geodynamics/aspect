/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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


#ifndef __aspect__prescribed_stokes_solution_function_h
#define __aspect__prescribed_stokes_solution_function_h

#include <aspect/prescribed_stokes_solution/interface.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/parsed_function.h>

namespace aspect
{
  namespace PrescribedStokesSolution
  {
    using namespace dealii;

    /**
     * A class that implements velocity and pressure solutions based on a
     * functional description provided in the input file.
     *
     * @ingroup PrescribedStokesSolution
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
         * Return the velocity and pressure as a function of position.
         */
        virtual
        void
        stokes_solution (const Point<dim> &p, Vector<double> &value) const;

        /**
         * A function that is called at the beginning of each time step to
         * indicate what the model time is for which the velocity and
         * pressure values will next be evaluated. For the current class,
         * the function passes to the parsed function what the current time is.
         */
        virtual
        void
        update ();

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
         * A function object representing the components of the velocity.
         */
        Functions::ParsedFunction<dim> prescribed_velocity_function;
        /**
         * A function object representing the pressure.
         */
        Functions::ParsedFunction<dim> prescribed_pressure_function;
    };
  }
}


#endif
