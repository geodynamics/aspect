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


#ifndef _aspect_prescribed_stokes_solution_function_h
#define _aspect_prescribed_stokes_solution_function_h

#include <aspect/prescribed_stokes_solution/interface.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/parsed_function.h>

namespace aspect
{
  namespace PrescribedStokesSolution
  {
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
        void
        stokes_solution (const Point<dim> &position, Vector<double> &value) const override;

        /**
         * A function that is called at the beginning of each time step to
         * indicate what the model time is for which the velocity and
         * pressure values will next be evaluated. For the current class,
         * the function passes to the parsed function what the current time is.
         */
        void
        update () override;

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
         * A function object representing the components of the velocity.
         */
        Functions::ParsedFunction<dim> prescribed_velocity_function;
        /**
         * A function object representing the pressure.
         */
        Functions::ParsedFunction<dim> prescribed_pressure_function;
        /**
         * A function object representing the fluid pressure (in models with melt transport).
         */
        Functions::ParsedFunction<dim> prescribed_fluid_pressure_function;
        /**
         * A function object representing the compaction pressure (in models with melt transport).
         */
        Functions::ParsedFunction<dim> prescribed_compaction_pressure_function;
        /**
         * A function object representing the components of the fluid velocity (in models with melt transport).
         */
        Functions::ParsedFunction<dim> prescribed_fluid_velocity_function;
    };
  }
}


#endif
