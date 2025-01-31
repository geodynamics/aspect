/*
  Copyright (C) 2012 - 2023 by the authors of the ASPECT code.

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


#ifndef _aspect_initial_prescribed_temperature_h
#define _aspect_initial_prescribed_temperature_h

#include <aspect/initial_temperature/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace InitialTemperature
  {
    /**
     * A class sets the initial temperature to the precribed temperature
     * outputs computed by the material model.
     *
     * @ingroup InitialTemperatures
     */
    template <int dim>
    class PrescribedTemperature : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Initialization function. This function is called once at the
         * beginning of the program after parse_parameters is run.
         *
         * This specific function makes sure that the objects that describe
         * initial conditions remain available throughout the run of the
         * program.
         */
        void
        initialize () override;

        /**
         * Return the initial temperature as a function of position.
         */
        double initial_temperature (const Point<dim> &position) const override;

      private:
        /**
         * A shared pointer to the initial composition object
         * that ensures that the current object can continue
         * to access the initial composition object beyond the
         * first time step.
         */
        std::shared_ptr<const aspect::InitialComposition::Manager<dim>> initial_composition;
    };
  }
}


#endif
