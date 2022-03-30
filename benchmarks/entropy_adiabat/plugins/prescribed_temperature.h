/*
  Copyright (C) 2012 - 2022 by the authors of the ASPECT code.

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
    using namespace dealii;

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
         * Return the initial temperature as a function of position.
         */
        double initial_temperature (const Point<dim> &position) const override;
    };
  }
}


#endif
