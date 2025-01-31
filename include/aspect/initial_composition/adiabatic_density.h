/*
  Copyright (C) 2017 - 2022 by the authors of the ASPECT code.

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


#ifndef _aspect_initial_composition_adiabatic_density_h
#define _aspect_initial_composition_adiabatic_density_h

#include <aspect/initial_composition/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>

namespace aspect
{
  namespace InitialComposition
  {
    /**
     * A class that implements initial conditions for the compositional fields
     * based on the adiabatic density profile. Note that only the field
     * of the type 'density' will be filled, for all other fields
     * this plugin returns 0.0.
     *
     * @ingroup InitialCompositionModels
     */
    template <int dim>
    class AdiabaticDensity : public Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        /**
         * Return the initial composition as the adiabatic density at this
         * point.
         */
        double initial_composition (const Point<dim> &position, const unsigned int n_comp) const override;
    };
  }
}


#endif
