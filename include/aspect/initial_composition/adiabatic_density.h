// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2017-2025 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


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
