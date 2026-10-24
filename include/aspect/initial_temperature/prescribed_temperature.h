// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2012-2025 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


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
