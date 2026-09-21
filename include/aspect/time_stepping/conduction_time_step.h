// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2018-2025 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#ifndef _aspect_time_stepping_conduction_time_step_h
#define _aspect_time_stepping_conduction_time_step_h

#include <aspect/time_stepping/interface.h>


namespace aspect
{
  namespace TimeStepping
  {
    /**
     * Compute the conduction time step based on the current solution and
     * return this as the time step.
     *
     * @ingroup TimeStepping
     */
    template <int dim>
    class ConductionTimeStep : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */
        ConductionTimeStep () = default;


        /**
         * @copydoc aspect::TimeStepping::Interface<dim>::execute()
         */
        double
        execute() override;

    };
  }
}


#endif
