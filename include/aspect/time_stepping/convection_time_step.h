// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2014-2025 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#ifndef _aspect_time_stepping_convection_time_step_h
#define _aspect_time_stepping_convection_time_step_h

#include <aspect/time_stepping/interface.h>


namespace aspect
{
  namespace TimeStepping
  {
    /**
     * Compute the convection time step based on the current solution and
     * return it.
     *
     * @ingroup TimeStepping
     */
    template <int dim>
    class ConvectionTimeStep : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */
        ConvectionTimeStep () = default;


        /**
         * @copydoc aspect::TimeStepping::Interface<dim>::execute()
         */
        double
        execute() override;

      private:
    };
  }
}


#endif
