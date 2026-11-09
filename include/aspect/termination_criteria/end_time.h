// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#ifndef _aspect_termination_criteria_end_time_h
#define _aspect_termination_criteria_end_time_h

#include <aspect/termination_criteria/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace TerminationCriteria
  {

    /**
     * A class that terminates the simulation when a specified end time is
     * reached.
     *
     * @ingroup TerminationCriteria
     */
    template <int dim>
    class EndTime : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:

        /**
         * Check this termination criterion and possibly reduce time step size
         */
        double check_for_last_time_step (const double time_step) const override;

        /**
         * Evaluate this termination criterion.
         *
         * @return Whether to terminate the simulation (true) or continue
         * (false).
         */
        bool
        execute () override;

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
         * This variable is read from the parameter file through a parameter
         * called 'End time'.
         */
        double end_time;
    };
  }
}

#endif
