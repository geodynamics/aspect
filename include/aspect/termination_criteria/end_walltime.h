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


#ifndef _aspect_termination_criteria_end_walltime_h
#define _aspect_termination_criteria_end_walltime_h

#include <aspect/termination_criteria/interface.h>
#include <aspect/simulator.h>

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
    class EndWalltime : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:

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
         * The maximum walltime duration in seconds. The program will be terminated
         * once this value is reached.
         *
         * This variable is read from the parameter file through a parameter
         * called 'Wall time'.
         */
        unsigned int walltime_duration;

        /**
         * The start_walltime is not stored into checkpoint, so it will get reset at restart.
         * Make start_walltime as a static variable that it gets initialized pretty much
         * right away as the program starts, rather than several seconds in once we get
         * to creating the plugin. It gives better estimate of the actual start time.
         */
        static std::time_t start_walltime;
    };
  }
}

#endif
