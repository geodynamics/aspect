// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2016-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#ifndef _aspect_postprocess_command_h
#define _aspect_postprocess_command_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>

#include <cstdlib>

namespace aspect
{
  namespace Postprocess
  {

    /**
     * A postprocessor that runs a shell command at the end of each time step.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class Command : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Declare new parameters
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Parse parameters
         */
        void
        parse_parameters (ParameterHandler &prm) override;

        /**
         * Execute the command
         */
        std::pair<std::string,std::string>
        execute (TableHandler &statistics) override;

      private:
        /**
         * This variable is read from the parameter file through a parameter called 'Command'.
         */
        std::string command;
        /**
         * This variable is read from the parameter file through a parameter called 'Terminate on failure'.
         */
        bool terminate_on_failure;
        /**
         * This variable is read from the parameter file through a parameter called 'Run on all processes'.
         */
        bool on_all_processes;
    };
  }
}


#endif
