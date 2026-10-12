// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2019-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#ifndef _aspect_termination_criteria_steady_temperature_h
#define _aspect_termination_criteria_steady_temperature_h

#include <aspect/termination_criteria/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace TerminationCriteria
  {
    /**
     * A class that implements a termination criterion based on the steady state
     * of the average temperature.
     *
     * @ingroup TerminationCriteria
     */
    template <int dim>
    class SteadyTemperature : public Interface<dim>, public SimulatorAccess<dim>
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

        void
        save (std::map<std::string, std::string> &status_strings) const override;

        void
        load (const std::map<std::string, std::string> &status_strings) override;

      private:
        /**
         * This variable is read from the parameter file through a parameter
         * called 'Time in steady state'.
         */
        double                                  necessary_time_in_steady_state;

        /**
         * This variable is read from the parameter file through a parameter
         * called 'Maximum relative deviation'.
         */
        double                                  allowed_relative_deviation;

        /**
         * A list of pairs (time, temperature) that we have computed at
         * previous time steps. This is used to determine when we have reached
         * steady state.
         */
        std::list<std::pair<double, double>>   time_temperature;
    };
  }
}

#endif
