// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2018-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#ifndef _aspect_time_stepping_repeat_on_cutback_h
#define _aspect_time_stepping_repeat_on_cutback_h

#include <aspect/time_stepping/interface.h>

namespace aspect
{
  namespace TimeStepping
  {
    /**
     * A class that implements a time stepping plugin to repeat a time step if the
     * next time step is significantly smaller than the last step.
     *
     * @ingroup TimeStepping
     */
    template <int dim>
    class RepeatOnCutback : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */
        RepeatOnCutback () = default;

        /**
         * @copydoc aspect::TimeStepping::Interface<dim>::execute()
         */
        double
        execute() override;

        /**
         * The main execute() function.
         */
        std::pair<Reaction, double>
        determine_reaction(const TimeStepInfo &info) override;

        static
        void
        declare_parameters (ParameterHandler &prm);

        void
        parse_parameters (ParameterHandler &prm) override;

      private:
        /**
         * Parameter to determine how much smaller the time step should be
         * repeated as.
         *
         * This variable is read from the parameter file through a parameter
         * called 'Cut back amount'.
         */
        double cut_back_amount;

        /**
         * Parameter that controls when to repeat a time step. If the newly
         * computed step size is smaller than the last step size multiplied by
         * this factor, the step is repeated.
         *
         * This variable is read from the parameter file through a parameter
         * called 'Relative repeat threshold'.
         */
        double repeat_threshold;
    };
  }
}


#endif
