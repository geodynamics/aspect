/*
  Copyright (C) 2018 - 2024 by the authors of the ASPECT code.

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


#ifndef _aspect_time_stepping_repeat_on_nonlinear_fail_h
#define _aspect_time_stepping_repeat_on_nonlinear_fail_h

#include <aspect/time_stepping/interface.h>

namespace aspect
{
  namespace TimeStepping
  {
    /**
     * A class that implements a time stepping plugin to repeat a time step if the
     * nonlinear solver failed to converge in the specified number of iterations.
     * The timestep size will be reduced by the given amount specified by the user
     * and hopefully results in the nonlinear solver converging. If necessary, the
     * timestep size will be repeatedly reduced.
     *
     * This class is automatically enabled if the user specifies that he/she wants
     * to do this.
     *
     * @ingroup TimeStepping
     */
    template <int dim>
    class RepeatOnNonlinearFail : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */
        RepeatOnNonlinearFail ();

        /**
         * @copydoc aspect::TimeStepping::Interface<dim>::execute()
         */
        double
        execute() override;

        /**
         * This function notifies the plugin that the nonlinear solver
         * failed.
         */
        void nonlinear_solver_has_failed() const;

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
         */
        double cut_back_factor;

        /**
         * Enabled by nonlinear_solver_has_failed() to signal that this
         * plugin needs to act in the current timestep;
         */
        mutable bool nonlinear_solver_just_failed;

        /**
         * How many times should we try cutting the timestep size before giving up?
         */
        unsigned int maximum_number_of_repeats;

        /**
         * How many times have we been repeating already in this timestep?
         */
        unsigned int current_number_of_repeats;
    };
  }
}


#endif
