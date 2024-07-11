/*
  Copyright (C) 2018 - 2020 by the authors of the ASPECT code.

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


#ifndef _aspect_time_stepping_repeat_on_cutback_h
#define _aspect_time_stepping_repeat_on_cutback_h

#include <aspect/time_stepping/interface.h>

namespace aspect
{
  namespace TimeStepping
  {
    using namespace dealii;

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
         */
        double cut_back_amount;

        /**
         * Parameter that controls when to repeat a time step. If the newly
         * computed step size is smaller than the last step size multiplied by
         * this factor, the step is repeated.
         */
        double repeat_threshold;
    };
  }
}


#endif
