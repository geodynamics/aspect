/*
  Copyright (C) 2011, 2012 by the authors of the ASPECT code.

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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/


#ifndef __aspect__termination_criteria_steady_rms_velocity_h
#define __aspect__termination_criteria_steady_rms_velocity_h

#include <aspect/termination_criteria/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace TerminationCriteria
  {

    /**
     * A class that implements a termination criterion based on steady state
     * of the root mean square of the velocity field.
     *
     * @ingroup TerminationCriteria
     */
    template <int dim>
    class SteadyRMSVelocity : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Evaluate this termination criterion.
         *
         * @return Whether to terminate the simulation (true) or continue
         * (false).
         */
        virtual
        bool
        execute (void);

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);

      private:
        double                                  time_length;
        double                                  relative_deviation;

        /**
         * A list of pairs (time, rms_velocity) that we have computed at
         * previous time steps. This is used to determine when we have reached
         * steady state.
         */
        std::list<std::pair<double, double> >   time_rmsvel;

    };
  }
}

#endif
