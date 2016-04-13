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

#ifndef __aspect__termination_criteria_end_step_h
#define __aspect__termination_criteria_end_step_h

#include <aspect/termination_criteria/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace TerminationCriteria
  {

    /**
     * A class that terminates the simulation when a specified end timestep is
     * reached.
     *
     * @ingroup TerminationCriteria
     */
    template <int dim>
    class EndStep : public Interface<dim>, public SimulatorAccess<dim>
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
        unsigned int end_step;
    };
  }
}

#endif
