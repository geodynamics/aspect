/*
  Copyright (C) 2011 - 2019- 2018 by the authors of the ASPECT code.

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
