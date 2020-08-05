/*
  Copyright (C) 2019 - 2020 by the authors of the ASPECT code.

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


#ifndef _aspect_time_stepping_interface_h
#define _aspect_time_stepping_interface_h

#include <aspect/global.h>
#include <aspect/simulator_access.h>
#include <aspect/termination_criteria/interface.h>

namespace aspect
{
  using namespace dealii;

  namespace TimeStepping
  {

    /**
     * A class to handle computation of the next time step (as desired by the user) and
     * checking if the simulation is finished.
     */
    template <int dim>
    class Manager : public SimulatorAccess<dim>
    {
      public:
        /**
         * Override initialize_simulator() so that we can also initialize the contained
         * termination_manager.
         */
        virtual void initialize_simulator (const Simulator<dim> &simulator_object) override;

        /**
         * Compute the size of the next time step potentially taking into account
         * the current solution (convection time step, conduction time step), settings
         * from parameters, and termination criteria (to hit the end time exactly).
         */
        double compute_time_step_size() const;

        /**
         * If true, the simulator should perform a checkpoint before terminating.
         */
        bool need_checkpoint_on_terminate() const;

        /**
         * Check if the simulation is ready to terminate sucessfully.
         */
        bool should_simulation_terminate_now() const;

        /**
         * Declare the parameters of all known termination criteria plugins,
         * as well as of ones this class has itself.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         * This determines which termination criteria objects will be created;
         * then let these objects read their parameters as well.
         */
        void
        parse_parameters (ParameterHandler &prm);

      private:

        /**
         * Whether to do a final checkpoint before termination. This is
         * specified in the parameters.
         */
        bool do_checkpoint_on_terminate;

        /**
         * The termination manager keeps track of the termination plugins and we use
         * it to determine the time_step size in the final time step.
         */
        TerminationCriteria::Manager<dim> termination_manager;
    };

  }
}

#endif
