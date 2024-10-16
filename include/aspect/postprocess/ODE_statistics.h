/*
  Copyright (C) 2024 by the authors of the ASPECT code.

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


#ifndef _aspect_postprocess_ODE_statistics_h
#define _aspect_postprocess_ODE_statistics_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace Postprocess
  {

    /**
     * A postprocessor that computes statistics about the ODEs solved in ASPECT.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class ODEStatistics : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Connect the callback functions to the respective signals.
         */
        void initialize() override;

        /**
         * Write the ODE statistics into the statistics file.
         */
        std::pair<std::string,std::string>
        execute (TableHandler &statistics) override;

      private:
        /**
         * This function clears all the saved solver information that is
         * stored in this class. It is executed every time the output is
         * written into the statistics object.
         */
        void
        clear_data();

        /**
         * Callback function that is connected to the post_ARKode_solve
         * signal to store the solver history.
         */
        void
        store_ODE_solver_history(const unsigned int iteration_count);

        /**
         * Variables that store the ODE solver history of the current
         * timestep, until they are written into the statistics object
         * upon the call to execute(). They are cleared after writing
         * the content.
         */
        unsigned int total_iteration_count;
        unsigned int number_of_solves;
    };
  }
}


#endif
