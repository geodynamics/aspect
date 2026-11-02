// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2025 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#ifndef _aspect_postprocess_entropy_statistics_h
#define _aspect_postprocess_entropy_statistics_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace Postprocess
  {

    /**
     * A postprocessor that computes statistics about the temperature equilibration
     * when using the entropy model with multiple components.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class EntropyStatistics : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Connect the callback functions to the respective signals.
         */
        void initialize() override;

        /**
         * Write the entropy statistics into the statistics file.
         */
        std::pair<std::string,std::string>
        execute (TableHandler &statistics) override;

      private:
        /**
         * Callback function that is connected to the post_multicomponent_equilibrium
         * signal to store the solver history.
         */
        void
        store_entropy_solver_history(unsigned int iteration_count);

        /**
         * Variables that store the multicomponent equilibration solver history of the current
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
