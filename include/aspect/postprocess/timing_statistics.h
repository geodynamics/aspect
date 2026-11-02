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


#ifndef _aspect_postprocess_timing_statistics_h
#define _aspect_postprocess_timing_statistics_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{

  namespace Postprocess
  {
    /**
     * A postprocessor that outputs timing information in
     * the statistics file, in particular the total wall time spent
     * in the different timing sections until the current timestep.
     * Note that this postprocessor cannot report accurate
     * timings for the postprocessing section as it is
     * executed before postprocessing is complete.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class TimingStatistics : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Write all timing statistics columns into the statistics object.
         */
        std::pair<std::string,std::string>
        execute (TableHandler &statistics) override;
    };
  }
}


#endif
