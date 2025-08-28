/*
  Copyright (C) 2025 - by the authors of the ASPECT code.

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
