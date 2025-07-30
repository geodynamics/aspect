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


#include <aspect/postprocess/timing_statistics.h>
#include <aspect/simulator.h>

namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    TimingStatistics<dim>::execute (TableHandler &statistics)
    {
      const auto &timer = this->get_computing_timer();
      const std::map<std::string, double> &timing_map = timer.get_summary_data(TimerOutput::total_wall_time);

      for (const auto &section: timing_map)
        {
          const std::string section_name = "Total wall time: " + section.first;
          statistics.add_value(section_name,
                               section.second);
          statistics.set_scientific(section_name, true);
        }

      return std::make_pair (std::string(),std::string());
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(TimingStatistics,
                                  "timing statistics",
                                  "A postprocessor that outputs timing information in "
                                  "the statistics file, in particular the total wall time spent "
                                  "in the different timing sections until the current timestep. "
                                  "Note that this postprocessor cannot report accurate "
                                  "timings for the postprocessing section as it is "
                                  "executed before postprocessing is complete.")
  }
}
