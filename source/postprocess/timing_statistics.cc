// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2025-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#include <aspect/postprocess/timing_statistics.h>

namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    TimingStatistics<dim>::execute (TableHandler &statistics)
    {
      // The timer only allows access to its data outside of
      // subsections, so we can't just ask for its summary data and
      // output that -- we need to leave the current section, get
      // access to the timing data, then re-enter the section. This is
      // awkward because it means that it will look like we entered
      // the section twice. Rather, just copy the whole timer object,
      // and leave the section in the copy, leaving the original one
      // untouched, and then output the information from the copied
      // object:
      TimerOutput copy_of_timer = this->get_computing_timer();
      copy_of_timer.leave_subsection("Postprocessing");
      const std::map<std::string, double> &timing_map
        = copy_of_timer.get_summary_data(TimerOutput::total_wall_time);

      for (const auto &section : timing_map)
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
