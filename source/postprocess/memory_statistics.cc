/*
 Copyright (C) 2016 - 2019 by the authors of the ASPECT code.

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


#include <aspect/postprocess/memory_statistics.h>

#include <aspect/simulator.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    MemoryStatistics<dim>::execute (TableHandler &statistics)
    {
      // memory consumption:
      const double mb = 1024*1024; // convert from bytes into mb
      statistics.add_value ("System matrix memory consumption (MB) ", this->get_system_matrix().memory_consumption()/mb);
      statistics.add_value ("Triangulation memory consumption (MB) ", this->get_triangulation().memory_consumption()/mb);
      statistics.add_value ("p4est memory consumption (MB) ", this->get_triangulation().memory_consumption_p4est()/mb);
      statistics.add_value ("DoFHandler memory consumption (MB) ", this->get_dof_handler().memory_consumption()/mb);
      statistics.add_value ("current_constraints memory consumption (MB) ", this->get_current_constraints().memory_consumption()/mb);
      statistics.add_value ("Solution vector memory consumption (MB) ", this->get_solution().memory_consumption()/mb);

      if (output_vmpeak)
        {
          // allow disabling of the output because this is not stable in automated tests:
          dealii::Utilities::System::MemoryStats stats;
          dealii::Utilities::System::get_memory_stats(stats);
          const double max_vmpeak = dealii::Utilities::MPI::max(stats.VmPeak/1024.0, this->get_mpi_communicator());
          statistics.add_value ("Peak virtual memory usage (VmPeak) (MB) ", max_vmpeak);
        }

      std::ostringstream output;
      output << std::fixed << std::setprecision(2) << this->get_system_matrix().memory_consumption()/mb << " MB";

      return std::pair<std::string, std::string> ("System matrix memory consumption: ",
                                                  output.str());

    }




    template <int dim>
    void
    MemoryStatistics<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Memory statistics");
        {
          prm.declare_entry ("Output peak virtual memory (VmPeak)", "true",
                             Patterns::Bool(),
                             "If set to 'true', also output the peak virtual memory "
                             "usage (computed as the maximum over all processors).");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    MemoryStatistics<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Memory statistics");
        {
          output_vmpeak = prm.get_bool ("Output peak virtual memory (VmPeak)");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(MemoryStatistics,
                                  "memory statistics",
                                  "A postprocessor that computes some statistics about "
                                  "the memory consumption. "
                                  "In particular, it computes the memory usage of the "
                                  "system matrix, triangulation, p4est, "
                                  "DoFHandler, current constraints, solution vector, "
                                  "and peak virtual memory usage, all in MB. "
                                  "It also outputs the memory usage of the system "
                                  "matrix to the screen.")
  }
}
