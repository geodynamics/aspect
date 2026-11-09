// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2016-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#include <aspect/postprocess/memory_statistics.h>
#include <aspect/simulator/solver/stokes_matrix_free.h>


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

      statistics.add_value ("System matrix memory consumption (MB)", this->get_system_matrix().memory_consumption()/mb);
      statistics.add_value ("Triangulation memory consumption (MB)", this->get_triangulation().memory_consumption()/mb);
      statistics.add_value ("p4est memory consumption (MB)", this->get_triangulation().memory_consumption_p4est()/mb);

      double dof_handler_mem = this->get_dof_handler().memory_consumption();
      double constraints_mem = this->get_current_constraints().memory_consumption();
      if (this->is_stokes_matrix_free())
        {
          dof_handler_mem += this->get_stokes_matrix_free().get_dof_handler_memory_consumption();
          constraints_mem += this->get_stokes_matrix_free().get_constraint_memory_consumption();
        }
      statistics.add_value ("DoFHandler memory consumption (MB)", dof_handler_mem/mb);
      statistics.add_value ("AffineConstraints<double> memory consumption (MB)", constraints_mem/mb);

      statistics.add_value ("Solution vector memory consumption (MB)", this->get_solution().memory_consumption()/mb);

      if (this->is_stokes_matrix_free())
        {
          const double mg_transfer_mem = this->get_stokes_matrix_free().get_mg_transfer_memory_consumption();
          statistics.add_value ("MGTransfer memory consumption (MB)", mg_transfer_mem/mb);

          const double cell_data_mem = this->get_stokes_matrix_free().get_cell_data_memory_consumption();
          statistics.add_value ("Matrix-free cell data memory consumption (MB)", cell_data_mem/mb);
        }

      if (output_vmpeak)
        {
          // allow disabling of the output because this is not stable in automated tests:
          dealii::Utilities::System::MemoryStats stats;
          dealii::Utilities::System::get_memory_stats(stats);
          const double max_vmpeak = dealii::Utilities::MPI::max(stats.VmPeak/1024.0, this->get_mpi_communicator());
          statistics.add_value ("Peak virtual memory usage (VmPeak) (MB)", max_vmpeak);
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
