/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

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

#include <aspect/adjoint/manager.h>
#include <aspect/adjoint/kernel_calculator.h>
#include <aspect/postprocess/adjoint_kernels.h>

#include <fstream>
#include <iomanip>
#include <map>
#include <string>

namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    AdjointKernels<dim>::execute(TableHandler &statistics)
    {
      const Adjoint::Manager<dim> &adjoint_manager = this->get_adjoint_manager();
      const Adjoint::KernelRepository<dim> &kernels = adjoint_manager.get_kernels();
      const std::map<std::string, double> objective_values =
        adjoint_manager.get_objective_values();

      statistics.add_value("Number of adjoint kernel contributions",
                           kernels.n_contributions());
      for (const auto &objective_value : objective_values)
        statistics.add_value("Adjoint objective " + objective_value.first,
                             objective_value.second);

      if (kernels.empty())
        return std::make_pair(std::string("Writing adjoint kernels"),
                              std::string("not available"));

      const unsigned int mpi_rank = Utilities::MPI::this_mpi_process(this->get_mpi_communicator());
      const std::string filename = this->get_output_directory()
                                   + "adjoint_kernels_rank_"
                                   + Utilities::int_to_string(mpi_rank, 5)
                                   + ".txt";

      std::ofstream out(filename);
      out << std::setprecision(16);
      out << "# time " << this->get_time() << "\n";
      out << "# timestep " << this->get_timestep_number() << "\n";
      for (const auto &objective_value : objective_values)
        out << "# objective_value\t" << objective_value.first << "\t" << objective_value.second << "\n";
      out << "# columns: objective term property active_cell_index cell_volume value\n";

      const Vector<double> &cell_volumes = kernels.cell_volumes();

      for (const auto &entry : kernels.contributions())
        {
          const Adjoint::KernelContributionKey &key = entry.first;
          const Vector<double> &values = entry.second;
          double local_integral = 0.0;

          for (const auto &cell : this->get_triangulation().active_cell_iterators())
            if (cell->is_locally_owned())
              {
                const unsigned int cell_index = cell->active_cell_index();
                local_integral += values(cell_index) * cell_volumes(cell_index);
              }

          const double integral = Utilities::MPI::sum(local_integral, this->get_mpi_communicator());
          out << "# contribution_integral\t"
              << key.objective_name << "\t"
              << key.physics_term_name << "\t"
              << Adjoint::property_name(key.property) << "\t"
              << integral << "\n";

          for (const auto &cell : this->get_triangulation().active_cell_iterators())
            if (cell->is_locally_owned())
              {
                const unsigned int cell_index = cell->active_cell_index();
                out << key.objective_name << "\t"
                    << key.physics_term_name << "\t"
                    << Adjoint::property_name(key.property) << "\t"
                    << cell_index << "\t"
                    << cell_volumes(cell_index) << "\t"
                    << values(cell_index) << "\n";
              }

          out << "\n";
        }

      return std::make_pair(std::string("Writing adjoint kernels"),
                            filename);
    }
  }
}

namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(AdjointKernels,
                                  "adjoint kernels",
                                  "A postprocessor that writes adjoint kernel contributions assembled by the instantaneous Stokes adjoint workflow.")
  }
}
