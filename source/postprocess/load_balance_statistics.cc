/*
  Copyright (C) 2016 - 2024 by the authors of the ASPECT code.

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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include <aspect/postprocess/load_balance_statistics.h>

#include <aspect/particle/manager.h>

namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    LoadBalanceStatistics<dim>::execute (TableHandler &statistics)
    {
      const unsigned int locally_owned_active_cells = this->get_triangulation().n_locally_owned_active_cells();
      const dealii::Utilities::MPI::MinMaxAvg cell_distribution
        = dealii::Utilities::MPI::min_max_avg(locally_owned_active_cells,
                                              this->get_mpi_communicator());

      statistics.add_value ("Minimal cells per process",
                            static_cast<unsigned int>(cell_distribution.min));
      statistics.add_value ("Maximal cells per process",
                            static_cast<unsigned int>(cell_distribution.max));
      statistics.add_value ("Average cells per process",
                            cell_distribution.avg);

      if (this->n_particle_managers() > 0)
        {
          types::particle_index locally_owned_particles = 0;
          for (unsigned int particle_manager_index = 0; particle_manager_index < this->n_particle_managers(); ++particle_manager_index)
            locally_owned_particles += this->get_particle_manager(particle_manager_index).
                                       get_particle_handler().n_locally_owned_particles();

          const dealii::Utilities::MPI::MinMaxAvg particles_per_process =
            dealii::Utilities::MPI::min_max_avg(locally_owned_particles, this->get_mpi_communicator());

          statistics.add_value ("Minimal particles per process",
                                static_cast<types::particle_index>(particles_per_process.min));
          statistics.add_value ("Maximal particles per process",
                                static_cast<types::particle_index>(particles_per_process.max));
          statistics.add_value ("Average particles per process",
                                particles_per_process.avg);

          const double local_ratio = (locally_owned_active_cells != 0)
                                     ?
                                     static_cast<double>(locally_owned_particles)
                                     / static_cast<double>(locally_owned_active_cells)
                                     :
                                     0.0;

          const dealii::Utilities::MPI::MinMaxAvg particle_to_cell_ratio
            = dealii::Utilities::MPI::min_max_avg(local_ratio,
                                                  this->get_mpi_communicator());

          statistics.add_value ("Minimal local particle to cell ratio", particle_to_cell_ratio.min);
          statistics.add_value ("Maximal local particle to cell ratio", particle_to_cell_ratio.max);
          statistics.add_value ("Average local particle to cell ratio", particle_to_cell_ratio.avg);
        }

      std::ostringstream output;
      output.precision(4);
      output << cell_distribution.min << '/'
             << cell_distribution.max << '/'
             << cell_distribution.avg;

      return std::pair<std::string, std::string> ("Cells per process min/max/avg:",
                                                  output.str());
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(LoadBalanceStatistics,
                                  "load balance statistics",
                                  "A postprocessor that computes statistics about "
                                  "the distribution of cells, and if present particles "
                                  "across subdomains. "
                                  "In particular, it computes maximal, average and "
                                  "minimal number of cells across all ranks. "
                                  "If there are particles it also computes the "
                                  "maximal, average, and minimum number of particles across "
                                  "all ranks, and maximal, average, and minimal ratio "
                                  "between local number of particles and local number "
                                  "of cells across all processes. All of these numbers "
                                  "can be useful to assess the load balance between "
                                  "different MPI ranks, as the difference between the "
                                  "minimal and maximal load should be as small as "
                                  "possible.")
  }
}
