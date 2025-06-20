/*
  Copyright (C) 2025 by the authors of the ASPECT code.

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

#include <aspect/postprocess/entropy_statistics.h>
#include <aspect/material_model/entropy_model.h>

namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    void
    EntropyStatistics<dim>::initialize()
    {
      const MaterialModel::EntropyModel<dim> &material_model = Plugins::get_plugin_as_type<const MaterialModel::EntropyModel<dim>>(this->get_material_model());

      material_model.post_multicomponent_equilibrium.connect(
        [&](const SimulatorAccess<dim> &/*simulator_access*/,
            const unsigned int iteration_count)
      {
        this->store_entropy_solver_history(iteration_count);
      });

      total_iteration_count = 0;
      number_of_solves = 0;
    }



    template <int dim>
    void
    EntropyStatistics<dim>::store_entropy_solver_history(unsigned int iteration_count)
    {
      total_iteration_count += iteration_count;
      number_of_solves += 1;
    }



    template <int dim>
    std::pair<std::string,std::string>
    EntropyStatistics<dim>::execute (TableHandler &statistics)
    {
      // Collect the data from all MPI ranks for calculating the average iteration count.
      std::vector<unsigned int> local_iteration_count = {number_of_solves, total_iteration_count};
      std::vector<unsigned int> iteration_count (2);

      Utilities::MPI::sum (local_iteration_count, this->get_mpi_communicator(), iteration_count);

      // average iteration per solve = total_iteration_count / number_of_solves
      const int average_iteration_count = iteration_count[0] > 0
                                          ?
                                          iteration_count[1] / iteration_count[0]
                                          :
                                          iteration_count[1];

      statistics.add_value("Average iterations for multicomponent entropy averaging",
                           average_iteration_count);

      total_iteration_count = 0;
      number_of_solves = 0;

      return std::make_pair (std::string(),std::string());
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(EntropyStatistics,
                                  "entropy statistics",
                                  "A postprocessor that computes the number of iterations used "
                                  "for equilibrating temperature when using the entropy material model "
                                  "with multiple components. It returns an average of iterations performed "
                                  "for every time step during the model evolution.")
  }
}
