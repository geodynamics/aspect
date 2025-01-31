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
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include <aspect/postprocess/particle_count_statistics.h>
#include <aspect/particle/manager.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    ParticleCountStatistics<dim>::execute (TableHandler &statistics)
    {
      unsigned int local_min_particles = std::numeric_limits<unsigned int>::max();
      unsigned int local_max_particles = 0;

      types::particle_index global_particles = 0;
      for (unsigned int particle_manager_index = 0; particle_manager_index < this->n_particle_managers(); ++particle_manager_index)
        global_particles += this->get_particle_manager(particle_manager_index).n_global_particles();

      // compute local min/max
      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          {
            unsigned int particles_in_cell = 0;
            for (unsigned int particle_manager_index = 0; particle_manager_index < this->n_particle_managers(); ++particle_manager_index)
              particles_in_cell += this->get_particle_manager(particle_manager_index).get_particle_handler().n_particles_in_cell(cell);

            local_min_particles = std::min(local_min_particles,particles_in_cell);
            local_max_particles = std::max(local_max_particles,particles_in_cell);
          }

      // now do the reductions over all processors
      const unsigned int global_min_particles = Utilities::MPI::min (local_min_particles,
                                                                     this->get_mpi_communicator());
      const unsigned int global_max_particles = Utilities::MPI::max (local_max_particles,
                                                                     this->get_mpi_communicator());

      // finally produce something for the statistics file
      statistics.add_value ("Minimal particles per cell: ", global_min_particles);
      statistics.add_value ("Average particles per cell: ", global_particles / this->get_triangulation().n_global_active_cells());
      statistics.add_value ("Maximal particles per cell: ", global_max_particles);

      std::ostringstream output;
      output << global_min_particles << ", "
             << global_particles / this->get_triangulation().n_global_active_cells() << ", "
             << global_max_particles;

      return std::pair<std::string, std::string> ("Particle count per cell min/avg/max:",
                                                  output.str());
    }



    template <int dim>
    std::list<std::string>
    ParticleCountStatistics<dim>::required_other_postprocessors() const
    {
      return {"particles"};
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(ParticleCountStatistics,
                                  "particle count statistics",
                                  "A postprocessor that computes some statistics about "
                                  "the particle distribution, if present in this simulation. "
                                  "In particular, it computes minimal, average and maximal "
                                  "values of particles per cell in the global domain.")
  }
}
