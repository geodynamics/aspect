/*
  Copyright (C) 2016 by the authors of the ASPECT code.

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

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

#include <aspect/postprocess/particles.h>
#include <aspect/particle/world.h>
#include <aspect/simulator.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    ParticleCountStatistics<dim>::execute (TableHandler &statistics)
    {
      const Postprocess::Particles<dim> *particle_postprocessor = this->template find_postprocessor<Postprocess::Particles<dim> >();

      AssertThrow(particle_postprocessor != 0,
                  ExcMessage("The <particles> postprocessor was not found in the list of "
                             "active postprocessors. You need to select this postprocessor to "
                             "be able to select the <particle count> visualization plugin."));

      const std::multimap<aspect::Particle::types::LevelInd, aspect::Particle::Particle<dim> > &particles =
        particle_postprocessor->get_particle_world().get_particles();

      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();

      unsigned int local_min_particles = std::numeric_limits<unsigned int>::max();
      unsigned int local_max_particles = 0;
      const Particle::types::particle_index global_particles = particle_postprocessor->get_particle_world().n_global_particles();

      // compute local min/max
      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
          {
            const aspect::Particle::types::LevelInd current_cell (cell->level(),cell->index());

            const unsigned int particles_in_cell = particles.count(current_cell);
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
