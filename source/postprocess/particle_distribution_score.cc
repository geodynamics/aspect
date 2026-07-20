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


#include <aspect/postprocess/particle_distribution_score.h>
#include <aspect/particle/manager.h>
#include <aspect/particle/histogram.h>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    ParticleDistributionScore<dim>::execute (TableHandler &statistics)
    {
      // These need to be vectors to account for multiple particle managers.
      std::vector<double> local_min_scores(this->n_particle_managers(),std::numeric_limits<double>::max());
      std::vector<double> local_max_scores(this->n_particle_managers(),0);
      std::vector<std::vector<double>> cell_scores(this->n_particle_managers());

      for (auto &scores: cell_scores)
        scores.reserve(this->get_triangulation().n_active_cells());
      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        {
          if (cell->is_locally_owned())
            {
              /*
              The method used to compute the density score is described in the ParticleHistogram class.
              */

              for (unsigned int particle_manager_index = 0; particle_manager_index < this->n_particle_managers(); ++particle_manager_index)
                {
                  const unsigned int particles_in_cell = this->get_particle_manager(particle_manager_index).get_particle_handler().n_particles_in_cell(cell);

                  if (particles_in_cell > 0)
                    {
                      Particle::ParticleHistogram<dim> histogram(granularity);

                      const typename Particle::ParticleHandler<dim>::particle_iterator_range particle_range =
                        this->get_particle_manager(particle_manager_index).get_particle_handler().particles_in_cell(cell);

                      const double distribution_score_current_cell
                        = histogram.evaluate_particle_range(particle_range,particles_in_cell);

                      cell_scores[particle_manager_index].push_back(distribution_score_current_cell);
                      local_max_scores[particle_manager_index] = std::max(local_max_scores[particle_manager_index], distribution_score_current_cell);
                      local_min_scores[particle_manager_index] = std::min(local_min_scores[particle_manager_index], distribution_score_current_cell);
                    }
                  else if (particles_in_cell == 0)
                    {
                      // The score should be bad if there are no particles in a cell
                      const double distribution_score_current_cell = 1.0;
                      cell_scores[particle_manager_index].push_back(distribution_score_current_cell);
                      local_max_scores[particle_manager_index] = std::max(local_max_scores[particle_manager_index], distribution_score_current_cell);
                      local_min_scores[particle_manager_index] = std::min(local_min_scores[particle_manager_index], distribution_score_current_cell);
                    }
                }
            }
        }

      std::ostringstream output;
      for (unsigned int particle_manager_index = 0; particle_manager_index < this->n_particle_managers(); ++particle_manager_index)
        {
          // Calculate the mean and standard deviation of cell scores across all processors
          const std::pair<double, double> mean_and_standard_deviation =
            Utilities::MPI::mean_and_standard_deviation(cell_scores[particle_manager_index].begin(),
                                                        cell_scores[particle_manager_index].end(),
                                                        this->get_mpi_communicator());
          // get final values for min and max score from all processors
          const double global_max_score = Utilities::MPI::max (local_max_scores[particle_manager_index], this->get_mpi_communicator());
          const double global_min_score = Utilities::MPI::min (local_min_scores[particle_manager_index], this->get_mpi_communicator());
          // write to statistics file for all particle managers
          const std::string particle_manager_index_prefix = (particle_manager_index == 0) ? "" : "Particles " + std::to_string(particle_manager_index+1) + ": ";
          statistics.add_value (particle_manager_index_prefix+"Minimal particle distribution score: ", global_min_score);
          statistics.add_value (particle_manager_index_prefix+"Average particle distribution score: ", mean_and_standard_deviation.first);
          statistics.add_value (particle_manager_index_prefix+"Maximal particle distribution score: ", global_max_score);
          statistics.add_value (particle_manager_index_prefix+"Cell Score Standard Deviation: ", mean_and_standard_deviation.second);

          // write screen output for the first particle manager
          if (particle_manager_index == 0)
            output << global_min_score << '/' << mean_and_standard_deviation.first << '/' << global_max_score << '/' << mean_and_standard_deviation.second;
        }
      return std::pair<std::string, std::string> ("Particle distribution score min/avg/max/stdev: ",
                                                  output.str());
    }



    template <int dim>
    std::list<std::string>
    ParticleDistributionScore<dim>::required_other_postprocessors() const
    {
      return {"particles"};
    }



    template <int dim>
    void
    ParticleDistributionScore<dim>::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Particle distribution score");
        {
          prm.declare_entry("Granularity","2",
                            Patterns::Integer (2),
                            "This parameter determines how many bins this postprocessor sorts particles "
                            "into. The ideal value for granularity depends on the maximum number of particles "
                            "in the cell. Generally higher values lead to more accuracy."
                           );
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }



    template <int dim>
    void
    ParticleDistributionScore<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Particle distribution score");
        {
          granularity = prm.get_integer ("Granularity");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }
  }
}



// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(ParticleDistributionScore,
                                  "particle distribution score",
                                  "This postprocessor is intended to help evaluate how much the density of particles "
                                  "varies within cells, with the goal of supporting the development of algorithms "
                                  "to select particles for deletion and addition during load balancing. "
                                  "Because particle deletion can happen in various ways when load balancing, "
                                  "there are times when unintended gradients of particle density can form, especially "
                                  "when particles are moving from more refined cells into coarser cells. "
                                  "The postprocessor calculates a value from 0 to 1 for every "
                                  "cell, with values closer to zero representing cells with a density which is more "
                                  "uniform across the cell's area and values closer to 1 representing a cell in which the particles are "
                                  "highly concentrated in one part of the cell. Essentially, the postprocessor is "
                                  "trying to describe numerically how much particle density varies within cells. "
                                  "It does this by sorting every particle in each cell into a bucket based on "
                                  "the particle's location, and comparing the resulting data structure to ideal "
                                  "and worst case versions.")
  }
}
