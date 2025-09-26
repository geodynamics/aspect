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


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    ParticleDistributionScore<dim>::execute (TableHandler &statistics)
    {
      double local_min_score = std::numeric_limits<double>::max();
      double local_max_score = 0;
      unsigned int number_of_cells = 0;
      std::vector<double> cell_scores;

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        {
          if (cell->is_locally_owned())
            {
              /*
                These are the steps taken to compute the density score:
                Create a table with perfect distribution containing an even number of particles in each bucket.
                Create a table with the "worst case" distribution in which all particles are in one bucket.
                Take the distance squared between the "perfect" table and the worst case table, treating both tables as vectors with granularity^dim elements.
                Create a table with the actual distribution.
                Take distance between "perfect" table and the actual table.
                Return the ratio between the observed distance and the worst case distance.
                A value of 1 represents the worst possible distribution, a value of 0 represents a completely even distribution.
              */
              unsigned int particles_in_cell = 0;
              for (unsigned int particle_manager_index = 0; particle_manager_index < this->n_particle_managers(); ++particle_manager_index)
                {
                  particles_in_cell += this->get_particle_manager(particle_manager_index).get_particle_handler().n_particles_in_cell(cell);
                }
              if (particles_in_cell > 0)
                {
                  ++number_of_cells;
                  const double particles_in_cell_double = static_cast<double>(particles_in_cell);
                  const double granularity_double = static_cast<double>(granularity);
                  const double ideal_n_particles_per_bucket = particles_in_cell_double/(Utilities::fixed_power<dim>(granularity_double));

                  /*
                  The "buckets_ideal" table contains doubles because the ideal
                  number of particles per bucket will not always be an integer.
                  */
                  Table<dim,double> buckets_ideal;
                  TableIndices<dim> bucket_sizes;
                  for (unsigned int i=0; i<dim; ++i)
                    bucket_sizes[i] = granularity;

                  buckets_ideal.reinit(bucket_sizes);
                  const double bucket_width = 1.0/granularity_double;
                  buckets_ideal.fill(ideal_n_particles_per_bucket);

                  Table<dim,unsigned int> buckets_actual;
                  buckets_actual.reinit(bucket_sizes);
                  sort_particles_into_buckets(cell,bucket_width,buckets_actual);

                  /*
                  In the worst case, all particles are in one bucket.
                  (granularity^dim)-1 is equal to the number of buckets
                  in the table minus 1 (the bucket all the particles are in).
                  */
                  const double worst_case_empty_buckets = static_cast<double>(Utilities::fixed_power<dim>(granularity_double))-1.0;
                  const double worst_case_error_squared =
                    ((particles_in_cell_double-ideal_n_particles_per_bucket)*(particles_in_cell_double-ideal_n_particles_per_bucket))+
                    ((ideal_n_particles_per_bucket*ideal_n_particles_per_bucket)*(worst_case_empty_buckets));

                  double actual_error_squared = 0;
                  for (unsigned int x=0; x<granularity; ++x)
                    {
                      for (unsigned int y=0; y<granularity; ++y)
                        {
                          TableIndices<dim> entry_index;
                          entry_index[0] = x;
                          entry_index[1] = y;
                          // do another loop if in 3d
                          if (dim == 3)
                            {
                              for (unsigned int z=0; z<granularity; ++z)
                                {
                                  entry_index[2] = z;
                                  const double value_ideal = buckets_ideal(entry_index);
                                  const double value_actual = static_cast<double>(buckets_actual(entry_index));
                                  actual_error_squared += (value_ideal - value_actual)*(value_ideal - value_actual);
                                }
                            }
                          else
                            {
                              const double value_ideal = buckets_ideal(entry_index);
                              const double value_actual = static_cast<double>(buckets_actual(entry_index));
                              actual_error_squared += (value_ideal - value_actual)*(value_ideal - value_actual);
                            }
                        }
                    }

                  // Take the ratio between the actual error and the worst case
                  // error, resulting in a score from 0 to 1 for the cell.
                  const double distribution_score_current_cell = actual_error_squared/worst_case_error_squared;

                  cell_scores.push_back(distribution_score_current_cell);

                  local_max_score = std::max(local_max_score, distribution_score_current_cell);
                  local_min_score = std::min(local_min_score, distribution_score_current_cell);
                }
              else if (particles_in_cell == 0)
                {
                  ++number_of_cells;
                  // The score should be bad if there are no particles in a cell
                  const double distribution_score_current_cell = 1.0;
                  cell_scores.push_back(distribution_score_current_cell);
                  local_max_score = std::max(local_max_score, distribution_score_current_cell);
                  local_min_score = std::min(local_min_score, distribution_score_current_cell);
                }
            }
        }

      // Calculate the mean and standard deviation of cell scores across all processors
      const std::pair<double, double> mean_and_standard_deviation =
        Utilities::MPI::mean_and_standard_deviation(cell_scores.begin(),
                                                    cell_scores.end(),
                                                    this->get_mpi_communicator());

      // get final values for min and max score from all processors
      const double global_max_score = Utilities::MPI::max (local_max_score, this->get_mpi_communicator());
      const double global_min_score = Utilities::MPI::min (local_min_score, this->get_mpi_communicator());

      // write to statistics file
      statistics.add_value ("Minimal particle distribution score: ", global_min_score);
      statistics.add_value ("Average particle distribution score: ", mean_and_standard_deviation.first);
      statistics.add_value ("Maximal particle distribution score: ", global_max_score);
      statistics.add_value ("Cell Score Standard Deviation: ", mean_and_standard_deviation.second);

      std::ostringstream output;
      output << global_min_score << '/' << mean_and_standard_deviation.first << '/' << global_max_score << '/' << mean_and_standard_deviation.second;

      return std::pair<std::string, std::string> ("Particle distribution score min/avg/max/stdev:",
                                                  output.str());
    }



    template <int dim>
    std::list<std::string>
    ParticleDistributionScore<dim>::required_other_postprocessors() const
    {
      return {"particles"};
    }



    template <int dim>
    void ParticleDistributionScore<dim>::sort_particles_into_buckets(
      const typename Triangulation<dim>::active_cell_iterator &cell,
      const double bucket_width,
      Table<dim,unsigned int> &buckets) const
    {
      for (unsigned int particle_manager_index = 0; particle_manager_index < this->n_particle_managers(); ++particle_manager_index)
        {
          // sort the particles within the current cell
          const typename Particle::ParticleHandler<dim>::particle_iterator_range particle_range =
            this->get_particle_manager(particle_manager_index).get_particle_handler().particles_in_cell(cell);

          for (const auto &particle: particle_range)
            {
              const double particle_x = particle.get_reference_location()[0];
              const double particle_y = particle.get_reference_location()[1];

              const double x_ratio = (particle_x) / (bucket_width);
              const double y_ratio = (particle_y) / (bucket_width);

              unsigned int x_index = static_cast<unsigned int>(std::floor(x_ratio));
              unsigned int y_index = static_cast<unsigned int>(std::floor(y_ratio));

              /*
              If a particle is exactly on the boundary of two cells its
              reference location will equal 1, and if this is the case,
              the "x/y/z_index" will be outside of the range of the table without
              these checks. The table has a number of entries equal to "granularity" in each dimension,
              and the table is indexed at 0, so if the "x/y/z_indez" equals "granularity" it
              will be out of range.
              */
              if (x_index == granularity)
                x_index = granularity-1;
              if (y_index == granularity)
                y_index = granularity-1;

              TableIndices<dim> entry_index;
              entry_index[0] = x_index;
              entry_index[1] = y_index;
              if (dim == 3)
                {
                  const double particle_z = particle.get_reference_location()[2];
                  const double z_ratio = (particle_z) / (bucket_width);
                  unsigned int z_index = static_cast<unsigned int>(std::floor(z_ratio));
                  if (z_index == granularity)
                    z_index = granularity-1;
                  entry_index[2] = z_index;
                }

              ++buckets(entry_index);
            }
        }
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
