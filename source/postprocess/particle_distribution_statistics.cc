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

#include <aspect/postprocess/particle_distribution_statistics.h>

namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    ParticleDistributionStatistics<dim>::execute (TableHandler &statistics)
    {
      unsigned int cells_with_particles = 0;
      double standard_deviation_sum = 0;
      double standard_deviation_min = std::numeric_limits<double>::max();
      double standard_deviation_max = std::numeric_limits<double>::min();
      double function_min_sum = 0;
      double function_max_sum = 0;
      double function_min_min = std::numeric_limits<double>::max();
      double function_max_max = std::numeric_limits<double>::min();


      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        {
          if (cell->is_locally_owned())
            {
              unsigned int particles_in_cell = 0;
              for (unsigned int particle_manager_index = 0; particle_manager_index < this->n_particle_managers(); ++particle_manager_index)
                particles_in_cell += this->get_particle_manager(particle_manager_index).get_particle_handler().n_particles_in_cell(cell);

              if (particles_in_cell > 0)
                {
                  ++cells_with_particles;

                  if (KDE_per_particle == false)
                    {
                      /*
                        Loop through every particle manager here, since the ParticlePDF class operates
                        on a single particle_iterator_range. The ParticlePDF class is written to operate
                        on a single particle_iterator_range at a time so that it can be called directly
                        by the particle_manager class to assist in particle load balancing.
                      */
                      for (unsigned int particle_manager_index = 0; particle_manager_index < this->n_particle_managers(); ++particle_manager_index)
                        {
                          ParticlePDF<dim> pdf(granularity,bandwidth,kernel_function);
                          const typename Particle::ParticleHandler<dim>::particle_iterator_range particle_range = this->get_particle_manager(particle_manager_index).get_particle_handler().particles_in_cell(cell);
                          const unsigned int this_manager_particles_in_cell = this->get_particle_manager(particle_manager_index).get_particle_handler().n_particles_in_cell(cell);

                          pdf.fill_from_particle_range(particle_range,this_manager_particles_in_cell);
                          pdf.compute_statistical_values();

                          standard_deviation_min = std::min(standard_deviation_min, pdf.get_standard_deviation());
                          standard_deviation_max = std::max(standard_deviation_max, pdf.get_standard_deviation());
                          standard_deviation_sum += pdf.get_standard_deviation();

                          function_min_sum += pdf.get_min();
                          function_max_sum += pdf.get_max();
                          function_min_min = std::min(function_min_min, pdf.get_min());
                          function_max_max = std::max(function_max_max, pdf.get_max());
                      }
                    }
                  else
                    {
                        for (unsigned int particle_manager_index = 0; particle_manager_index < this->n_particle_managers(); ++particle_manager_index)
                        {
                          ParticlePDF<dim> pdf(bandwidth,kernel_function);
                          const typename Particle::ParticleHandler<dim>::particle_iterator_range particle_range = this->get_particle_manager(particle_manager_index).get_particle_handler().particles_in_cell(cell);
                          const unsigned int this_manager_particles_in_cell = this->get_particle_manager(particle_manager_index).get_particle_handler().n_particles_in_cell(cell);

                          pdf.fill_from_particle_range(particle_range,this_manager_particles_in_cell);
                          pdf.compute_statistical_values();

                          standard_deviation_min = std::min(standard_deviation_min, pdf.get_standard_deviation());
                          standard_deviation_max = std::max(standard_deviation_max, pdf.get_standard_deviation());
                          standard_deviation_sum += pdf.get_standard_deviation();

                          function_min_sum += pdf.get_min();
                          function_max_sum += pdf.get_max();
                          function_min_min = std::min(function_min_min, pdf.get_min());
                          function_max_max = std::max(function_max_max, pdf.get_max());
                        }
                    }
                }
            }
        }

      // Get final values from all processors
      const double global_standard_deviation_max = Utilities::MPI::max (standard_deviation_max, this->get_mpi_communicator());
      const double global_standard_deviation_min = Utilities::MPI::min (standard_deviation_min, this->get_mpi_communicator());
      const double global_cells_with_particles = Utilities::MPI::sum (cells_with_particles, this->get_mpi_communicator());
      const double global_standard_deviation_sum = Utilities::MPI::sum (standard_deviation_sum, this->get_mpi_communicator());
      const double global_standard_deviation_mean = global_standard_deviation_sum/global_cells_with_particles;
      const double global_function_min_min = Utilities::MPI::min(function_min_min,this->get_mpi_communicator());
      const double global_function_max_max = Utilities::MPI::min(function_max_max,this->get_mpi_communicator());

      // Get the average of the functions max and min values
      const double global_function_min_sum = Utilities::MPI::sum (function_min_sum, this->get_mpi_communicator());
      const double global_function_max_sum = Utilities::MPI::sum (function_max_sum, this->get_mpi_communicator());
      const double global_function_min_mean = global_function_min_sum/global_cells_with_particles;
      const double global_function_max_mean = global_function_max_sum/global_cells_with_particles;

      // Write to statistics file
      statistics.add_value ("Minimum PDF standard deviation ", global_standard_deviation_min);
      statistics.add_value ("Mean of PDF standard deviation: ", global_standard_deviation_mean);
      statistics.add_value ("Maximum PDF standard deviation: ", global_standard_deviation_max);
      statistics.add_value ("Mean of PDF minimum values: ", global_function_min_mean);
      statistics.add_value ("Mean PDF maximum values: ", global_function_max_mean);
      statistics.add_value ("Minimum of PDF minimum values: ", global_function_min_min);
      statistics.add_value ("Maximum of PDF maximum values: ", global_function_max_max);


      std::ostringstream output;
      output << global_standard_deviation_min << "/" << global_standard_deviation_mean << "/" << global_standard_deviation_max << ", "
             << global_function_min_mean << "/" << global_function_max_mean << ", " << global_function_min_min << "/" << global_function_max_max;


      return std::pair<std::string, std::string> ("Particle Distribution Stats (stddev min/mean/max, mean min/max, absolute min/max):",
                                                  output.str());
    }



    template <int dim>
    std::list<std::string>
    ParticleDistributionStatistics<dim>::required_other_postprocessors() const
    {
      return {"particles"};
    }



    template <int dim>
    void
    ParticleDistributionStatistics<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Particle distribution statistics");
        {
          prm.declare_entry("Kernel function","Uniform",
                            Patterns::Selection("Uniform|Gaussian|Triangular|CutOffFunctionW1"),
                            "The kernel smoothing function to use for kernel density estimation."
                           );
          prm.declare_entry("KDE granularity","2",
                            Patterns::Integer (1),
                            "The granularity parameter determines how many discrete inputs exist for "
                            "the probability density function generated by the kernel density estimator. "
                            "The domain of the function is multidimensional so the granularity value determines "
                            "the range of inputs in each dimension. For example, a granularity value of 2 "
                            "results in a PDF which is defined for the inputs 0-1 in each of its dimensions. ");
          prm.declare_entry("Kernel bandwidth","0.3",
                            Patterns::Double (0.001),
                            "The bandwidth parameter sets the bandwidth used for the kernel function. "
                            "The size of the bandwidth determines the output of the kernel function each, "
                            "time it is called. Larger bandwidth result in a point-density function which "
                            "has more smoothing because there is more overlap between the domains of the "
                            "individual kernel functions. The opposite is true for smaller bandwidth values."
                           );
          prm.declare_entry("Use KDE per particle","false",
                            Patterns::Bool(),
                            "If `KDE_per_particle` is true, the point-density function is defined at the position "
                            "of every particle in the cell. If it is false, the point-density-function is defined "
                            "on a regular grid throughout the cell."
                           );
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }



    template <int dim>
    void
    ParticleDistributionStatistics<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Particle distribution statistics");
        {
          KDE_per_particle = prm.get_bool("Use KDE per particle");
          granularity = prm.get_integer("KDE granularity");
          bandwidth = prm.get_double("Kernel bandwidth");
          std::string kernel_function_string = prm.get("Kernel function");

          if (kernel_function_string =="Triangular")
            {
              kernel_function = ParticlePDF<dim>::KernelFunctions::triangular;
            }
          else if (kernel_function_string =="Gaussian")
            {
              kernel_function = ParticlePDF<dim>::KernelFunctions::gaussian;
            }
          else if (kernel_function_string == "CutoffFunctionW1")
            {
              kernel_function = ParticlePDF<dim>::KernelFunctions::cutoff_function_w1_dealii;
            }
          else
            {
              kernel_function = ParticlePDF<dim>::KernelFunctions::uniform;
            }
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
    ASPECT_REGISTER_POSTPROCESSOR(ParticleDistributionStatistics,
                                  "Particle Distribution Statistics",
                                  "A postprocessor that computes some statistics about "
                                  "the particle distribution within grid cells. In particular "
                                  "it calculates a point-density function for every cell and derives "
                                  "the maximum, minimum, and standard deviation values for every cell. "
                                  "The postprocessor reports the average of these values from every cell. "
                                  "It also reports the absolute maximum and minimum values across all cells. "
                                  "These sorts of statistics are useful to determine whether schemes that "
                                  "move, add, remove, or otherwise change the number of particles associated "
                                  "with each cell result in a roughly uniform distribution of particles in "
                                  "each cell, or whether particles tend to cluster in some parts of cells "
                                  "leaving other parts mostly empty. For example, comparing the maximum "
                                  "standard deviations of different load balancing schemes applied to a given "
                                  "test case illuminates which load balancing method creates more clustered particles. "
                                  "The maximum and minimum values of these point-density functions are also useful "
                                  "in the same way. These statistical values are computed from the point-density function "
                                  "and make up a quantitative description of particle clustering which is intended to "
                                  "supplement qualitative descriptions of particle clustering. The goal behind "
                                  "describing particle clustering numerically is to assist in developing new methods "
                                  "to add and delete particles which maintain roughly uniform particle density within cells.")
  }
}
