/*
  Copyright (C) 2015 - 2016 by the authors of the ASPECT code.

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

#include <aspect/particle/generator/probability_density_function.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/geometry_info.h>

#include <boost/lexical_cast.hpp>


namespace aspect
{
  namespace Particle
  {
    namespace Generator
    {
      template <int dim>
      void
      ProbabilityDensityFunction<dim>::generate_particles(std::multimap<types::LevelInd, Particle<dim> > &particles)
      {
        // Get the local accumulated probabilities for every cell
        const std::vector<double> accumulated_cell_weights = compute_local_accumulated_cell_weights();

        // Sum the local integrals over all nodes
        double local_weight_integral = accumulated_cell_weights.back();
        const double global_weight_integral = Utilities::MPI::sum (local_weight_integral,
                                                                   this->get_mpi_communicator());

        AssertThrow(global_weight_integral > std::numeric_limits<double>::min(),
                    ExcMessage("The integral of the user prescribed probability "
                               "density function over the domain equals zero, "
                               "ASPECT has no way to determine the cell of "
                               "generated particles. Please ensure that the "
                               "provided function is positive in at least a "
                               "part of the domain, also check the syntax of "
                               "the function."));

        // Determine the starting weight of this process, which is the sum of
        // the weights of all processes with a lower rank
        double local_start_weight = 0.0;
        MPI_Scan(&local_weight_integral, &local_start_weight, 1, MPI_DOUBLE, MPI_SUM, this->get_mpi_communicator());
        local_start_weight -= local_weight_integral;

        // Calculate start id and number of local particles
        const types::particle_index start_id = llround(static_cast<double> (n_particles)  * local_start_weight / global_weight_integral);
        const types::particle_index n_local_particles = llround(static_cast<double> (n_particles) * local_weight_integral / global_weight_integral);

        // Uniform distribution on the interval [0,local_weight_integral). This
        // will be used to select cells for all particles.
        boost::random::uniform_real_distribution<double> uniform_distribution(0.0, local_weight_integral);

        // Loop over all particles to create locally and pick their cells
        std::vector<unsigned int> particles_per_cell(this->get_triangulation().n_locally_owned_active_cells(),0);
        for (types::particle_index current_particle_index = 0; current_particle_index < n_local_particles; ++current_particle_index)
          {
            // Draw the random number that determines the cell of the particle
            const double random_weight = uniform_distribution(this->random_number_generator);

            const std::vector<double>::const_iterator selected_cell = std::lower_bound(accumulated_cell_weights.begin(),
                                                                                       accumulated_cell_weights.end(),
                                                                                       random_weight);
            const unsigned int cell_index = std::distance(accumulated_cell_weights.begin(),selected_cell);

            ++particles_per_cell[cell_index];
          }

        generate_particles_in_subdomain(particles_per_cell,start_id,n_local_particles,particles);
      }

      template <int dim>
      std::vector<double>
      ProbabilityDensityFunction<dim>::compute_local_accumulated_cell_weights () const
      {
        std::vector<double> accumulated_cell_weights;
        accumulated_cell_weights.reserve(this->get_triangulation().n_locally_owned_active_cells());
        double accumulated_cell_weight = 0.0;

        // compute the integral weight by quadrature
        typename DoFHandler<dim>::active_cell_iterator
        cell = this->get_dof_handler().begin_active(),
        endc = this->get_dof_handler().end();
        for (; cell!=endc; ++cell)
          if (cell->is_locally_owned())
            {
              // get_cell_weight makes sure to return positive values
              accumulated_cell_weight += get_cell_weight(cell);
              accumulated_cell_weights.push_back(accumulated_cell_weight);
            }
        return accumulated_cell_weights;
      }

      template <int dim>
      double
      ProbabilityDensityFunction<dim>::get_cell_weight (typename DoFHandler<dim>::active_cell_iterator &cell) const
      {
        // Evaluate function at all cell midpoints, sort cells according to weight
        const QMidpoint<dim> quadrature_formula;

        // In the simplest case we do not even need a FEValues object, because
        // using cell->center() and cell->measure() would be equivalent. This
        // fails however for higher-order mappings like we use.
        FEValues<dim> fe_values (this->get_mapping(),
                                 this->get_fe(),
                                 quadrature_formula,
                                 update_quadrature_points |
                                 update_JxW_values);

        fe_values.reinit (cell);
        const std::vector<Point<dim> > position = fe_values.get_quadrature_points();
        const double quadrature_point_weight = function.value(position[0]);

        AssertThrow(quadrature_point_weight >= 0.0,
                    ProbabilityFunctionNegative<dim>(position[0]));

        return quadrature_point_weight * fe_values.JxW(0);
      }

      template <int dim>
      void
      ProbabilityDensityFunction<dim>::generate_particles_in_subdomain (const std::vector<unsigned int> &particles_per_cell,
                                                                        const types::particle_index first_particle_index,
                                                                        const types::particle_index n_local_particles,
                                                                        std::multimap<types::LevelInd, Particle<dim> > &particles)
      {
        // Generate particles per cell
        unsigned int cell_index = 0;
        unsigned int current_particle_index = first_particle_index;

        // We first store the generated particles in a vector. Since they are
        // generated cell-by-cell, they will already be sorted in the correct
        // order to be later transferred to the multimap with O(N) complexity.
        // If we would insert them into the multimap one-by-one it would
        // increase the complexity to O(N log(N)).
        std::vector<std::pair<types::LevelInd, Particle<dim> > > local_particles;
        local_particles.reserve(n_local_particles);
        for (typename DoFHandler<dim>::active_cell_iterator cell = this->get_dof_handler().begin_active();
             cell!=this->get_dof_handler().end();
             ++cell)
          if (cell->is_locally_owned())
            {
              for (unsigned int i = 0; i < particles_per_cell[cell_index]; ++i)
                {
                  local_particles.push_back(this->generate_particle(cell,current_particle_index));
                  ++current_particle_index;
                }
              ++cell_index;
            }

        particles.insert(local_particles.begin(),local_particles.end());
      }



      template <int dim>
      void
      ProbabilityDensityFunction<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Particles");
          {
            prm.declare_entry ("Number of particles", "1000",
                               Patterns::Double (0),
                               "Total number of particles to create (not per processor or per element). "
                               "The number is parsed as a floating point number (so that one can "
                               "specify, for example, '1e4' particles) but it is interpreted as "
                               "an integer, of course.");

            prm.enter_subsection("Generator");
            {
              prm.enter_subsection("Probability density function");
              {
                Functions::ParsedFunction<dim>::declare_parameters (prm, 1);
              }
              prm.leave_subsection();
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }


      template <int dim>
      void
      ProbabilityDensityFunction<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Particles");
          {
            n_particles    = static_cast<types::particle_index>(prm.get_double ("Number of particles"));

            prm.enter_subsection("Generator");
            {
              prm.enter_subsection("Probability density function");
              {
                try
                  {
                    function.parse_parameters (prm);
                  }
                catch (...)
                  {
                    std::cerr << "ERROR: FunctionParser failed to parse\n"
                              << "\t'Particle.Generator.Probability density function'\n"
                              << "with expression\n"
                              << "\t'" << prm.get("Function expression") << "'";
                    throw;
                  }
              }
              prm.leave_subsection();
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Particle
  {
    namespace Generator
    {
      ASPECT_REGISTER_PARTICLE_GENERATOR(ProbabilityDensityFunction,
                                         "probability density function",
                                         "Generate a random distribution of "
                                         "particles over the entire simulation domain. "
                                         "The probability density is prescribed in the "
                                         "form of a user-prescribed function. The "
                                         "format of this function follows the syntax "
                                         "understood by the muparser library, see "
                                         "Section~\\ref{sec:muparser-format}. The "
                                         "return value of the function is always "
                                         "checked to be a non-negative probability "
                                         "density but it can be zero in "
                                         "parts of the domain.")
    }
  }
}
