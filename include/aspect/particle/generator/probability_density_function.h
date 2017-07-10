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

#ifndef _aspect_particle_generator_probability_density_function_h
#define _aspect_particle_generator_probability_density_function_h

#include <aspect/particle/generator/interface.h>

#include <deal.II/base/parsed_function.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#include <boost/random.hpp>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

namespace aspect
{
  namespace Particle
  {
    namespace Generator
    {
      /**
       * Exception
       */
      template <int dim>
      DeclException1 (ProbabilityFunctionNegative,
                      Point<dim>,
                      << "Your probability density function in the particle generator "
                      "returned a negative probability density for the following position: "
                      << arg1 << ". Please check your function expression.");


      /**
       * Generates a random distribution of particles over the simulation
       * domain. The particle density is determined by a user-defined
       * probability density function in the parameter file. This is done using
       * a "roulette wheel" style selection. Every cell is weighted by the value
       * of the provided function at its center multiplied with the cell
       * volume. Then a map between the accumulated cell weight
       * and the cell index of the current cell is constructed. Consequently,
       * a random number between zero and the global integral of the
       * probability density function uniquely defines one particular cell.
       * Afterwards, every process generates n_global_particles random numbers,
       * but only generates a particle if it is the owner of the active cell
       * that is associated with this random number.
       *
       * @ingroup ParticleGenerators
       */
      template <int dim>
      class ProbabilityDensityFunction : public Interface<dim>
      {
        public:
          /**
           * Generate a set of particles in the current
           * particle world. The particle density is set by an analytically
           * prescribed density function that is set as an input parameter.
           * This function builds a list of probabilities for all local cells
           * and then calls generate_particles_in_subdomain() to generate
           * the local particles.
           *
           * @param [in,out] particles A multimap between cells and their
           * particles. This map will be filled in this function.
           */
          virtual
          void
          generate_particles(std::multimap<types::LevelInd, Particle<dim> > &particles);

          /**
           * Declare the parameters this class takes through input files.
           */
          static
          void
          declare_parameters (ParameterHandler &prm);

          /**
           * Read the parameters this class declares from the parameter file.
           */
          virtual
          void
          parse_parameters (ParameterHandler &prm);

        protected:

          /**
           * Returns the weight of one cell, which is interpreted as the probability
           * to generate particles in this cell.
           */
          virtual
          double
          get_cell_weight (typename DoFHandler<dim>::active_cell_iterator &cell) const;

        private:
          /**
           * Number of particles to create
           */
          types::particle_index n_particles;

          /**
           * A function object representing the particle location probability
           * density.
           */
          Functions::ParsedFunction<dim> function;

          /**
           * Generate a set of particles distributed within the local domain
           * according to the probability density function.
           *
           * @param [in] particles_per_cell A vector with n_locally_owned_cells entries
           * that determines how many particles are generated in each cell.
           * @param [in] local_start_id The starting ID to assign to generated particles of the local process.
           * @param [in] n_local_particles The total number of particles to generate locally.
           * @param [out] particles A map between cells and all generated particles.
           *
           */
          void
          generate_particles_in_subdomain (const std::vector<unsigned int> &particles_per_cell,
                                           const types::particle_index local_start_id,
                                           const types::particle_index n_local_particles,
                                           std::multimap<types::LevelInd, Particle<dim> > &particles);

          /**
           * This function loops over all active cells in the local subdomain
           * and returns a vector of accumulated cell weights with the size
           * n_locally_owned_active_cells(). This vector is calculated
           * by looping over all locally owned cells and accumulating the
           * return value of get_cell_weight(cell).
           */
          std::vector<double>
          compute_local_accumulated_cell_weights () const;
      };

    }
  }
}

#endif
