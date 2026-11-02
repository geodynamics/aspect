// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2015-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

#ifndef _aspect_particle_generator_probability_density_function_h
#define _aspect_particle_generator_probability_density_function_h

#include <aspect/particle/generator/interface.h>

#include <deal.II/base/parsed_function.h>

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
           * Generate a set of particles in the given
           * particle handler. The particle density is set by an analytically
           * prescribed density function that is set as an input parameter.
           *
           * @param [in,out] particle_handler The particle handler into which
           * the generated particles should be inserted.
           */
          void
          generate_particles(Particles::ParticleHandler<dim> &particle_handler) override;

          /**
           * Declare the parameters this class takes through input files.
           */
          static
          void
          declare_parameters (ParameterHandler &prm);

          /**
           * Read the parameters this class declares from the parameter file.
           */
          void
          parse_parameters (ParameterHandler &prm) override;

        private:
          /**
           * Number of particles to create
           *
           * This variable is read from the parameter file through a parameter
           * called 'Number of particles'.
           */
          types::particle_index n_particles;

          /**
           * If true, particle numbers per cell are calculated randomly
           * according to their respective probability density. If false,
           * first determine how many particles each cell should have based
           * on the integral of the density over each of the cells, and then
           * once we know how many particles we want on each cell, choose their
           * locations randomly within each cell.
           * This variable is read from the parameter file through a parameter
           * called 'Random cell selection'.
           */
          bool random_cell_selection;

          /**
           * The seed for the random number generator that controls the
           * particle generation.
           *
           * This variable is read from the parameter file through a parameter
           * called 'Random number seed'.
           */
          unsigned int random_number_seed;

          /**
           * A function object representing the particle location probability
           * density.
           */
          Functions::ParsedFunction<dim> function;
      };

    }
  }
}

#endif
