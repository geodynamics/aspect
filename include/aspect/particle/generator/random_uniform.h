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

#ifndef _aspect_particle_generator_random_uniform_h
#define _aspect_particle_generator_random_uniform_h

#include <aspect/particle/generator/interface.h>

namespace aspect
{
  namespace Particle
  {
    namespace Generator
    {
      /**
       * Generate random uniform distribution of particles over entire simulation domain.
       *
       * @ingroup ParticleGenerators
       */
      template <int dim>
      class RandomUniform : public Interface<dim>
      {
        public:
          /**
           * Generate a set of particles in the given
           * particle handler. The particles are generated
           * with a uniform density in the whole domain.
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
      };

    }
  }
}

#endif
