/*
 Copyright (C) 2015 - 2022 by the authors of the ASPECT code.

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

#ifndef _aspect_particle_generator_random_uniform_h
#define _aspect_particle_generator_random_uniform_h

#include <aspect/particle/generator/probability_density_function.h>

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

          // avoid -Woverloaded-virtual
          // TODO: remove this using directive once the following deprecated
          // function in the interface class has been removed:
          // generate_particles(std::multimap<Particles::internal::LevelInd, Particle<dim>> &particles)
          using Generator::Interface<dim>::generate_particles;

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
           */
          types::particle_index n_particles;

          /**
           * If true, particle numbers per cell are calculated randomly
           * according to their respective probability density. If false,
           * first determine how many particles each cell should have based
           * on the integral of the density over each of the cells, and then
           * once we know how many particles we want on each cell, choose their
           * locations randomly within each cell.
           */
          bool random_cell_selection;

          /**
           * The seed for the random number generator that controls the
           * particle generation.
           */
          unsigned int random_number_seed;
      };

    }
  }
}

#endif
