/*
 Copyright (C) 2015 by the authors of the ASPECT code.

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
 along with ASPECT; see the file doc/COPYING.  If not see
 <http://www.gnu.org/licenses/>.
 */

#ifndef __aspect__particle_generator_random_uniform_h
#define __aspect__particle_generator_random_uniform_h

#include <aspect/particle/generator/interface.h>

#include <boost/random.hpp>

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
           * Constructor.
           */
          RandomUniform();

          /**
           * Generate a uniformly randomly distributed set of particles in the current triangulation.
           *
           * @return A multimap containing cells and their contained particles.
           */
          virtual
          std::multimap<types::LevelInd, Particle<dim> >
          generate_particles();


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

        private:
          /**
           * Random number generator and an object that describes a
           * uniform distribution on the interval [0,1]. These
           * will be used to generate particle locations at random.
           */
          boost::mt19937            random_number_generator;
          boost::uniform_01<double> uniform_distribution_01;

          /**
           * Number of initial particles to create
           */
          types::particle_index n_tracers;


          /**
           * Generate a set of particles uniformly randomly distributed within the
           * specified triangulation. This is done using "roulette wheel" style
           * selection weighted by cell volume. We do cell-by-cell assignment of
           * particles because the decomposition of the mesh may result in a highly
           * non-rectangular local mesh which makes uniform particle distribution difficult.
           *
           * @param [in] num_particles The number of particles to generate in this subdomain
           * @param [in] start_id The starting ID to assign to generated particles
           *
           * @return A multimap containing cells and their contained particles.
           *
           */
          std::multimap<types::LevelInd, Particle<dim> >
          uniform_random_particles_in_subdomain (const types::particle_index num_particles,
                                                 const types::particle_index start_id);
      };

    }
  }
}

#endif
