/*
 Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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

#ifndef __aspect__particle_generator_random_function_h
#define __aspect__particle_generator_random_function_h

#include <aspect/particle/generator/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/particle/definitions.h>

#include <deal.II/base/parsed_function.h>

#include <boost/random.hpp>

namespace aspect
{
  namespace Particle
  {
    namespace Generator
    {
      /**
       * Generate random distribution of particles over entire simulation domain.
       * The density of particles is determined by a user-defined probability
       * density function.
       */
      template <int dim>
      class RandomFunction : public Interface<dim>, public SimulatorAccess<dim>
      {
        public:
          /**
           * Constructor.
           *
           */
          RandomFunction();

          /**
           * Generate a uniformly randomly distributed set of particles in the current triangulation.
           */
          virtual
          void
          generate_particles(const double total_num_particles,
                             Particle::World<dim> &world);

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
           * A function object representing the temperature.
           */
          Functions::ParsedFunction<dim> function;

          /**
           * Generate a set of particles uniformly randomly distributed within the
           * specified triangulation. This is done using "roulette wheel" style
           * selection weighted by cell volume. We do cell-by-cell assignment of
           * particles because the decomposition of the mesh may result in a highly
           * non-rectangular local mesh which makes uniform particle distribution difficult.
           *
           * @param [in] cells Map between accumulated cell weight and cell index
           * @param [in] global_weight The integrated probability function
           * @param [in] start_weight The starting weight of the first cell of the local process.
           * @param [in] num_particles The total number of particles to generate.
           * @param [in] start_id The starting ID to assign to generated particles of the local process.
           * @param [inout] world The particle world the particles will exist in
           *
           */
          void uniform_random_particles_in_subdomain (const std::map<double,LevelInd> &cells,
                                                      const double global_weight,
                                                      const double start_weight,
                                                      const unsigned int num_particles,
                                                      const unsigned int start_id,
                                                      Particle::World<dim> &world);
      };

    }
  }
}

#endif
