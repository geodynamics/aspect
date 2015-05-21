/*
 Copyright (C) 2011 - 2014 by the authors of the ASPECT code.

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

#ifndef __aspect__particle_generator_uniform_box_h
#define __aspect__particle_generator_uniform_box_h

#include <aspect/particle/generator/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace Particle
  {
    namespace Generator
    {
      // Generate random uniform distribution of particles over entire simulation domain
      template <int dim>
      class UniformBox : public Interface<dim>, public SimulatorAccess<dim>
      {
        public:
          /**
           * Constructor.
           *
           * @param[in] The MPI communicator for synchronizing particle generation.
           */
          UniformBox();

          /**
           * Generate a uniformly randomly distributed set of particles in the current triangulation.
           */
          // TODO: fix the particle system so it works even with processors assigned 0 cells
          virtual
          void
          generate_particles(Particle::World<dim> &world,
                             const double total_num_particles);

          /**
           * TODO: Update comments
           * Generate a set of particles uniformly distributed within the
           * specified triangulation.
           * selection weighted by cell volume. We do cell-by-cell assignment of
           * particles because the decomposition of the mesh may result in a highly
           * non-rectangular local mesh which makes uniform particle distribution difficult.
           *
           * @param [in] world The particle world the particles will exist in
           * @param [in] num_particles The number of particles to generate in this subdomain
           * @param [in] start_id The starting ID to assign to generated particles
           */
          void uniformly_distributed_particles_in_subdomain (Particle::World<dim> &world,
                                                      const unsigned int num_particles,
                                                      const unsigned int start_id);

          void
          generate_particle(Particle::World<dim> &world,
                            const Point<dim> &position,
                            const unsigned int id);

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
          Point<dim> P_min;
          Point<dim> P_max;
      };


    }
  }
}

#endif
