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

#ifndef __aspect__particle_integrator_euler_h
#define __aspect__particle_integrator_euler_h

#include <aspect/particle/integrator/interface.h>

namespace aspect
{
  namespace Particle
  {
    namespace Integrator
    {
      /**
       * Explicit Euler scheme integrator, where y_{n+1} = y_n + dt * v(y_n).
       * This requires only one step per integration, and doesn't involve any
       * additional integration steps (it is a one-step method), therefore
       * no data needs to be stored between subsequent steps.
       *
       * @ingroup ParticleIntegrators
       */
      template <int dim>
      class Euler : public Interface<dim>
      {
        public:
          /**
           * Perform an integration step of moving the particles of one cell
           * by the specified timestep dt. This function implements an explicit
           * euler integration scheme.
           *
           * @param [in] begin_particle An iterator to the first particle to be moved.
           * @param [in] end_particle An iterator to the last particle to be moved.
           * @param [in] old_velocities The velocities at t_n, i.e. before the
           * particle movement, for all particles between @p begin_particle
           * and @p end_particle at their current position.
           * @param [in] velocities The velocities at the particle positions
           * at t_{n+1}, i.e. after the particle movement. Note that this is
           * the velocity at the old positions, but at the new time. It is the
           * responsibility of this function to compute the new location of
           * the particles.
           * @param [in] dt The length of the integration timestep.
           */
          virtual
          void
          local_integrate_step(const typename std::multimap<types::LevelInd, Particle<dim> >::iterator &begin_particle,
                               const typename std::multimap<types::LevelInd, Particle<dim> >::iterator &end_particle,
                               const std::vector<Tensor<1,dim> > &old_velocities,
                               const std::vector<Tensor<1,dim> > &velocities,
                               const double dt);
      };

    }
  }
}

#endif
