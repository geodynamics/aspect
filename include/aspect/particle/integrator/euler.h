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

#ifndef _aspect_particle_integrator_euler_h
#define _aspect_particle_integrator_euler_h

#include <aspect/particle/integrator/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace Particle
  {
    namespace Integrator
    {
      /**
       * Explicit Euler scheme integrator.
       * This requires only one step per integration, and doesn't involve any
       * additional integration steps (it is a one-step method), therefore
       * no data needs to be stored between subsequent steps.
       *
       * @ingroup ParticleIntegrators
       */
      template <int dim>
      class Euler : public Interface<dim>, public SimulatorAccess<dim>
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
          void
          local_integrate_step(const typename ParticleHandler<dim>::particle_iterator &begin_particle,
                               const typename ParticleHandler<dim>::particle_iterator &end_particle,
                               const std::vector<Tensor<1,dim>> &old_velocities,
                               const std::vector<Tensor<1,dim>> &velocities,
                               const double dt) override;

          /**
           * We need to tell the property manager how many intermediate properties this integrator requires,
           * so that it can allocate sufficient space for each particle. However, the integrator is not
           * created at the time the property manager is set up and we can not reverse the order of creation,
           * because the integrator needs to know where to store its properties, which requires the property manager
           * to be finished setting up properties. Luckily the number of properties is constant, so we can make it
           * a static property of this class. Therefore, the property manager can access this variable even
           * before any object is constructed.
           *
           * The forward euler integrator does not need any intermediate storage space.
           */
          static const unsigned int n_integrator_properties = 0;
      };

    }
  }
}

#endif
