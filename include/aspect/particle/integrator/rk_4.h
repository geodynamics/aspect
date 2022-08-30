/*
 Copyright (C) 2015 - 2021 by the authors of the ASPECT code.

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

#ifndef _aspect_particle_integrator_rk_4_h
#define _aspect_particle_integrator_rk_4_h

#include <aspect/particle/integrator/interface.h>

namespace aspect
{
  namespace Particle
  {
    namespace Integrator
    {
      /**
       * Runge Kutta fourth order integrator. This scheme requires
       * storing the original location and intermediate k1, k2, k3
       * values in the particle properties.
       *
       * @ingroup ParticleIntegrators
       */
      template <int dim>
      class RK4 : public Interface<dim>, public SimulatorAccess<dim>
      {
        public:
          RK4();

          /**
           * Look up where the RK4 data is stored. Done once and cached to
           * avoid repeated lookups.
           */
          virtual
          void
          initialize () override;

          /**
           * Perform an integration step of moving the particles of one cell
           * by the specified timestep dt. This class implements a Runge-
           * Kutta integration scheme that is fourth order accurate
           * in space.
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
           * This function is called at the end of every integration step.
           * For the current class it will either increment the step variable
           * and return true to request another integration step, or reset the
           * step variable to 0 and return false if we are already in step 4.
           *
           * @return This function returns true if the integrator requires
           * another integration step. The particle integration will continue
           * to start new integration steps until this function returns false.
           */
          bool new_integration_step() override;

          /**
           * We need to tell the property manager how many intermediate properties this integrator requires,
           * so that it can allocate sufficient space for each particle. However, the integrator is not
           * created at the time the property manager is set up and we can not reverse the order of creation,
           * because the integrator needs to know where to store its properties, which requires the property manager
           * to be finished setting up properties. Luckily the number of properties is constant, so we can make it
           * a static property of this class. Therefore, the property manager can access this variable even
           * before any object is constructed.
           *
           * The Runge-Kutta 4 integrator requires 4 tensors with dim components each.
           */
          static const unsigned int n_integrator_properties = 4*dim;

        private:
          /**
           * The current integration step, i.e. for RK4 a number between 0
           * and 3.
           */
          unsigned int integrator_substep;

          /**
           * The location of the 4 RK4 data fields stored in the particle properties.
           */
          std::array<unsigned int,4> property_index_k;
      };
    }
  }
}

#endif
