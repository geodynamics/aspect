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

#ifndef __aspect__particle_integrator_rk_2_h
#define __aspect__particle_integrator_rk_2_h

#include <aspect/particle/integrator/interface.h>

#include <aspect/simulator_access.h>


namespace aspect
{
  namespace Particle
  {
    namespace Integrator
    {
      /**
       * Runge Kutta second order integrator, where y_{n+1} = y_n + dt*v(0.5*k_1), k_1 = dt*v(y_n).
       * This scheme requires storing the original location, and the read/write_data functions reflect this.
       *
       * @ingroup ParticleIntegrators
       */
      template <int dim>
      class RK2 : public Interface<dim>, public SimulatorAccess<dim>
      {
        public:
          RK2();

          /**
           * Perform an integration step of moving the particles of one cell
           * by the specified timestep dt. This function implements a
           * second-order accurate Runge-Kutta integration scheme.
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

          /**
           * This function is called at the end of every integration step.
           * For the current class it will either increment the step variable
           * and return true to request another integration step, or reset the
           * step variable to 0 and return false if we are already in step 2.
           *
           * @return This function returns true if the integrator requires
           * another integration step. The particle integration will continue
           * to start new integration steps until this function returns false.
           */
          virtual bool new_integration_step();

          /**
           * Return data length of the integration related data required for
           * communication in terms of number of bytes. When data about
           * particles is transported from one processor to another, or stored
           * on disk for snapshots, integrators get the chance to store
           * whatever information they need with each particle. This function
           * returns how many pieces of additional information a concrete
           * integrator class needs to store for each particle.
           *
           * @return The number of bytes required to store the relevant
           * integrator data for one particle.
           */
          virtual unsigned int get_data_size() const;

          /**
           * Read integration related data for a particle specified by particle_id
           * from the data array. This function is called after transferring
           * a particle to the local domain during an integration step.
           *
           * @param [in] data A pointer into the data array. The pointer
           * marks the position where this function starts reading.
           * @param [in] particle_id The id number of the particle to read the data
           * for.
           * @return The updated position of the pointer into the data array.
           * The return value is @p data advanced by get_data_size() bytes.
           */
          virtual
          const void *
          read_data(const void *data,
                    const types::particle_index particle_id);

          /**
           * Write integration related data to a vector for a particle
           * specified by particle_id.
           *
           * @param [in,out] data A pointer to the array of data to write
           * integrator data into. The pointer should be advanced by
           * get_data_size() bytes within this function.
           * @param [in] particle_id The id number of the particle to write the data
           * for.
           */
          virtual
          void *
          write_data(void *data,
                     const types::particle_index particle_id) const;

        private:
          /**
           * The current integration step, i.e for RK2 a number that is either
           * 0 or 1.
           */
          unsigned int integrator_substep;

          /**
           * The particle location before the first integration step. This is
           * used in the second step and transferred to another process if
           * the tracer leaves the domain during the first step.
           */
          std::map<types::particle_index, Point<dim> >   loc0;

      };

    }
  }
}

#endif
