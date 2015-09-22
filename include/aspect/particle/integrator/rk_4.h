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

#ifndef __aspect__particle_integrator_rk_4_h
#define __aspect__particle_integrator_rk_4_h

#include <aspect/particle/integrator/interface.h>

namespace aspect
{
  namespace Particle
  {
    namespace Integrator
    {
      /**
      * Runge Kutta fourth order integrator, where y_{n+1} = y_n + (1/6)*k1 + (1/3)*k2 + (1/3)*k3 + (1/6)*k4
      * and k1, k2, k3, k4 are defined as usual.
      * This scheme requires storing the original location and intermediate k1, k2, k3 values,
      * so the read/write_data functions reflect this.
      */
      template <int dim>
      class RK4Integrator : public Interface<dim>
      {
        public:
          RK4Integrator();

          /**
           * Perform an integration step of moving the particles of one cell
           * by the specified timestep dt. Implementations of this function
           * must update the particle location. Between calls to this function
           * the velocity at the updated particle positions is evaluated and
           * passed to integrate_step during the next call.
           *
           * @param [in] begin_particle An iterator to the first particle to be moved.
           * @param [in] end_particle An iterator to the last particle to be moved.
           * @param [in] old_velocities The velocities at t_n, i.e. before the
           * particle movement.
           * @param [in] velocities The velocities at t_{n+1}, i.e. after the
           * particle movement.
           * @param [in] dt The length of the integration timestep.
           */
          virtual
          void
          local_integrate_step(const typename std::multimap<LevelInd, Particle<dim> >::iterator &begin_particle,
                               const typename std::multimap<LevelInd, Particle<dim> >::iterator &end_particle,
                               const std::vector<Tensor<1,dim> > &old_velocities,
                               const std::vector<Tensor<1,dim> > &velocities,
                               const double dt);

          /**
           * This function is called at the end of every integration step.
           * In case of multi-step integrators is signals the beginning of a
           * new integration step.
           */
          virtual void advance_step();

          /**
           * @return This function returns true if another integration step is required
           * by the integrator. Its default implementation returns false.
           */
          virtual bool continue_integration() const;

          /**
           * Return data length of the integration related data required for
           * communication in terms of number of doubles.
           *
           * @return The number of doubles required to store the relevant
           * integrator data for one particle.
           */
          virtual unsigned int data_length() const;

          /**
           * Read integration related data for a particle specified by id_num
           * from the data array.
           *
           * @param [in] data The array of double data to read from.
           * @param [in] id_num The id number of the particle to read the data
           * for.
           */
          virtual void read_data(const void *&data,
                                 const particle_index id_num);

          /**
           * Write integration related data to a vector for a particle
           * specified by id_num.
           *
           * @param [in,out] data The vector of doubles to write integrator
           * data into.
           * @param [in] id_num The id number of the particle to write the data
           * for.
           */
          virtual void write_data(void *&data,
                                  const particle_index id_num) const;

        private:
          /**
           * The current integration step.
           */
          unsigned int                     step;

          /**
           * The particle location before the first integration step. This is
           * used in the following steps and transferred to another process if
           * the tracer leaves the domain during one of the steps.
           */
          std::map<particle_index, Point<dim> >    loc0;

          /**
           * The intermediate values of the RK4 scheme. These are
           * used in the following steps and transferred to another process if
           * the tracer leaves the domain during one of the steps.
           */
          std::map<particle_index, Tensor<1,dim> > k1, k2, k3;

      };
    }
  }
}

#endif
