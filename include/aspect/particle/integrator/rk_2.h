/*
 Copyright (C) 2015 - 2024 by the authors of the ASPECT code.

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

#ifndef _aspect_particle_integrator_rk_2_h
#define _aspect_particle_integrator_rk_2_h

#include <aspect/particle/integrator/interface.h>

#include <aspect/simulator_access.h>


namespace aspect
{
  namespace Particle
  {
    namespace Integrator
    {
      /**
       * Runge Kutta second order integrator.
       * This scheme requires storing the original location in the particle properties.
       *
       * @ingroup ParticleIntegrators
       */
      template <int dim>
      class RK2 : public Interface<dim>, public SimulatorAccess<dim>
      {
        public:
          RK2();

          /**
           * Look up where the RK2 data is stored. Done once and cached to
           * avoid repeated lookups.
           */
          void
          initialize () override;

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
           * step variable to 0 and return false if we are already in step 2.
           *
           * @return This function returns true if the integrator requires
           * another integration step. The particle integration will continue
           * to start new integration steps until this function returns false.
           */
          bool new_integration_step() override;

          /**
           * Return a list of boolean values indicating which solution vectors
           * are required for the integration. The first entry indicates if
           * the particle integrator requires the solution vector at the old
           * old time (k-1), the second entry indicates if the particle integrator
           * requires the solution vector at the old time (k), and the third entry
           * indicates if the particle integrator requires the solution vector
           * at the new time (k+1).
           *
           * The RK2 integrator requires the solution vector at the
           * old time (k) for the first integration step, and the solution
           * vector at both the old and new time for the second integration step
           * (if higher_order_in_time is set to true).
           */
          std::array<bool, 3> required_solution_vectors() const override;

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

          /**
           * We need to tell the property manager how many intermediate properties this integrator requires,
           * so that it can allocate sufficient space for each particle. However, the integrator is not
           * created at the time the property manager is set up and we can not reverse the order of creation,
           * because the integrator needs to know where to store its properties, which requires the property manager
           * to be finished setting up properties. Luckily the number of properties is constant, so we can make it
           * a static property of this class. Therefore, the property manager can access this variable even
           * before any object is constructed.
           *
           * The Runge-Kutta 2 integrator requires a single point with dim components.
           */
          static constexpr unsigned int n_integrator_properties = dim;

        private:
          /**
           * The current integration step, i.e for RK2 a number that is either
           * 0 or 1.
           */
          unsigned int integrator_substep;

          /**
           * The location of the RK2 data that is stored in the particle properties.
           */
          unsigned int property_index_old_location;

          /**
           * Whether to evaluate old and current velocity to compute a solution
           * that is higher order accurate in time.
           */
          bool higher_order_in_time;
      };

    }
  }
}

#endif
