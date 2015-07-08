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

#ifndef __aspect__particle_integrator_rk_2_multistep_h
#define __aspect__particle_integrator_rk_2_multistep_h

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
      * This scheme is identical to the rk2 integrator, but instead of storing the old location
      * in the integrator, it is stored as particle property."
      */
      template <int dim>
      class RK2IntegratorMultiStep : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
      {
        public:
          RK2IntegratorMultiStep();

          virtual bool integrate_step(typename std::multimap<LevelInd, Particle<dim> > &particles,
                                      const std::vector<Tensor<1,dim> > &old_velocities,
                                      const std::vector<Tensor<1,dim> > &velocities,
                                      const double dt);

          virtual unsigned int data_length() const;

          virtual void read_data(std::vector<double>::const_iterator &data,
                                 const double &id_num);

          virtual void write_data(std::vector<double>::iterator &data,
                                  const double &id_num) const;

        private:
          unsigned int                    step;
      };
    }
  }
}

#endif
