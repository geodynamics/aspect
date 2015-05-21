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

namespace aspect
{
  namespace Particle
  {
    namespace Integrator
    {
      /**
      * Runge Kutta second order integrator, where y_{n+1} = y_n + dt*v(0.5*k_1), k_1 = dt*v(y_n).
      * This scheme requires storing the original location, and the read/write_data functions reflec
      * This version is different in that particle movement "lags" behind by one timestep, and there
      */
      template <int dim>
      class RK2IntegratorMultiStep : public Interface<dim>
      {
        public:
          RK2Integrator();
          virtual bool integrate_step(typename std::multimap<LevelInd, BaseParticle<dim> > &particles, const double dt);
          virtual void add_mpi_types(std::vector<MPIDataInfo> &data_info);
          virtual unsigned int data_len() const;
          virtual unsigned int read_data(const std::vector<double> &data, const unsigned int &pos, const double &id_num);
          virtual void write_data(std::vector<double> &data, const double &id_num) const;

        private:
          unsigned int                    step;
          std::map<double, Point<dim> >   loc0;

      };
    }
  }
}

#endif
