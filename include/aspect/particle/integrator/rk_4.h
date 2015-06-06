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
          virtual bool integrate_step(typename std::multimap<LevelInd, Particle<dim> > &particles,
                                      const std::vector<Tensor<1,dim> > &old_velocities,
                                      const std::vector<Tensor<1,dim> > &velocities,
                                      const double dt);

          virtual unsigned int data_length() const;
          virtual unsigned int read_data(const std::vector<double> &data, const unsigned int &pos, const double &id_num);
          virtual void write_data(std::vector<double> &data, const double &id_num) const;

        private:
          unsigned int                     step;
          std::map<double, Point<dim> >    loc0;
          std::map<double, Tensor<1,dim> > k1, k2, k3;

      };
    }
  }
}

#endif
