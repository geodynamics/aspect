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
       * Euler scheme integrator, where y_{n+1} = y_n + dt * v(y_n).
       * This requires only one step per integration, and doesn't involve any extra data.
       */
      template <int dim>
      class EulerIntegrator : public Interface<dim>
      {
        public:
          virtual bool integrate_step(typename std::multimap<LevelInd, BaseParticle<dim> > &particles,
                                      const std::vector<Tensor<1,dim> > &old_velocities,
                                      const std::vector<Tensor<1,dim> > &,
                                      const double dt);

          virtual unsigned int data_length() const;
          virtual unsigned int read_data(const std::vector<double> &, const unsigned int &pos, const double &);
          virtual void write_data(std::vector<double> &, const double &) const;
      };

    }
  }
}

#endif
