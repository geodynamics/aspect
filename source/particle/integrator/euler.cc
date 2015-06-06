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

#include <aspect/particle/integrator/euler.h>

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
      bool
      EulerIntegrator<dim>::integrate_step(typename std::multimap<LevelInd, Particle<dim> > &particles,
                                           const std::vector<Tensor<1,dim> > &old_velocities,
                                           const std::vector<Tensor<1,dim> > &,
                                           const double dt)
      {
        typename std::multimap<LevelInd, Particle<dim> >::iterator it = particles.begin();
        typename std::vector<Tensor<1,dim> >::const_iterator vel = old_velocities.begin();

        for (; it!=particles.end(), vel!=old_velocities.end(); ++it,++vel)
          {
            const Point<dim> loc = it->second.get_location();
            it->second.set_location(loc + dt*(*vel));
          }

        return false;
      }

      template <int dim>
      unsigned int
      EulerIntegrator<dim>::data_length() const
      {
        return 0;
      }
      template <int dim>
      unsigned int
      EulerIntegrator<dim>::read_data(const std::vector<double> &,
                                      const unsigned int &pos,
                                      const double &)
      {
        return pos;
      }
      template <int dim>
      void
      EulerIntegrator<dim>::write_data(std::vector<double> &,
                                       const double &) const
      {
      }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Particle
  {
    namespace Integrator
    {
      ASPECT_REGISTER_PARTICLE_INTEGRATOR(EulerIntegrator,
                                          "euler",
                                          "Euler scheme integrator, where y_{n+1} = y_n + dt * v(y_n). "
                                          "This requires only one step per integration, and doesn't "
                                          "involve any extra data.")
    }
  }
}
