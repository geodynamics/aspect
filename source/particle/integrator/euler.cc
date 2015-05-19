/*
  Copyright (C) 2011 - 2014 by the authors of the ASPECT code.

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
          EulerIntegrator<dim>::integrate_step(typename std::multimap<LevelInd, BaseParticle<dim> > &particles,
                                               const double dt)
          {
            typename std::multimap<LevelInd, BaseParticle<dim> >::iterator       it;
            Point<dim>                          loc, vel;

            for (it=particles.begin(); it!=particles.end(); ++it)
              {
                loc = it->second.get_location();
                vel = it->second.get_velocity();
                it->second.set_location(loc + dt*vel);
              }

            return false;
          }

        template <int dim>
          void
          EulerIntegrator<dim>::add_mpi_types(std::vector<MPIDataInfo> &data_info)
          {}
        template <int dim>
          unsigned int
          EulerIntegrator<dim>::data_len() const
          {
            return 0;
          }
        template <int dim>
          unsigned int
          EulerIntegrator<dim>::read_data(const std::vector<double> &data, const unsigned int &pos, const double &id_num)
          {
            return pos;
          }
        template <int dim>
          void
          EulerIntegrator<dim>::write_data(std::vector<double> &data, const double &id_num) const
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
