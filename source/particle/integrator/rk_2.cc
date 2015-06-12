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

#include <aspect/particle/integrator/rk_2.h>

namespace aspect
{
  namespace Particle
  {
    namespace Integrator
    {
      /**
       * Runge Kutta second order integrator, where y_{n+1} = y_n + dt*v(0.5*k_1), k_1 = dt*v(y_n).
       * This scheme requires storing the original location, and the read/write_data functions reflect this.
       */
      template <int dim>
      RK2Integrator<dim>::RK2Integrator(void)
      {
        step = 0;
        loc0.clear();
      }


      template <int dim>
      bool
      RK2Integrator<dim>::integrate_step(typename std::multimap<LevelInd, Particle<dim> > &particles,
                                         const std::vector<Tensor<1,dim> > &old_velocities,
                                         const std::vector<Tensor<1,dim> > &velocities,
                                         const double dt)
      {
        typename std::multimap<LevelInd, Particle<dim> >::iterator it = particles.begin();
        typename std::vector<Tensor<1,dim> >::const_iterator old_vel = old_velocities.begin();
        typename std::vector<Tensor<1,dim> >::const_iterator vel = velocities.begin();

        for (; it!=particles.end(), vel!=velocities.end(), old_vel!=old_velocities.end(); ++it,++vel,++old_vel)
          {
            const double id_num = it->second.get_id();
            const Point<dim> loc = it->second.get_location();
            if (step == 0)
              {
                loc0[id_num] = loc;
                it->second.set_location(loc + 0.5 * dt * (*old_vel));
              }
            else if (step == 1)
              {
                it->second.set_location(loc0[id_num] + dt * (*old_vel + *vel)/2.0);
              }
            else
              {
                Assert(false,
                       ExcMessage("The RK2 integrator should never continue after two integration steps."));
              }
          }

        if (step == 1) loc0.clear();
        step = (step+1)%2;

        // Continue until we're at the last step
        return (step != 0);
      }

      template <int dim>
      unsigned int
      RK2Integrator<dim>::data_length() const
      {
        return dim;
      }

      template <int dim>
      void
      RK2Integrator<dim>::read_data(std::vector<double>::const_iterator &data,
                                    const double &id_num)
      {
        // Read location data
        for (unsigned int i=0; i<dim; ++i)
          {
            loc0[id_num](i) = *data++;
          }
      }

      template <int dim>
      void
      RK2Integrator<dim>::write_data(std::vector<double>::iterator &data, const double &id_num) const
      {
        // Write location data
        const typename std::map<double, Point<dim> >::const_iterator it = loc0.find(id_num);
        for (unsigned int i=0; i<dim; ++i,++data)
          {
            *data = it->second(i);
          }
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
      ASPECT_REGISTER_PARTICLE_INTEGRATOR(RK2Integrator,
                                          "rk2",
                                          "Runge Kutta second order integrator, where "
                                          "y_{n+1} = y_n + dt*v(0.5*k_1), k_1 = dt*v(y_n). "
                                          "This scheme requires storing the original location, "
                                          "and the read/write_data functions reflect this.")
    }
  }
}

