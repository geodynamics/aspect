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

#include <aspect/particle/integrator/rk_4.h>

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
          RK4Integrator<dim>::RK4Integrator(void)
          {
            step = 0;
            loc0.clear();
            k1.clear();
            k2.clear();
            k3.clear();
          }

        template <int dim>
          bool
          RK4Integrator<dim>::integrate_step(typename std::multimap<LevelInd, BaseParticle<dim> > &particles,
                                             const double                                          dt)
          {
            typename std::multimap<LevelInd, BaseParticle<dim> >::iterator       it;
            Point<dim>                          loc, vel, k4;
            double                              id_num;

            for (it=particles.begin(); it!=particles.end(); ++it)
              {
                id_num = it->second.get_id();
                loc = it->second.get_location();
                vel = it->second.get_velocity();
                if (step == 0)
                  {
                    loc0[id_num] = loc;
                    k1[id_num] = dt*vel;
                    it->second.set_location(loc + 0.5*k1[id_num]);
                  }
                else if (step == 1)
                  {
                    k2[id_num] = dt*vel;
                    it->second.set_location(loc0[id_num] + 0.5*k2[id_num]);
                  }
                else if (step == 2)
                  {
                    k3[id_num] = dt*vel;
                    it->second.set_location(loc0[id_num] + k3[id_num]);
                  }
                else if (step == 3)
                  {
                    k4 = dt*vel;
                    it->second.set_location(loc0[id_num] + (k1[id_num]+2.0*k2[id_num]+2.0*k3[id_num]+k4)/6.0);
                  }
                else
                  {
                    // Error!
                  }
              }

            step = (step+1)%4;
            if (step == 0)
              {
                loc0.clear();
                k1.clear();
                k2.clear();
                k3.clear();
              }

            // Continue until we're at the last step
            return (step != 0);
          }

        template <int dim>
          void
          RK4Integrator<dim>::add_mpi_types(std::vector<MPIDataInfo> &data_info)
          {
            // Add the loc0, k1, k2, and k3 data
            data_info.push_back(MPIDataInfo("loc0", dim));
            data_info.push_back(MPIDataInfo("k1", dim));
            data_info.push_back(MPIDataInfo("k2", dim));
            data_info.push_back(MPIDataInfo("k3", dim));
          }

        template <int dim>
          unsigned int
          RK4Integrator<dim>::data_len() const
          {
            return 4*dim;
          }

        template <int dim>
          unsigned int
          RK4Integrator<dim>::read_data(const std::vector<double> &data, const unsigned int &pos, const double &id_num)
          {
            unsigned int    i, p = pos;

            // Read location data
            for (i=0; i<dim; ++i)
              {
                loc0[id_num](i) = data[p++];
              }
            // Read k1, k2 and k3
            for (i=0; i<dim; ++i)
              {
                k1[id_num](i) = data[p++];
              }
            for (i=0; i<dim; ++i)
              {
                k2[id_num](i) = data[p++];
              }
            for (i=0; i<dim; ++i)
              {
                k3[id_num](i) = data[p++];
              }

            return p;
          }

        template <int dim>
          void
          RK4Integrator<dim>::write_data(std::vector<double> &data, const double &id_num) const
          {
            typename std::map<double, Point<dim> >::const_iterator it;
            unsigned int  i;

            // Write location data
            it = loc0.find(id_num);
            for (i=0; i<dim; ++i)
              {
                data.push_back(it->second(i));
              }
            // Write k1, k2 and k3
            it = k1.find(id_num);
            for (i=0; i<dim; ++i)
              {
                data.push_back(it->second(i));
              }
            it = k2.find(id_num);
            for (i=0; i<dim; ++i)
              {
                data.push_back(it->second(i));
              }
            it = k3.find(id_num);
            for (i=0; i<dim; ++i)
              {
                data.push_back(it->second(i));
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
    ASPECT_REGISTER_PARTICLE_INTEGRATOR(RK4Integrator,
                                               "rk4",
                                               "Runge Kutta fourth order integrator, where "
                                               "y_{n+1} = y_n + (1/6)*k1 + (1/3)*k2 + (1/3)*k3 + (1/6)*k4 "
                                               "and k1, k2, k3, k4 are defined as usual. "
                                               "This scheme requires storing the original location "
                                               "and intermediate k1, k2, k3 values, so the "
                                               "read/write_data functions reflect this.")
    }
  }
}
