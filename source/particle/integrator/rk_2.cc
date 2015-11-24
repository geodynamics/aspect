/*
  Copyright (C) 2015 by the authors of the ASPECT code.

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
      template <int dim>
      RK2<dim>::RK2()
      {
        step = 0;
        loc0.clear();
      }

      template <int dim>
      void
      RK2<dim>::local_integrate_step(const typename std::multimap<types::LevelInd, Particle<dim> >::iterator &begin_particle,
                                     const typename std::multimap<types::LevelInd, Particle<dim> >::iterator &end_particle,
                                     const std::vector<Tensor<1,dim> > &old_velocities,
                                     const std::vector<Tensor<1,dim> > &velocities,
                                     const double dt)
      {
        typename std::multimap<types::LevelInd, Particle<dim> >::iterator it = begin_particle;
        typename std::vector<Tensor<1,dim> >::const_iterator old_vel = old_velocities.begin();
        typename std::vector<Tensor<1,dim> >::const_iterator vel = velocities.begin();

        for (; it!=end_particle, vel!=velocities.end(), old_vel!=old_velocities.end(); ++it,++vel,++old_vel)
          {
            const types::particle_index id_num = it->second.get_id();
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
      }

      template <int dim>
      bool
      RK2<dim>::new_integration_step()
      {
        if (step == 1) loc0.clear();
        step = (step+1)%2;

        // Continue until we're at the last step
        return (step != 0);
      }

      template <int dim>
      unsigned int
      RK2<dim>::get_data_size() const
      {
        // TODO: If we finished integration we do not need to transfer integrator data. Return 0 in that case.
        return dim * sizeof(double);
      }

      template <int dim>
      void
      RK2<dim>::read_data(const void *&data,
                          const types::particle_index id_num)
      {
        // TODO: If we finished integration we do not need to transfer integrator data. Return early in that case.

        const double *integrator_data = static_cast<const double *> (data);

        // Read location data
        for (unsigned int i=0; i<dim; ++i)
          loc0[id_num](i) = *integrator_data++;

        data = static_cast<const void *> (integrator_data);
      }

      template <int dim>
      void
      RK2<dim>::write_data(void *&data,
                           const types::particle_index id_num) const
      {
        // TODO: If we finished integration we do not need to transfer integrator data. Return early in that case.

        double *integrator_data = static_cast<double *> (data);

        // Write location data
        const typename std::map<types::particle_index, Point<dim> >::const_iterator it = loc0.find(id_num);
        for (unsigned int i=0; i<dim; ++i,++integrator_data)
          *integrator_data = it->second(i);

        data = static_cast<void *> (integrator_data);
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
      ASPECT_REGISTER_PARTICLE_INTEGRATOR(RK2,
                                          "rk2",
                                          "Runge Kutta second order integrator, where "
                                          "y_{n+1} = y_n + dt*v(0.5*k_1), k_1 = dt*v(y_n). "
                                          "This scheme requires storing the original location, "
                                          "and the read/write_data functions reflect this.")
    }
  }
}

