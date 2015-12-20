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
        :
        integrator_substep(0)
      {}

      template <int dim>
      void
      RK2<dim>::local_integrate_step(const typename std::multimap<types::LevelInd, Particle<dim> >::iterator &begin_particle,
                                     const typename std::multimap<types::LevelInd, Particle<dim> >::iterator &end_particle,
                                     const std::vector<Tensor<1,dim> > &old_velocities,
                                     const std::vector<Tensor<1,dim> > &velocities,
                                     const double dt)
      {
        Assert(static_cast<unsigned int> (std::distance(begin_particle, end_particle)) == old_velocities.size(),
               ExcMessage("The particle integrator expects the old velocity vector to be of equal size "
                          "to the number of particles to advect. For some unknown reason they are different, "
                          "most likely something went wrong in the calling function."));

        Assert(old_velocities.size() == velocities.size(),
               ExcMessage("The particle integrator expects the velocity vector to be of equal size "
                          "to the number of particles to advect. For some unknown reason they are different, "
                          "most likely something went wrong in the calling function."));

        typename std::vector<Tensor<1,dim> >::const_iterator old_velocity = old_velocities.begin();
        typename std::vector<Tensor<1,dim> >::const_iterator velocity = velocities.begin();

        for (typename std::multimap<types::LevelInd, Particle<dim> >::iterator it = begin_particle;
             it != end_particle; ++it, ++velocity, ++old_velocity)
          {
            const types::particle_index particle_id = it->second.get_id();
            const Point<dim> loc = it->second.get_location();
            if (integrator_substep == 0)
              {
                loc0[particle_id] = loc;
                it->second.set_location(loc + 0.5 * dt * (*old_velocity));
              }
            else if (integrator_substep == 1)
              {
                it->second.set_location(loc0[particle_id] + dt * (*old_velocity + *velocity) / 2.0);
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
        if (integrator_substep == 1) loc0.clear();
        integrator_substep = (integrator_substep + 1) % 2;

        // Continue until we're at the last step
        return (integrator_substep != 0);
      }

      template <int dim>
      unsigned int
      RK2<dim>::get_data_size() const
      {
        // If integration is finished, we do not need to transfer integrator
        // data to other processors, because it will be deleted soon anyway.
        // Skip the MPI transfer in this case.
        if (integrator_substep == 1)
          return 0;

        return dim * sizeof(double);
      }

      template <int dim>
      const void *
      RK2<dim>::read_data(const void *data,
                          const types::particle_index particle_id)
      {
        // If integration is finished, we do not need to transfer integrator
        // data to other processors, because it will be deleted soon anyway.
        // Skip the MPI transfer in this case.
        if (integrator_substep == 1)
          return data;

        const double *integrator_data = static_cast<const double *> (data);

        // Read location data
        for (unsigned int i=0; i<dim; ++i)
          loc0[particle_id](i) = *integrator_data++;

        return static_cast<const void *> (integrator_data);
      }

      template <int dim>
      void *
      RK2<dim>::write_data(void *data,
                           const types::particle_index particle_id) const
      {
        // If integration is finished, we do not need to transfer integrator
        // data to other processors, because it will be deleted soon anyway.
        // Skip the MPI transfer in this case.
        if (integrator_substep == 1)
          return data;

        double *integrator_data = static_cast<double *> (data);

        // Write location data
        const typename std::map<types::particle_index, Point<dim> >::const_iterator it = loc0.find(particle_id);
        for (unsigned int i=0; i<dim; ++i,++integrator_data)
          *integrator_data = it->second(i);

        return static_cast<void *> (integrator_data);
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
                                          "$y_{n+1} = y_n + dt*v(0.5*k_1)$, $k_1 = dt*v(y_n)$.")
    }
  }
}

