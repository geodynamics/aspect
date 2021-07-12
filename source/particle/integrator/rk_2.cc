/*
  Copyright (C) 2015 - 2021 by the authors of the ASPECT code.

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
 along with ASPECT; see the file LICENSE.  If not see
 <http://www.gnu.org/licenses/>.
 */

#include <aspect/particle/integrator/rk_2.h>
#include <aspect/particle/property/interface.h>
#include <aspect/particle/world.h>

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
      RK2<dim>::initialize ()
      {
        const auto &property_information = this->get_particle_world().get_property_manager().get_data_info();
        property_index_old_location = property_information.get_position_by_field_name("internal: integrator properties");
      }



      template <int dim>
      void
      RK2<dim>::local_integrate_step(const typename ParticleHandler<dim>::particle_iterator &begin_particle,
                                     const typename ParticleHandler<dim>::particle_iterator &end_particle,
                                     const std::vector<Tensor<1,dim>> &old_velocities,
                                     const std::vector<Tensor<1,dim>> &velocities,
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

        typename std::vector<Tensor<1,dim>>::const_iterator old_velocity = old_velocities.begin();
        typename std::vector<Tensor<1,dim>>::const_iterator velocity = velocities.begin();

        for (typename ParticleHandler<dim>::particle_iterator it = begin_particle;
             it != end_particle; ++it, ++velocity, ++old_velocity)
          {
            ArrayView<double> properties = it->get_properties();

            if (integrator_substep == 0)
              {
                const Tensor<1,dim> k1 = dt * (*old_velocity);
                const Point<dim> loc0 = it->get_location();

                for (unsigned int i=0; i<dim; ++i)
                  properties[property_index_old_location + i] = loc0[i];

                it->set_location(loc0 + 0.5 * k1);
              }
            else if (integrator_substep == 1)
              {
                const Tensor<1,dim> k2 = dt * (*old_velocity + *velocity) / 2.0;
                Point<dim> loc0;

                for (unsigned int i=0; i<dim; ++i)
                  loc0[i] = properties[property_index_old_location + i];

                it->set_location(loc0 + k2);
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
        integrator_substep = (integrator_substep + 1) % 2;

        // Continue until we're at the last step
        return (integrator_substep != 0);
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
                                          "Second Order Runge Kutta integrator "
                                          "$y_{n+1} = y_n + \\Delta t\\, v(t_{n+1/2}, y_{n} + \\frac{1}{2} k_1)$ "
                                          "where $k_1 = \\Delta t\\, v(t_{n}, y_{n})$")
    }
  }
}
