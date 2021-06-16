/*
  Copyright (C) 2015 - 2018 by the authors of the ASPECT code.

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

#include <aspect/particle/integrator/rk_4.h>

namespace aspect
{
  namespace Particle
  {
    namespace Integrator
    {
      template <int dim>
      RK4<dim>::RK4()
        :
        integrator_substep(0)
      {}



      template <int dim>
      void
      RK4<dim>::initialize ()
      {
        const auto &property_information = this->get_particle_world().get_property_manager().get_data_info();

        property_k[0] = property_information.get_position_by_field_name("integrator properties");
        property_k[1] = property_k[0] + dim;
        property_k[2] = property_k[1] + dim;
        property_k[3] = property_k[2] + dim;
      }



      template <int dim>
      void
      RK4<dim>::local_integrate_step(const typename ParticleHandler<dim>::particle_iterator &begin_particle,
                                     const typename ParticleHandler<dim>::particle_iterator &end_particle,
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

        for (typename ParticleHandler<dim>::particle_iterator it = begin_particle;
             it != end_particle; ++it, ++velocity, ++old_velocity)
          {
            ArrayView<double> properties = it->get_properties();

            if (integrator_substep == 0)
              {
                const Tensor<1,dim> k1 = dt * (*old_velocity);
                const Point<dim> loc0 = it->get_location();

                for (unsigned int i=0; i<dim; ++i)
                  {
                    properties[property_k[0] + i] = loc0[i];
                    properties[property_k[1] + i] = k1[i];
                  }

                it->set_location(loc0 + 0.5 * k1);
              }
            else if (integrator_substep == 1)
              {
                const Tensor<1,dim> k2 = dt * ((*old_velocity) + (*velocity)) / 2.0;
                Point<dim> loc0;

                for (unsigned int i=0; i<dim; ++i)
                  {
                    loc0[i] = properties[property_k[0] + i];
                    properties[property_k[2] + i] = k2[i];
                  }
                it->set_location(loc0 + 0.5 * k2);
              }
            else if (integrator_substep == 2)
              {
                const Tensor<1,dim> k3 = dt * (*old_velocity + *velocity) / 2.0;
                Point<dim> loc0;

                for (unsigned int i=0; i<dim; ++i)
                  {
                    loc0[i] = properties[property_k[0] + i];
                    properties[property_k[3] + i] = k3[i];
                  }

                it->set_location(loc0 + k3);
              }
            else if (integrator_substep == 3)
              {
                const Tensor<1,dim> k4 = dt * (*velocity);
                Point<dim> loc0, k1, k2, k3;

                for (unsigned int i=0; i<dim; ++i)
                  {
                    loc0[i] = properties[property_k[0] + i];
                    k1[i] = properties[property_k[1] + i];
                    k2[i] = properties[property_k[2] + i];
                    k3[i] = properties[property_k[3] + i];
                  }

                it->set_location(loc0 + (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0);
              }
            else
              {
                Assert(false,
                       ExcMessage("The RK4 integrator should never continue after four integration steps."));
              }
          }
      }



      template <int dim>
      bool
      RK4<dim>::new_integration_step()
      {
        integrator_substep = (integrator_substep+1)%4;

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
      ASPECT_REGISTER_PARTICLE_INTEGRATOR(RK4,
                                          "rk4",
                                          "Runge Kutta fourth order integrator, where "
                                          "$y_{n+1} = y_n + \\frac{1}{6} k_1 + \\frac{1}{3} k_2 "
                                          "+ \\frac{1}{3} k_3 + \\frac{1}{6} k_4$ "
                                          "and $k_1$, $k_2$, $k_3$, $k_4$ are defined as usual.")
    }
  }
}
