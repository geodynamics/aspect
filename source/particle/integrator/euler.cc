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
      void
      Euler<dim>::local_integrate_step(const typename std::multimap<types::LevelInd, Particle<dim> >::iterator &begin_particle,
                                       const typename std::multimap<types::LevelInd, Particle<dim> >::iterator &end_particle,
                                       const std::vector<Tensor<1,dim> > &old_velocities,
                                       const std::vector<Tensor<1,dim> > &,
                                       const double dt)
      {
        Assert(static_cast<unsigned int> (std::distance(begin_particle, end_particle)) == old_velocities.size(),
               ExcMessage("The particle integrator expects the velocity vector to be of equal size "
                          "to the number of particles to advect. For some unknown reason they are different, "
                          "most likely something went wrong in the calling function."));

        typename std::vector<Tensor<1,dim> >::const_iterator old_velocity = old_velocities.begin();

        for (typename std::multimap<types::LevelInd, Particle<dim> >::iterator it = begin_particle;
             it != end_particle; ++it, ++old_velocity)
          {
            const Point<dim> loc = it->second.get_location();
            it->second.set_location(loc + dt * (*old_velocity));
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
      ASPECT_REGISTER_PARTICLE_INTEGRATOR(Euler,
                                          "euler",
                                          "Explicit Euler scheme integrator, where y_{n+1} = y_n + dt * v(y_n). "
                                          "This requires only one integration substep per timestep.")
    }
  }
}
