/*
  Copyright (C) 2015 - 2024 by the authors of the ASPECT code.

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

#include <aspect/particle/integrator/euler.h>
#include <aspect/geometry_model/interface.h>

namespace aspect
{
  namespace Particle
  {
    namespace Integrator
    {
      template <int dim>
      void
      Euler<dim>::local_integrate_step(const typename ParticleHandler<dim>::particle_iterator &begin_particle,
                                       const typename ParticleHandler<dim>::particle_iterator &end_particle,
                                       const std::vector<Tensor<1,dim>> &old_velocities,
                                       const std::vector<Tensor<1,dim>> &,
                                       const double dt)
      {
        Assert(static_cast<unsigned int> (std::distance(begin_particle, end_particle)) == old_velocities.size(),
               ExcMessage("The particle integrator expects the velocity vector to be of equal size "
                          "to the number of particles to advect. For some unknown reason they are different, "
                          "most likely something went wrong in the calling function."));

        const auto cell = begin_particle->get_surrounding_cell();
        bool at_periodic_boundary = false;
        if (this->get_triangulation().get_periodic_face_map().empty() == false)
          for (const auto face_index: cell->face_indices())
            if (cell->at_boundary(face_index))
              if (cell->has_periodic_neighbor(face_index))
                {
                  at_periodic_boundary = true;
                  break;
                }

        typename std::vector<Tensor<1,dim>>::const_iterator old_velocity = old_velocities.begin();

        for (typename ParticleHandler<dim>::particle_iterator it = begin_particle;
             it != end_particle; ++it, ++old_velocity)
          {
#if DEAL_II_VERSION_GTE(9, 6, 0)
            // Get a reference to the particle location, so that we can update it in-place
            Point<dim> &location = it->get_location();
#else
            Point<dim> location = it->get_location();
#endif
            location += dt * (*old_velocity);

            if (at_periodic_boundary)
              this->get_geometry_model().adjust_positions_for_periodicity(location);

#if !DEAL_II_VERSION_GTE(9, 6, 0)
            it->set_location(location);
#endif
          }
      }



      template <int dim>
      std::array<bool, 3>
      Euler<dim>::required_solution_vectors() const
      {
        return {{false, true, false}};
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
                                          "Explicit Euler scheme integrator, where "
                                          "$y_{n+1} = y_n + \\Delta t \\, v(y_n)$. "
                                          "This requires only one integration substep per timestep.")
    }
  }
}
