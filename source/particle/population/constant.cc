// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2015-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

#include <aspect/particle/population/constant.h>
#include <aspect/particle/manager.h>


namespace aspect
{
  namespace Particle
  {
    namespace Population
    {
      template <int dim>
      std::pair<unsigned int, unsigned int>
      Constant<dim>::particle_count_bounds(const typename DoFHandler<dim>::active_cell_iterator &) const
      {
        const Particle::Manager<dim> &manager = this->get_particle_manager(this->get_particle_manager_index());
        return {manager.get_min_particles_per_cell(), manager.get_max_particles_per_cell()};
      }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Particle
  {
    namespace Population
    {
      ASPECT_REGISTER_PARTICLE_POPULATION(Constant,
                                          "constant",
                                          "Prescribes a constant minimum and maximum number "
                                          "of particles allowed in each cell. "
                                          "The corresponding input parameters are read from the "
                                          "general ``Particles/'' section in the input file. "
                                          "Minimum and maximum particle counts are only enforced "
                                          "if particles are added or removed as part of the "
                                          "particle load balancing.")
    }
  }
}
