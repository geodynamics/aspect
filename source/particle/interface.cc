/*
  Copyright (C) 2024 by the authors of the ASPECT code.

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

#include <aspect/particle/interface.h>
#include <deal.II/base/numbers.h>

namespace aspect
{
  namespace Particle
  {
    ParticleInterfaceBase::ParticleInterfaceBase()
      : particle_manager_index(dealii::numbers::invalid_unsigned_int)
    {}



    void
    ParticleInterfaceBase::set_particle_manager_index(const unsigned int particle_manager_index)
    {
      this->particle_manager_index = particle_manager_index;
    }



    unsigned int
    ParticleInterfaceBase::get_particle_manager_index() const
    {
      return particle_manager_index;
    }
  }
}
