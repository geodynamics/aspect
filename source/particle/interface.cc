// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2024 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

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
