// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2012-2024 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

#ifndef _aspect_particle_world_h
#define _aspect_particle_world_h

#include <aspect/particle/manager.h>

namespace aspect
{
  namespace Particle
  {
    /**
     * The World class is a deprecated name for the Manager class. It is
     * provided for backwards compatibility only and will be removed in
     * the future. Please use the Manager class instead.
     *
     * @deprecated: This class is deprecated and will be removed in the future.
     * Please use the Particle::Manager class instead.
     */
    template <int dim>
    using World = Manager<dim>;
  }
}

#endif
