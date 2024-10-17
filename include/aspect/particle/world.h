/*
 Copyright (C) 2012 - 2024 by the authors of the ASPECT code.

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
