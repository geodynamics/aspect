/*
 Copyright (C) 2015 - 2022 by the authors of the ASPECT code.

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

#ifndef _aspect_particle_interface_h
#define _aspect_particle_interface_h

#include "aspect/plugins.h"

namespace aspect
{
  namespace Particle
  {
    /**
     * An abstract base class defining methods and member variables
     * shared along all the particle plugins. This includes the generator,
     * integrator, interpolator and the property classes.
     *
     * @ingroup Particle
     */
    class ParticleInterfaceBase : public Plugins::InterfaceBase
    {
      public:
        /**
         * @brief Set which particle world the plugin belong to.
         *
         * @param particle_world_index The index of the particle world this plugin belongs to.
         */
        void set_particle_world_index(unsigned int particle_world_index);

        /**
         * @brief Gets which particle world the plugin belong to.
         */
        unsigned int get_particle_world_index() const;

      private:
        /**
         * Stores the index to the particle world in the
         * particle worlds vector, to which the plugin belongs.
         */
        unsigned int particle_world_index;
    };

  }
}


#endif
