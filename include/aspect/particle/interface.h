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

#ifndef _aspect_particle_interface_h
#define _aspect_particle_interface_h

#include <aspect/plugins.h>

namespace aspect
{
  namespace Particle
  {
    /**
     * A base class defining methods and member variables
     * shared along all the particle plugins. This includes the generator,
     * integrator, interpolator and the property classes.
     *
     * @ingroup Particle
     */
    class ParticleInterfaceBase : public Plugins::InterfaceBase
    {
      public:
        /**
         * Constructor.
         */
        ParticleInterfaceBase();

        /**
         * @brief Set which particle manager the plugin belongs to.
         *
         * @param particle_manager_index The index of the particle manager this plugin belongs to.
         */
        void set_particle_manager_index(const unsigned int particle_manager_index);

        /**
         * @brief Gets which particle manager the plugin belong to.
         */
        unsigned int get_particle_manager_index() const;

      private:
        /**
         * Stores the index to the particle manager, to which the plugin belongs.
         */
        unsigned int particle_manager_index;
    };

  }
}


#endif
