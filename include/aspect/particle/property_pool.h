/*
 Copyright (C) 2016 by the authors of the ASPECT code.

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

#ifndef _aspect_particle_property_pool_h
#define _aspect_particle_property_pool_h

#include <aspect/global.h>

#include <deal.II/base/array_view.h>

namespace aspect
{
  namespace Particle
  {
    using namespace dealii;

    /**
     * This class manages the memory space in which particles store their
     * properties. Because this is dynamic memory and every particle needs the
     * same amount it is more efficient to let this be handled by a central
     * manager that does not need to allocate/deallocate memory every time a
     * particle is constructed/destroyed.
     */
    class PropertyPool
    {
      public:
        /**
         * Typedef for the handle that is returned to the particles, and that
         * uniquely identifies the slot of memory that is reserved for this
         * particle.
         */
        typedef double *Handle;

        /**
         * Define a default (invalid) value for handles.
         */
        static const Handle invalid_handle;

        /**
         * Constructor. Stores the number of properties per reserved slot.
         */
        PropertyPool (const unsigned int n_properties_per_slot);

        /**
         * Returns a new handle that allows accessing the reserved block
         * of memory.
         */
        Handle allocate_properties_array ();

        /**
         * Mark the properties corresponding to the handle @p handle as
         * deleted. After calling this function calling get_properties() with
         * the freed handle causes undefined behavior.
         */
        void deallocate_properties_array (const Handle handle);

        /**
         * Return an ArrayView to the properties that correspond to the given
         * handle @p handle.
         */
        ArrayView<double> get_properties (const Handle handle);

        /**
         * Reserves the dynamic memory needed for storing the properties of
         * @p size particles.
         */
        void reserve(const std::size_t size);

      private:
        /**
         * The number of properties that are reserved per particle.
         */
        const unsigned int n_properties;
    };

  }
}

#endif
