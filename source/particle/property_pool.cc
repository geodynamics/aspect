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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/

#include <aspect/particle/property_pool.h>
#include <aspect/particle/particle.h>

#include <aspect/global.h>
#include <aspect/utilities.h>
#include <aspect/compat.h>

namespace aspect
{
  namespace Particle
  {
#if DEAL_II_VERSION_GTE(8,5,0)
    const PropertyPool::Handle PropertyPool::invalid_handle = numbers::invalid_unsigned_int;
#else
    const PropertyPool::Handle PropertyPool::invalid_handle = NULL;
#endif


    PropertyPool::PropertyPool (const unsigned int n_properties_per_slot)
      :
      n_properties (n_properties_per_slot)
#if DEAL_II_VERSION_GTE(8,5,0)
      ,
      memory_pool(),
      translation_map(),
      free_slots()
#endif
    {}



    PropertyPool::Handle
    PropertyPool::allocate_properties_array ()
    {
#if DEAL_II_VERSION_GTE(8,5,0)
      if (n_properties == 0)
        return numbers::invalid_unsigned_int;

      if (free_slots.is_empty())
        reserve(free_slots.size() * 2 + 1);

      const Handle handle = free_slots.pop_back();
      translation_map.insert(translation_entry(handle, handle * n_properties));

      return handle;
#else
      return new double[n_properties];
#endif
    }



    void
    PropertyPool::deallocate_properties_array (Handle handle)
    {
#if DEAL_II_VERSION_GTE(8,5,0)
      free_slots.add_index(handle);
#else
      delete[] handle;
#endif
    }



    ArrayView<double>
    PropertyPool::get_properties (const Handle handle)
    {
#if DEAL_II_VERSION_GTE(8,5,0)
      if (n_properties == 0)
        return ArrayView<double>(NULL,0);

      Assert(handle < free_slots.size(),
             ExcMessage("The particle property pool was asked to access a "
                        "region of memory that is not allocated for particle properties."));
      Assert(!free_slots.is_element(handle),
             ExcMessage("The particle property pool was asked to access a "
                        "region of memory that is allocated for particle properties, "
                        " but was reported as not in use."));

      // Here we make use of the fact that Handle is nothing else than an index
      // of an IndexSet (i.e. a unsigned number between 0 and free_slots.size())
      return ArrayView<double>(&memory_pool[translation_map.left.at(handle)],n_properties);
#else
      return ArrayView<double>(handle, n_properties);
#endif
    }


    void
    PropertyPool::reserve(const std::size_t size)
    {
#if DEAL_II_VERSION_GTE(8,5,0)
      if (size > free_slots.size())
        {
          memory_pool.resize(size * n_properties);

          IndexSet new_free_slots (size);

          for (IndexSet::IntervalIterator r = free_slots.begin_intervals(); r != free_slots.end_intervals(); ++r)
            new_free_slots.add_range (*(r->begin()),r->last()+1);  // a previously free range is still free

          new_free_slots.add_range (free_slots.size(), size);  // all new slots are free

          free_slots = new_free_slots;
        }
#else
      (void)size;
#endif
    }
  }
}
