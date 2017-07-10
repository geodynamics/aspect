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

#include <aspect/particle/property_pool.h>
#include <aspect/particle/particle.h>

namespace aspect
{
  namespace Particle
  {
    const PropertyPool::Handle PropertyPool::invalid_handle = NULL;


    PropertyPool::PropertyPool (const unsigned int n_properties_per_slot)
      :
      n_properties (n_properties_per_slot)
    {}



    PropertyPool::Handle
    PropertyPool::allocate_properties_array ()
    {
      return new double[n_properties];
    }



    void
    PropertyPool::deallocate_properties_array (Handle handle)
    {
      delete[] handle;
    }



    ArrayView<double>
    PropertyPool::get_properties (const Handle handle)
    {
      return ArrayView<double>(handle, n_properties);
    }


    void
    PropertyPool::reserve(const std::size_t size)
    {
      (void)size;
    }
  }
}
