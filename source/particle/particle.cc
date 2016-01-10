/*
  Copyright (C) 2015 by the authors of the ASPECT code.

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

#include <aspect/particle/particle.h>

namespace aspect
{
  namespace Particle
  {
    template <int dim>
    Particle<dim>::Particle (const Point<dim> &new_loc,
                             const types::particle_index &new_id)
      :
      location (new_loc),
      id (new_id),
      properties ()
    {
    }

    template <int dim>
    Particle<dim>::Particle ()
      :
      location (),
      id (0),
      properties()
    {
    }


    template <int dim>
    Particle<dim>::Particle (const void *&data,
                             const unsigned int data_size)
    {
      // The data_size includes the space for position and id, so the number
      // of properties is the total size minus the space for position and id
      // divided by the size of one double (currently we only allow doubles as
      // tracer properties).
      const unsigned int property_size = data_size - dim * sizeof(double) - sizeof(types::particle_index);
      properties.resize(property_size / sizeof(double));

      const types::particle_index *id_data = static_cast<const types::particle_index *> (data);
      id = *id_data++;
      const double *pdata = reinterpret_cast<const double *> (id_data);

      for (unsigned int i = 0; i < dim; ++i)
        location(i) = *pdata++;

      for (unsigned int i = 0; i < properties.size(); ++i)
        properties [i] = *pdata++;

      data = static_cast<const void *> (pdata);
    }


    template <int dim>
    Particle<dim>::~Particle ()
    {
    }


    template <int dim>
    void
    Particle<dim>::set_location (const Point<dim> &new_loc)
    {
      location = new_loc;
    }

    template <int dim>
    void
    Particle<dim>::set_n_property_components (const unsigned int n_components)
    {
      properties.resize(n_components);
    }

    template <int dim>
    void
    Particle<dim>::write_data (void *&data) const
    {
      types::particle_index *id_data  = static_cast<types::particle_index *> (data);
      *id_data = id;
      ++id_data;
      double *pdata = reinterpret_cast<double *> (id_data);

      // Write location data
      for (unsigned int i = 0; i < dim; ++i,++pdata)
        *pdata =location(i);

      // Write property data
      for (unsigned int i = 0; i < properties.size(); ++i,++pdata)
        *pdata = properties[i];

      data = static_cast<void *> (pdata);
    }


    template <int dim>
    const Point<dim> &
    Particle<dim>::get_location () const
    {
      return location;
    }

    template <int dim>
    types::particle_index
    Particle<dim>::get_id () const
    {
      return id;
    }

    template <int dim>
    void
    Particle<dim>::set_properties (const std::vector<double> &new_properties)
    {
      properties = new_properties;
    }

    template <int dim>
    const std::vector<double> &
    Particle<dim>::get_properties () const
    {
      return properties;
    }

    template <int dim>
    std::vector<double> &
    Particle<dim>::get_properties ()
    {
      return properties;
    }
  }
}


// explicit instantiation of the functions we implement in this file
namespace aspect
{
  namespace Particle
  {
#define INSTANTIATE(dim) \
  template class Particle<dim>;

    ASPECT_INSTANTIATE(INSTANTIATE)
  }
}

