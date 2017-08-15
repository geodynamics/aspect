/*
  Copyright (C) 2015 - 2016 by the authors of the ASPECT code.

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

#include <aspect/particle/particle.h>

namespace aspect
{
  namespace Particle
  {
    template <int dim, int spacedim>
    Particle<dim,spacedim>::Particle ()
      :
      location (),
      reference_location(),
      id (0),
      property_pool(NULL),
      properties(PropertyPool::invalid_handle)
    {
    }


    template <int dim, int spacedim>
    Particle<dim,spacedim>::Particle (const Point<spacedim> &location,
                                      const Point<dim> &reference_location,
                                      const types::particle_index id)
      :
      location (location),
      reference_location (reference_location),
      id (id),
      property_pool(NULL),
      properties (PropertyPool::invalid_handle)
    {
    }

    template <int dim, int spacedim>
    Particle<dim,spacedim>::Particle (const Particle<dim,spacedim> &particle)
      :
      location (particle.get_location()),
      reference_location (particle.get_reference_location()),
      id (particle.get_id()),
      property_pool(particle.property_pool),
      properties ((property_pool != NULL) ? property_pool->allocate_properties_array() : PropertyPool::invalid_handle)
    {
      if (property_pool != NULL)
        {
          const ArrayView<double> my_properties = property_pool->get_properties(properties);

          if (my_properties.size() != 0)
            {
              const ArrayView<const double> their_properties = particle.get_properties();

              std::copy(&their_properties[0],&their_properties[0]+their_properties.size(),&my_properties[0]);
            }
        }
    }


    template <int dim, int spacedim>
    Particle<dim,spacedim>::Particle (const void *&data,
                                      PropertyPool &new_property_pool)
    {
      const types::particle_index *id_data = static_cast<const types::particle_index *> (data);
      id = *id_data++;
      const double *pdata = reinterpret_cast<const double *> (id_data);

      for (unsigned int i = 0; i < dim; ++i)
        location(i) = *pdata++;

      for (unsigned int i = 0; i < dim; ++i)
        reference_location(i) = *pdata++;

      property_pool = &new_property_pool;
      properties = property_pool->allocate_properties_array();

      // See if there are properties to load
      const ArrayView<double> particle_properties = property_pool->get_properties(properties);
      for (unsigned int i = 0; i < particle_properties.size(); ++i)
        particle_properties[i] = *pdata++;

      data = static_cast<const void *> (pdata);
    }

#ifdef DEAL_II_WITH_CXX11

    template <int dim, int spacedim>
    Particle<dim,spacedim>::Particle (Particle<dim,spacedim> &&particle)
      :
      location (particle.location),
      reference_location(particle.reference_location),
      id (particle.id),
      property_pool(particle.property_pool),
      properties(particle.properties)
    {
      particle.properties = PropertyPool::invalid_handle;
    }

    template <int dim, int spacedim>
    Particle<dim,spacedim> &
    Particle<dim,spacedim>::operator=(const Particle<dim,spacedim> &particle)
    {
      if (this != &particle)
        {
          location = particle.location;
          reference_location = particle.reference_location;
          id = particle.id;
          property_pool = particle.property_pool;

          if (property_pool != NULL)
            {
              properties = property_pool->allocate_properties_array();
              const ArrayView<const double> their_properties = particle.get_properties();

              if (their_properties.size() != 0)
                {
                  const ArrayView<double> my_properties = property_pool->get_properties(properties);
                  std::copy(&their_properties[0],&their_properties[0]+their_properties.size(),&my_properties[0]);
                }
            }
          else
            properties = PropertyPool::invalid_handle;
        }
      return *this;
    }

    template <int dim, int spacedim>
    Particle<dim,spacedim> &
    Particle<dim,spacedim>::operator=(Particle<dim,spacedim> &&particle)
    {
      if (this != &particle)
        {
          location = particle.location;
          reference_location = particle.reference_location;
          id = particle.id;
          property_pool = particle.property_pool;
          properties = particle.properties;
          particle.properties = PropertyPool::invalid_handle;
        }
      return *this;
    }
#endif

    template <int dim, int spacedim>
    Particle<dim,spacedim>::~Particle ()
    {
      if (properties != PropertyPool::invalid_handle)
        property_pool->deallocate_properties_array(properties);
    }

    template <int dim, int spacedim>
    void
    Particle<dim,spacedim>::write_data (void *&data) const
    {
      types::particle_index *id_data  = static_cast<types::particle_index *> (data);
      *id_data = id;
      ++id_data;
      double *pdata = reinterpret_cast<double *> (id_data);

      // Write location data
      for (unsigned int i = 0; i < dim; ++i,++pdata)
        *pdata = location(i);

      // Write reference location data
      for (unsigned int i = 0; i < dim; ++i,++pdata)
        *pdata = reference_location(i);

      // Write property data
      const ArrayView<double> particle_properties = property_pool->get_properties(properties);
      for (unsigned int i = 0; i < particle_properties.size(); ++i,++pdata)
        *pdata = particle_properties[i];

      data = static_cast<void *> (pdata);
    }

    template <int dim, int spacedim>
    void
    Particle<dim,spacedim>::set_location (const Point<spacedim> &new_loc)
    {
      location = new_loc;
    }

    template <int dim, int spacedim>
    const Point<spacedim> &
    Particle<dim,spacedim>::get_location () const
    {
      return location;
    }

    template <int dim, int spacedim>
    void
    Particle<dim,spacedim>::set_reference_location (const Point<dim> &new_loc)
    {
      reference_location = new_loc;
    }

    template <int dim, int spacedim>
    const Point<dim> &
    Particle<dim,spacedim>::get_reference_location () const
    {
      return reference_location;
    }

    template <int dim, int spacedim>
    types::particle_index
    Particle<dim,spacedim>::get_id () const
    {
      return id;
    }

    template <int dim, int spacedim>
    void
    Particle<dim,spacedim>::set_property_pool (PropertyPool &new_property_pool)
    {
      property_pool = &new_property_pool;
    }

    template <int dim, int spacedim>
    void
    Particle<dim,spacedim>::set_properties (const std::vector<double> &new_properties)
    {
      if (properties == PropertyPool::invalid_handle)
        properties = property_pool->allocate_properties_array();

      const ArrayView<double> old_properties = property_pool->get_properties(properties);

      std::copy(new_properties.begin(),new_properties.end(),&old_properties[0]);
    }

    template <int dim, int spacedim>
    const ArrayView<const double>
    Particle<dim,spacedim>::get_properties () const
    {
      Assert(property_pool != NULL,
             ExcInternalError());

      return property_pool->get_properties(properties);
    }


    template <int dim, int spacedim>
    const ArrayView<double>
    Particle<dim,spacedim>::get_properties ()
    {
      Assert(property_pool != NULL,
             ExcInternalError());

      return property_pool->get_properties(properties);
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

