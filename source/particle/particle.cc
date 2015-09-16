/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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
    inline
    Particle<dim>::Particle (const Point<dim> &new_loc,
                             const double &new_id)
      :
      location (new_loc),
      id (new_id),
      val ()
    {
    }

    template <int dim>
    inline
    Particle<dim>::Particle ()
      :
      location (),
      id (0),
      val()
    {
    }


    template <int dim>
    inline
    Particle<dim>::Particle (const void *&data,
                             const unsigned int data_len)
      :
      val(data_len-dim-1)
    {
      const unsigned int *id_data = static_cast<const unsigned int *> (data);
      id = *id_data++;
      const double *pdata = reinterpret_cast<const double *> (id_data);

      for (unsigned int i = 0; i < dim; ++i)
        location(i) = *pdata++;

      for (unsigned int i = 0; i < val.size(); ++i)
        val [i] = *pdata++;

      data = static_cast<const void *> (pdata);
    }


    template <int dim>
    inline
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
    Particle<dim>::set_data_len (const unsigned int data_len)
    {
      val.resize(data_len - dim - 1);
    }

    template <int dim>
    void
    Particle<dim>::write_data (void *&data) const
    {
      unsigned int *id_data  = static_cast<unsigned int *> (data);
      *id_data = id;
      ++id_data;
      double *pdata = reinterpret_cast<double *> (id_data);

      // Write location data
      for (unsigned int i = 0; i < dim; ++i,++pdata)
        *pdata =location(i);

      // Write property data
      for (unsigned int i = 0; i < val.size(); ++i,++pdata)
        *pdata = val[i];

      data = static_cast<void *> (pdata);
    }

    template <int dim>
    unsigned int
    Particle<dim>::data_len () const
    {
      return dim + 1 + val.size();
    }


    template <int dim>
    const Point<dim> &
    Particle<dim>::get_location () const
    {
      return location;
    }

    template <int dim>
    unsigned int
    Particle<dim>::get_id () const
    {
      return id;
    }

    template <int dim>
    void
    Particle<dim>::set_properties (const std::vector<double> &new_properties)
    {
      val = new_properties;
    }

    template <int dim>
    const std::vector<double> &
    Particle<dim>::get_properties () const
    {
      return val;
    }

    template <int dim>
    std::vector<double> &
    Particle<dim>::get_properties ()
    {
      return val;
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

