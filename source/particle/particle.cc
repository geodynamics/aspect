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

#include <list>

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
      is_local (true),
      val ()
    {
    }

    template <int dim>
    inline
    Particle<dim>::Particle ()
      :
      location (),
      id (0),
      is_local (true),
      val()
    {
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
    Particle<dim>::write_data (std::vector<double> &data) const
    {
      // Write location data
      for (unsigned int i = 0; i < dim; ++i)
        {
          data.push_back(location(i));
        }

      data.push_back(id);

      for (unsigned int i = 0; i < val.size(); i++)
        {
          data.push_back(val [i]);
        }
    }

    template <int dim>
    unsigned int
    Particle<dim>::read_data(const std::vector<double> &data, const unsigned int pos)
    {
      unsigned int p = pos;
      // Read location data
      for (unsigned int i = 0; i < dim; ++i)
        {
          location (i) = data[p++];
        }

      id = data[p++];

      for (unsigned int i = 0; i < val.size(); ++i)
        {
          val [i] = data[p++];
        }

      return p;
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
    double
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

    template <int dim>
    bool
    Particle<dim>::local () const
    {
      return is_local;
    }

    template <int dim>
    void
    Particle<dim>::set_local (bool new_local)
    {
      is_local = new_local;
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

