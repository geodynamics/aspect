/*
 Copyright (C) 2011, 2012, 2013 by the authors of the ASPECT code.

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
/*  $Id$  */

#include <aspect/particle/particle.h>

namespace aspect
{
  namespace Particle
  {
    // Base class of particles - represents a particle with position, velocity, and an ID number
    template <int dim>
    inline
    BaseParticle<dim>::BaseParticle (const Point<dim>& new_loc,
                                  const double& new_id)
                                  :
                                  location (new_loc),
                                  _id (new_id),
                                  is_local (true),
                                  check_vel (true)
                                  {
                                  }

    // Base class of particles - represents a particle with position, velocity, and an ID number
    template <int dim>
    inline
    BaseParticle<dim>::BaseParticle ()
    :
    location (),
    velocity (),
    _id (0),
    is_local (true),
    check_vel (true)
    {
    }


    // Base class of particles - represents a particle with position, velocity, and an ID number
    template <int dim>
    inline
    BaseParticle<dim>::~BaseParticle ()
    {
    }


    // Base class of particles - represents a particle with position, velocity, and an ID number
    template <int dim>
    void
    BaseParticle<dim>::set_location (const Point<dim> &new_loc)
    {
      location = new_loc;
    }


    // Base class of particles - represents a particle with position, velocity, and an ID number
    template <int dim>
    char*
    BaseParticle<dim>::write_data (ParticleDataFormat format,
                                   char* data) const
    {
      char* p = data;
      unsigned int i;
      // Then write our data in the appropriate format
      switch (format)
      {
      case MPI_DATA:
      case HDF5_DATA:
        // Write location data
        for (i = 0; i < dim; ++i)
          {
            double val = location (i);
            memcpy (p, &val, sizeof(double));
            p += sizeof(double);
          }
        // Write velocity data
        for (i = 0; i < dim; ++i)
          {
            double val = velocity (i);
            memcpy (p, &val, sizeof(double));
            p += sizeof(double);
          }
        double val = _id;
        memcpy (p, &val, sizeof(double));
        p += sizeof(double);
        break;
      }
      return p;
                                }

    // Base class of particles - represents a particle with position, velocity, and an ID number
    template <int dim>
    const char*
    BaseParticle<dim>::read_data (ParticleDataFormat format,
                               const char* data)
    {
      const char* p = data;
      unsigned int i;
      switch (format)
      {
      case MPI_DATA:
      case HDF5_DATA:
        // Read location data
        for (i = 0; i < dim; ++i)
          {
            double val;
            memcpy (&val, p, sizeof(double));
            location (i) = val;
            p += sizeof(double);
          }
        // Write velocity data
        for (i = 0; i < dim; ++i)
          {
            double val;
            memcpy (&val, p, sizeof(double));
            velocity (i) = val;
            p += sizeof(double);
          }
        double val;
        memcpy (&val, p, sizeof(double));
        _id = val;
        p += sizeof(double);
        break;
      }
      return p;
                               }

    // Base class of particles - represents a particle with position, velocity, and an ID number
    template <int dim>
    unsigned int
    BaseParticle<dim>::data_len (ParticleDataFormat format)
    {
      switch (format)
      {
      case MPI_DATA:
      case HDF5_DATA:
        return (dim + dim + 1) * sizeof(double);
      }
      return 0;
    }


    // explicit instantiation
    template class BaseParticle<2>;
    template class BaseParticle<3>;
  }
}
