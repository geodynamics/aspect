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
                                  id (new_id),
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
    id (0),
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
        double val = id;
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
        id = val;
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


    // Base class of particles - represents a particle with position, velocity, and an ID number
    template <int dim>
    Point<dim>
    BaseParticle<dim>::get_location () const
    {
      return location;
    }

    // Base class of particles - represents a particle with position, velocity, and an ID number
    template <int dim>
    void
    BaseParticle<dim>::set_velocity (Point<dim> new_vel)
    {
      velocity = new_vel;
    }

    // Base class of particles - represents a particle with position, velocity, and an ID number
    template <int dim>
    Point<dim>
    BaseParticle<dim>::get_velocity () const
    {
      return velocity;
    }

    // Base class of particles - represents a particle with position, velocity, and an ID number
    template <int dim>
    double
    BaseParticle<dim>::get_id () const
    {
      return id;
    }

    // Base class of particles - represents a particle with position, velocity, and an ID number
    template <int dim>
    bool
    BaseParticle<dim>::local () const
    {
      return is_local;
    }

    // Base class of particles - represents a particle with position, velocity, and an ID number
    template <int dim>
    void
    BaseParticle<dim>::set_local (bool new_local)
    {
      is_local = new_local;
    }

    // Base class of particles - represents a particle with position, velocity, and an ID number
    template <int dim>
    bool
    BaseParticle<dim>::vel_check () const
    {
      return check_vel;
    }

    // Base class of particles - represents a particle with position, velocity, and an ID number
    template <int dim>
    void
    BaseParticle<dim>::set_vel_check (bool new_vel_check)
    {
      check_vel = new_vel_check;
    }

    // Base class of particles - represents a particle with position, velocity, and an ID number
    template <int dim>
    void
    BaseParticle<dim>::add_mpi_types (std::vector<MPIDataInfo>& data_info)
    {
      // Add the position, velocity, ID
      data_info.push_back (
          MPIDataInfo ("pos", dim, MPI_DOUBLE, sizeof(double)));
      data_info.push_back (
          MPIDataInfo ("velocity", dim, MPI_DOUBLE, sizeof(double)));
      data_info.push_back (MPIDataInfo ("id", 1, MPI_DOUBLE, sizeof(double)));
    }



    // explicit instantiation
    template class BaseParticle<2>;
    template class BaseParticle<3>;
  }
}
