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

#ifndef __aspect__particle_particle_h
#define __aspect__particle_particle_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace Particle
  {
    // Typedef of cell level/index pair
    typedef std::pair<int, int> LevelInd;

    class MPIDataInfo
    {
      public:
        std::string     name;
        unsigned int    n_elements;
        MPI_Datatype    data_type;
        unsigned int    size_in_bytes;

        MPIDataInfo(std::string name,
                    unsigned int num_elems,
                    MPI_Datatype data_type,
                    unsigned int elem_size_bytes)
        :
          name(name),
          n_elements(num_elems),
          data_type(data_type),
          size_in_bytes(elem_size_bytes) {};
    };

    enum ParticleDataFormat
    {
      MPI_DATA,
      HDF5_DATA
    };

    // Base class of particles - represents a particle with position, velocity, and an ID number
    template <int dim>
    class BaseParticle
    {
      private:
        // Current particle location
        Point<dim>      location;

        // Current particle velocity
        Point<dim>      velocity;

        // Globally unique ID of particle
//TODO: make this an unsigned int. but this needs adjustment in data_len, read_data, write_data and all of the places that call write_data and then parse the output...
        double          id;

        // Whether this particle is in the local subdomain or not
        bool            is_local;

        // Whether to check the velocity of this particle
        // This is used for integration schemes which require multiple
        // integration steps for some particles, but not for others
        bool            check_vel;

      public:
        BaseParticle ();

        BaseParticle (const Point<dim>& new_loc,
                      const double& new_id);

        virtual
        ~BaseParticle ();

        static unsigned int
        data_len (ParticleDataFormat format);

        virtual const char*
        read_data (ParticleDataFormat format,
                   const char* data);


        virtual char*
        write_data (ParticleDataFormat format,
                    char* data) const;

        void
        set_location (const Point<dim> &new_loc);

        Point<dim>
        get_location () const;

        void
        set_velocity (Point<dim> new_vel);
        Point<dim>
        get_velocity () const;

        double
        get_id () const;

        bool
        local () const;

        void
        set_local (bool new_local);

        bool
        vel_check () const;

        void
        set_vel_check (bool new_vel_check);

        static void
        add_mpi_types (std::vector<MPIDataInfo>& data_info);
    };

    // A particle with associated values, such as scalars, vectors or tensors
    template <int dim, int data_dim>
    class DataParticle : public BaseParticle<dim>
    {
      private:
        double      _val[data_dim];

      private:
        DataParticle()
        {
          for (unsigned int i=0; i<data_dim; ++i) _val[i] = 0;
        };

        DataParticle(const Point<dim> &new_loc, const double &new_id) : BaseParticle<dim>(new_loc, new_id)
        {
          for (unsigned int i=0; i<data_dim; ++i) _val[i] = 0;
        };

        static unsigned int data_len(ParticleDataFormat format)
        {
          unsigned int        base_size  = BaseParticle<dim>::data_len(format);

          switch (format)
            {
              case MPI_DATA:
              case HDF5_DATA:
                return base_size + data_dim * sizeof(double);
            }
          return 0;
        };

        virtual const char *read_data(ParticleDataFormat format, const char *data)
        {
          const char          *p = data;
          unsigned int  i;

          // Read the parent data first
          p = BaseParticle<dim>::read_data(format, data);

          // Then read our data in the appropriate format
          switch (format)
            {
              case MPI_DATA:
                for (i=0; i<data_dim; ++i)
                  {
                    double val;
                    memcpy (&val, p, sizeof(double));
                    _val[i] = val;
                    p += sizeof(double);
                  }
                break;
              case HDF5_DATA:
                break;
            }

          return p;
        };


        virtual char *write_data(ParticleDataFormat format, char *data) const
        {
          char          *p = data;
          unsigned int  i;

          // Write the parent data first
          p = BaseParticle<dim>::write_data(format, data);

          // Then write our data in the appropriate format
          switch (format)
            {
              case MPI_DATA:
                for (i=0; i<data_dim; ++i)
                  {
                    memcpy (p, _val[i], sizeof(double));
                    p += sizeof(double);
                  }
                break;
              case HDF5_DATA:
                break;
            }

          return p;
        };

        // Returns a vector from the first dim components of _val
        Point<dim>
        get_vector () const;
        
        // Sets the first dim components of _val to the specified vector value
        void set_vector(Point<dim> new_vec)
        {
          AssertThrow(data_dim>=dim, std::out_of_range("set_vector"));
          for (unsigned int i=0; i<dim; ++i)
            _val[i] = new_vec(i);
        }

        double &operator[](const unsigned int &ind)
        {
          AssertThrow(data_dim>ind, std::out_of_range("DataParticle[]"));
          return _val[ind];
        }

        double operator[](const unsigned int &ind) const
        {
          AssertThrow(data_dim>ind, std::out_of_range("DataParticle[]"));
          return _val[ind];
        }

        static void add_mpi_types(std::vector<MPIDataInfo> &data_info)
        {
          // Set up the parent types first
          BaseParticle<dim>::add_mpi_types(data_info);

          // Then add our own
          data_info.push_back(MPIDataInfo("data", data_dim, MPI_DOUBLE, sizeof(double)));
        };
    };

    // A particle with associated values, such as scalars, vectors or tensors
    template <int dim, int data_dim>
    inline Point<dim>
    DataParticle<dim,data_dim>::get_vector () const
    {
      AssertThrow(data_dim >= dim, std::out_of_range ("get_vector"));
      Point < dim > p;
      for (unsigned int i = 0; i < dim; ++i)
        p (i) = _val[i];
    }

   }
}

#endif
