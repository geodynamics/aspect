/*
 Copyright (C) 2011, 2012 by the authors of the ASPECT code.

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
#include <aspect/simulator.h>

namespace aspect
{
  namespace Particle
  {
    // Typedef of cell level/index pair
    typedef std::pair<int, int> LevelInd;

    class MPIDataInfo
    {
      public:
        std::string     _name;
        unsigned int    _num_elems;
        MPI_Datatype    _data_type;
        unsigned int    _elem_size_bytes;

        MPIDataInfo(std::string name, unsigned int num_elems, MPI_Datatype data_type, unsigned int elem_size_bytes) : _name(name), _num_elems(num_elems), _data_type(data_type), _elem_size_bytes(elem_size_bytes) {};
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
        Point<dim>      _loc;

        // Current particle velocity
        Point<dim>      _vel;

        // Globally unique ID of particle
        double          _id;

        // Whether this particle is in the local subdomain or not
        bool            _local;

        // Whether to check the velocity of this particle
        // This is used for integration schemes which require multiple
        // integration steps for some particles, but not for others
        bool            _check_vel;

      public:
        BaseParticle(void) : _loc(), _vel(), _id(0), _local(true), _check_vel(true) {};

        BaseParticle(const Point<dim> &new_loc, const double &new_id) : _loc(new_loc), _id(new_id), _local(true), _check_vel(true) {};

        virtual ~BaseParticle(void) {};
        
        static unsigned int data_len(ParticleDataFormat format)
        {
          switch (format)
            {
              case MPI_DATA:
              case HDF5_DATA:
                return (dim+dim+1)*sizeof(double);
            }
          return 0;
        };
        virtual char *read_data(ParticleDataFormat format, char *data)
        {
          char            *p = data;
          unsigned int    i;

          switch (format)
            {
              case MPI_DATA:
              case HDF5_DATA:
                // Read location data
                for (i=0; i<dim; ++i)
                  {
                    _loc(i) = ((double *)p)[0];
                    p += sizeof(double);
                  }
                // Write velocity data
                for (i=0; i<dim; ++i)
                  {
                    _vel(i) = ((double *)p)[0];
                    p += sizeof(double);
                  }
                _id = ((double *)p)[0];
                p += sizeof(double);
                break;
            }

          return p;
        };
        virtual char *write_data(ParticleDataFormat format, char *data) const
        {
          char          *p = data;
          unsigned int  i;

          // Then write our data in the appropriate format
          switch (format)
            {
              case MPI_DATA:
              case HDF5_DATA:
                // Write location data
                for (i=0; i<dim; ++i)
                  {
                    ((double *)p)[0] = _loc(i);
                    p += sizeof(double);
                  }
                // Write velocity data
                for (i=0; i<dim; ++i)
                  {
                    ((double *)p)[0] = _vel(i);
                    p += sizeof(double);
                  }
                ((double *)p)[0] = _id;
                p += sizeof(double);
                break;
            }

          return p;
        };

        void set_location(Point<dim> new_loc)
        {
          _loc = new_loc;
        };
        Point<dim> location(void) const
        {
          return _loc;
        };

        void set_velocity(Point<dim> new_vel)
        {
          _vel = new_vel;
        };
        Point<dim> velocity(void) const
        {
          return _vel;
        };

        double id_num(void) const
        {
          return _id;
        };

        bool local(void) const
        {
          return _local;
        };
        void set_local(bool new_local)
        {
          _local = new_local;
        };

        bool vel_check(void) const
        {
          return _check_vel;
        };
        void set_vel_check(bool new_vel_check)
        {
          _check_vel = new_vel_check;
        };

        static void add_mpi_types(std::vector<MPIDataInfo> &data_info)
        {
          // Add the position, velocity, ID
          data_info.push_back(MPIDataInfo("pos", dim, MPI_DOUBLE, sizeof(double)));
          data_info.push_back(MPIDataInfo("velocity", dim, MPI_DOUBLE, sizeof(double)));
          data_info.push_back(MPIDataInfo("id", 1, MPI_DOUBLE, sizeof(double)));
        };
    };

    // A particle with associated values, such as scalars, vectors or tensors
    template <int dim, int data_dim>
    class DataParticle : public BaseParticle<dim>
    {
      private:
        double      _val[data_dim];

      public:
        DataParticle(void)
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
        virtual char *read_data(ParticleDataFormat format, char *data)
        {
          char          *p = data;
          unsigned int  i;

          // Read the parent data first
          p = BaseParticle<dim>::read_data(format, data);

          // Then read our data in the appropriate format
          switch (format)
            {
              case MPI_DATA:
                for (i=0; i<data_dim; ++i)
                  {
                    _val[i] = ((double *)p)[0];
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
                    ((double *)p)[0] = _val[i];
                    p += sizeof(double);
                  }
                break;
              case HDF5_DATA:
                break;
            }

          return p;
        };

        // Returns a vector from the first dim components of _val
        Point<dim> get_vector(void) const
        {
          AssertThrow(data_dim>=dim, std::out_of_range("get_vector"));
          Point<dim>  p;
          for (unsigned int i=0; i<dim; ++i) p(i) = _val[i];
        };
        // Sets the first dim components of _val to the specified vector value
        void set_vector(Point<dim> new_vec)
        {
          AssertThrow(data_dim>=dim, std::out_of_range("set_vector"));
          for (unsigned int i=0; i<dim; ++i) _val[i] = new_vec(i);
        };

        double &operator[](const unsigned int &ind)
        {
          AssertThrow(data_dim>ind, std::out_of_range("DataParticle[]"));
          return _val[ind];
        };
        const double operator[](const unsigned int &ind) const
        {
          AssertThrow(data_dim>ind, std::out_of_range("DataParticle[]"));
          return _val[ind];
        };

        static void add_mpi_types(std::vector<MPIDataInfo> &data_info)
        {
          // Set up the parent types first
          BaseParticle<dim>::add_mpi_types(data_info);

          // Then add our own
          data_info.push_back(MPIDataInfo("data", data_dim, MPI_DOUBLE, sizeof(double)));
        };
    };
  }
}

#endif
