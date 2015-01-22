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

#ifndef __aspect__particle_particle_h
#define __aspect__particle_particle_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace Particle
  {
    /**
     * Typedef of cell level/index pair
     */
    typedef std::pair<int, int> LevelInd;

    class MPIDataInfo
    {
      public:
        std::string     name;
        unsigned int    n_elements;

        MPIDataInfo(std::string name,
                    unsigned int num_elems)
          :
          name(name),
          n_elements(num_elems) {};
    };

    /**
     * Base class of particles - represents a particle with position,
     * velocity, and an ID number. This class can be extended to include data
     * related to a particle. An example of this is shown in the DataParticle
     * class.
     */
    template <int dim>
    class BaseParticle
    {
      private:
        /**
         * Current particle location
         */
        Point<dim>      location;

        /**
         * Current particle velocity
         */
        Point<dim>      velocity;

        /**
         * Globally unique ID of particle
         */
        double          id;

        /**
         * Whether this particle is in the local subdomain or not
         */
        bool            is_local;

        /**
         * Whether to check the velocity of this particle This is used for
         * integration schemes which require multiple integration steps for
         * some particles, but not for others
         */
        bool            check_vel;

      public:
        /**
         * Empty constructor for BaseParticle, creates a particle at the
         * origin with zero velocity.
         */
        BaseParticle ();

        /**
         * Constructor for BaseParticle, creates a particle with the specified
         * ID at the specified location with zero velocity. Note that Aspect
         * does not check for duplicate particle IDs so the user must be sure
         * the IDs are unique over all processes.
         *
         * @param[in] new_loc Initial location of particle.
         * @param[in] new_id Globally unique ID number of particle.
         */
        BaseParticle (const Point<dim> &new_loc,
                      const double &new_id);

        /**
         * Destructor for BaseParticle
         */
        virtual
        ~BaseParticle ();

        /**
         * Get the number of doubles required to represent this particle for
         * communication.
         *
         * @return Number of doubles required to represent this particle
         */
        static unsigned int
        data_len ();

        /**
         * Read the particle data from the specified vector of doubles.
         *
         * @param [in] data The vector of double data to read from.
         * @param [in] pos The position in the data vector to start reading
         * from.
         * @return The position in the vector of the next unread double.
         */
        virtual unsigned int read_data(const std::vector<double> &data, const unsigned int &pos);

        /**
         * Write particle data to a vector of doubles.
         *
         * @param [in,out] data The vector of doubles to write integrator data
         * into.
         */
        virtual void write_data(std::vector<double> &data) const;

        /**
         * Set the location of this particle. Note that this does not check
         * whether this is a valid location in the simulation domain.
         *
         * @param [in] new_loc The new location for this particle.
         */
        void
        set_location (const Point<dim> &new_loc);

        /**
         * Get the location of this particle.
         *
         * @return The location of this particle.
         */
        Point<dim>
        get_location () const;

        /**
         * Set the velocity of this particle.
         *
         * @param [in] new_vel The new velocity for this particle.
         */
        void
        set_velocity (Point<dim> new_vel);

        /**
         * Get the velocity of this particle.
         *
         * @return The velocity of this particle.
         */
        Point<dim>
        get_velocity () const;

        /**
         * Get the ID number of this particle.
         *
         * @return The id of this particle.
         */
        double
        get_id () const;

        /**
         * Check whether the particle is marked as being local to this
         * subdomain. Note that this function does not actually perform the
         * check for locality.
         *
         * @return Whether the particle is marked as local.
         */
        bool
        local () const;

        /**
         * Mark the particle as being local of not. Note that this function
         * does not perform the check for locality.
         *
         * @param[in] new_local Whether to mark the particle as local.
         */
        void
        set_local (bool new_local);

        /**
         * Whether to check the particle velocity at its current location.
         * This is used for integrators where the particle velocity may not
         * need to be checked every step.
         *
         * @return Whether to check the particle velocity
         */
        bool
        vel_check () const;

        /**
         * Mark whether to check the particle velocity.
         *
         * @param[in] new_vel_check Whether to check the particle velocity.
         */
        void
        set_vel_check (bool new_vel_check);

        /**
         * Add the MPI data description for this particle type to the vector.
         *
         * @param[in,out] data_info Vector to which MPI data description is
         * appended.
         */
        static void
        add_mpi_types (std::vector<MPIDataInfo> &data_info);
    };

    /**
     * DataParticle provides an example of how to extend the BaseParticle
     * class to include related particle data. This allows users to attach
     * scalars/vectors/tensors/etc to particles and ensure they are
     * transmitted correctly over MPI and written to output files.
     */
    template <int dim, int data_dim>
    class DataParticle : public BaseParticle<dim>
    {
      private:
        double      val[data_dim];

      private:
        DataParticle()
        {
          for (unsigned int i=0; i<data_dim; ++i) val[i] = 0;
        };

        DataParticle(const Point<dim> &new_loc, const double &new_id) : BaseParticle<dim>(new_loc, new_id)
        {
          for (unsigned int i=0; i<data_dim; ++i) val[i] = 0;
        };

        static unsigned int data_len()
        {
          return BaseParticle<dim>::data_len() + data_dim;
        };

        virtual unsigned int read_data(const std::vector<double> &data, const unsigned int &pos)
        {
          unsigned int  p;

          // Read the parent data first
          p = BaseParticle<dim>::read_data(data, pos);

          // Then read the DataParticle data
          for (unsigned i=0; i<data_dim; ++i)
            {
              val[i] = data.at(p++);
            }

          return p;
        };


        virtual void write_data(std::vector<double> &data) const
        {
          // Write the parent data first
          BaseParticle<dim>::write_data(data);

          // Then write the DataParticle data
          for (unsigned i=0; i<data_dim; ++i)
            {
              data.push_back(val[i]);
            }
        };

        /**
         * Returns a vector from the first dim components of val
         *
         * @return vector representation of first dim components of the
         * particle data
         */
        Point<dim>
        get_vector () const;

        /**
         * Sets the first dim components of val to the specified vector value
         *
         * @param [in] new_vec Vector to set the DataParticle data to
         */
        void set_vector(Point<dim> new_vec)
        {
          AssertThrow(data_dim>=dim, std::out_of_range("set_vector"));
          for (unsigned int i=0; i<dim; ++i)
            val[i] = new_vec(i);
        }

        /**
         * Return a reference to an element of the DataParticle data
         *
         * @param [in] ind Index of the data array
         * @return Reference to double value at the requested index
         */
        double &operator[](const unsigned int &ind)
        {
          AssertThrow(data_dim>ind, std::out_of_range("DataParticle[]"));
          return val[ind];
        }

        /**
         * Return the value of an element of the DataParticle data
         *
         * @param [in] ind Index of the data array
         * @return Value at the requested index
         */
        double operator[](const unsigned int &ind) const
        {
          AssertThrow(data_dim>ind, std::out_of_range("DataParticle[]"));
          return val[ind];
        }

        /**
         * Set up the MPI data type information for the DataParticle type
         *
         * @param [in,out] data_info Vector to append MPIDataInfo objects to
         */
        static void add_mpi_types(std::vector<MPIDataInfo> &data_info)
        {
          // Set up the parent types first
          BaseParticle<dim>::add_mpi_types(data_info);

          // Then add our own
          data_info.push_back(MPIDataInfo("data", data_dim));
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
        p (i) = val[i];
    }

  }
}

#endif

