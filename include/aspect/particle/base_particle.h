/*
 Copyright (C) 2011 - 2014 by the authors of the ASPECT code.

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

#ifndef __aspect__particle_base_particle_h
#define __aspect__particle_base_particle_h

#include <aspect/global.h>
#include <deal.II/base/point.h>

namespace aspect
{
  namespace Particle
  {
    using namespace dealii;

    /**
     * Base class of particles - represents a particle with position,
     * velocity, an ID number and an unknown number of properties. This class
     * can be extended to include data related to a particle by the property
     * manager
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
         * TODO: Remove, if possible
         */
        Point<dim>      velocity;

        /**
         * Globally unique ID of particle
         * TODO: Integer?
         */
        double          id;

        /**
         * Whether this particle is in the local subdomain or not. This variable
         * is only needed during that part of the algorithm that transfers
         * particles between processors.
         */
        bool            is_local;

        /**
         * Whether to check the velocity of this particle This is used for
         * integration schemes which require multiple integration steps for
         * some particles, but not for others
         * TODO: this is not needed and should be stored in the integrator
         */
        bool            check_vel;

        /**
         * The serialized vector of all tracer properties
         */
        std::vector<double>      val;


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
        unsigned int
        data_len () const;

        /**
          * Get the number of doubles required to represent this particle for
          * communication.
          *
          * @return Number of doubles required to represent this particle
          */
        void
        set_data_len (const unsigned int data_len);

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
         * Set the properties of this particle.
         *
         * @param [in] new_properties The new properties for this particle.
         */
        void
        set_properties (std::vector<double> new_properties);

        /**
         * Get the properties of this particle.
         *
         * @return The properties of this particle.
         */
        const
        std::vector<double>
        get_properties () const;

        /**
         * Get write-access to properties of this particle.
         *
         * @return The properties of this particle.
         */
        std::vector<double>&
        get_properties ();

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
    };

  }
}

#endif

