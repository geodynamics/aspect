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

#ifndef __aspect__particle_particle_h
#define __aspect__particle_particle_h

#include <aspect/global.h>

#include <deal.II/base/point.h>

#include <boost/serialization/vector.hpp>

namespace aspect
{
  namespace Particle
  {
    using namespace dealii;

    /**
     * Base class of particles - represents a particle with position,
     * an ID number and an variable number of properties. This class
     * can be extended to include data related to a particle by the property
     * manager.
     */
    template <int dim>
    class Particle
    {
      public:
        /**
         * Empty constructor for Particle, creates a particle at the
         * origin.
         */
        Particle ();

        /**
         * Constructor for Particle, creates a particle with the specified
         * ID at the specified location. Note that Aspect
         * does not check for duplicate particle IDs so the generator must
         * make sure the IDs are unique over all processes.
         *
         * @param[in] new_loc Initial location of particle.
         * @param[in] new_id Globally unique ID number of particle.
         */
        Particle (const Point<dim> &new_loc,
                  const double &new_id);

        /**
         * Constructor for Particle, creates a particle from a data vector.
         * This constructor is usually called after sending a particle to a
         * different process.
         *
         * @param[in,out] begin_data First data component.
         * @param[in] data_len Number of components of the begin_data vector
         * that will be read in by this particle.
         */
        Particle (std::vector<double>::const_iterator &begin_data,
                  const unsigned int data_len);

        /**
         * Destructor for Particle
         */
        virtual
        ~Particle ();

        /**
         * Get the number of doubles required to represent this particle for
         * communication. This includes the base properties like position
         * and id as well as the additional user requested properties.
         *
         * @return Number of doubles required to represent this particle
         */
        unsigned int
        data_len () const;

        /**
          * Set the number of doubles required to represent this particle for
          * communication. This includes the base properties position
          * and id as well as the additional user requested properties.
          *
          * @param [in] data_len Number of doubles to represent this particle
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
        virtual unsigned int read_data(const std::vector<double> &data, const unsigned int pos);

        /**
         * Write particle data to a vector of doubles. The vector is extended
         * by the data written by this particle via push_back().
         *
         * @param [in,out] data The vector of doubles to write integrator data
         * into.
         */
        virtual void write_data(std::vector<double> &data) const;

        /**
         * Write particle data to a vector of doubles. The vector is expected
         * to be large enough to take the data, and the input iterator should
         * point to the first element in which the data should be written. This
         * function avoids the resizing of the previous write_data function.
         *
         * @param [in,out] data The vector of doubles to write integrator data
         * into. This iterator points to the first element, in which the data
         * should be written.
         */
        virtual void write_data(std::vector<double>::iterator &data) const;


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
        const Point<dim> &
        get_location () const;

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
        set_properties (const std::vector<double> &new_properties);

        /**
         * Get write-access to properties of this particle.
         *
         * @return The properties of this particle.
         */
        std::vector<double> &
        get_properties ();

        /**
         * Get the properties of this particle.
         *
         * @return The properties of this particle.
         */
        const std::vector<double> &
        get_properties () const;

        /**
         * Serialize the contents of this class.
         */
        template <class Archive>
        void serialize (Archive &ar, const unsigned int version);

      private:
        /**
         * Current particle location
         */
        Point<dim>      location;

        /**
         * Globally unique ID of particle
         * TODO: Integer?
         */
        double          id;

        /**
         * The serialized vector of all tracer properties
         */
        std::vector<double>      val;
    };

    /* -------------------------- inline and template functions ---------------------- */

    template <int dim>
    template <class Archive>
    void Particle<dim>::serialize (Archive &ar, const unsigned int)
    {
      ar &location
      & id
      & val
      ;
    }
  }
}

#endif

