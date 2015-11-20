/*
 Copyright (C) 2015 by the authors of the ASPECT code.

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
#include <deal.II/base/types.h>

#include <boost/serialization/vector.hpp>

namespace aspect
{
  namespace Particle
  {
    using namespace dealii;
    /**
     * A namespace for all type definitions related to particles.
     */
    namespace types
    {
      using namespace dealii::types;

      /**
       * Typedef of cell level/index pair. TODO: replace this by the
       * active_cell_index from deal.II 8.3 onwards.
       */
      typedef std::pair<int, int> LevelInd;

      /* Type definitions */

#ifdef DEAL_II_WITH_64BIT_INDICES
      /**
       * The type used for indices of tracers. While in
       * sequential computations the 4 billion indices of 32-bit unsigned integers
       * is plenty, parallel computations using hundreds of processes can overflow
       * this number and we need a bigger index space. We here utilize the same
       * build variable that controls the dof indices of deal.II because the number
       * of degrees of freedom and the number of tracers are typically on the same
       * order of magnitude.
       *
       * The data type always indicates an unsigned integer type.
       */
      typedef unsigned long long int particle_index;

      /**
       * An identifier that denotes the MPI type associated with
       * types::global_dof_index.
       */
#  define ASPECT_TRACER_INDEX_MPI_TYPE MPI_UNSIGNED_LONG_LONG
#else
      /**
       * The type used for indices of tracers. While in
       * sequential computations the 4 billion indices of 32-bit unsigned integers
       * is plenty, parallel computations using hundreds of processes can overflow
       * this number and we need a bigger index space. We here utilize the same
       * build variable that controls the dof indices of deal.II because the number
       * of degrees of freedom and the number of tracers are typically on the same
       * order of magnitude.
       *
       * The data type always indicates an unsigned integer type.
       */
      typedef unsigned int particle_index;

      /**
       * An identifier that denotes the MPI type associated with
       * types::global_dof_index.
       */
#  define ASPECT_TRACER_INDEX_MPI_TYPE MPI_UNSIGNED
#endif
    }

    /**
     * Base class of particles - represents a particle with position,
     * an ID number and a variable number of properties. This class
     * can be extended to include data related to a particle by the property
     * manager.
     *
     * @ingroup Particle
     *
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
                  const types::particle_index &new_id);

        /**
         * Constructor for Particle, creates a particle from a data vector.
         * This constructor is usually called after sending a particle to a
         * different process.
         *
         * @param[in,out] begin_data A pointer to a memory location from which
         * to read the information that completely describes a particle.
         * This class then de-serializes its data from this memory location,
         * using @p data_size as the length of the memory block from which to
         * read the data. @p begin_data is advanced by @p data_size bytes in
         * this constructor.

         * @param[in] data_size Size in bytes of the begin_data array
         * that will be read in by this particle.
         */
        Particle (const void *&begin_data,
                  const unsigned int data_size);

        /**
         * Destructor for Particle
         */
        virtual
        ~Particle ();

        /**
          * Resize the properties member variable to hold the number of doubles
          * required to represent this particle's properties. This function
          * is called by the PropertyManager class to prepare this particle
          * for the initialization of its properties.
          *
          * @param [in] n_components Number of additional property components
          * that this particle holds.
          */
        void
        set_n_property_components (const unsigned int n_components);

        /**
         * Write particle data into a data array. The array is expected
         * to be large enough to take the data, and the void pointer should
         * point to the first element in which the data should be written. This
         * function is meant for serializing all particle properties and
         * afterwards de-serializing the properties by calling the appropriate
         * constructor Particle(void *&data, const unsigned int data_size);
         *
         * @param [in,out] data The memory location to write particle data
         * into. This pointer points to the begin of the memory, in which the
         * data will be written and it will be advanced by the serialized size
         * of this particle.
         */
        void
        write_data(void *&data) const;

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
        types::particle_index
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
        Point<dim>             location;

        /**
         * Globally unique ID of particle
         */
        types::particle_index  id;

        /**
         * The vector of all tracer properties
         */
        std::vector<double>    properties;
    };

    /* -------------------------- inline and template functions ---------------------- */

    template <int dim>
    template <class Archive>
    void Particle<dim>::serialize (Archive &ar, const unsigned int)
    {
      ar &location
      & id
      & properties
      ;
    }
  }
}

#endif

