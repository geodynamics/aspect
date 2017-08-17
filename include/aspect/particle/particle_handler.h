/*
 Copyright (C) 2017 by the authors of the ASPECT code.

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
 along with ASPECT; see the file LICENSE.  If not see
 <http://www.gnu.org/licenses/>.
 */

#ifndef _aspect_particle_particle_handler_h
#define _aspect_particle_particle_handler_h

#include <aspect/global.h>
#include <aspect/particle/particle.h>
#include <aspect/particle/particle_accessor.h>
#include <aspect/particle/particle_iterator.h>

#include <aspect/particle/property_pool.h>

#include <deal.II/base/subscriptor.h>
#include <deal.II/base/array_view.h>
#include <deal.II/base/smartpointer.h>

#include <boost/serialization/map.hpp>
#include <boost/range/iterator_range.hpp>

namespace aspect
{
  namespace Particle
  {
    using namespace dealii;

    template <int> class World;

    /**
     * This class manages the storage and handling of particles. It provides
     * the data structures necessary to store particles efficiently, accessor
     * functions to iterate over particles and find particles, and algorithms
     * to distribute particles in parallel domains.
     *
     * @ingroup Particle
     */
    template <int dim, int spacedim=dim>
    class ParticleHandler: public Subscriptor
    {
      public:
        /**
         * A type that can be used to iterate over all particles in the domain.
         */
        typedef ParticleIterator<dim,spacedim> particle_iterator;

        /**
         * Default constructor.
         */
        ParticleHandler();

        /**
         * Constructor that initializes the particle handler with respect to
         * a given triangulation and MPI communicator. Pointers to the
         * triangulation and the communicator are stored inside of the particle
         *
         */
        ParticleHandler(const parallel::distributed::Triangulation<dim,spacedim> &tria,
                        const MPI_Comm mpi_communicator);

        /**
         * Destructor.
         */
        ~ParticleHandler();

        /**
         * Initialize the particle handler. This function does not clear the
         * internal data structures, it just sets the connections to the
         * MPI communicator and the triangulation. This is useful after de-
         * serialization of a particle handler.
         */
        void initialize(const parallel::distributed::Triangulation<dim,spacedim> &tria,
                        const MPI_Comm mpi_communicator);

        /**
         * Clear all particle related data.
         */
        void clear();

        /**
         * Return an iterator to the first particle.
         */
        ParticleHandler<dim,spacedim>::particle_iterator begin() const;

        /**
         * Return an iterator to the first particle.
         */
        particle_iterator begin();

        /**
         * Return an iterator past the end of the particles.
         */
        particle_iterator end() const;

        /**
         * Return an iterator past the end of the particles.
         */
        particle_iterator end();

        /**
         * Return a pair of particle iterators that mark the begin and end of
         * the particles in a particular cell. The last iterator is the first
         * particle that is no longer in the cell.
         */
        boost::iterator_range<particle_iterator>
        particle_range_in_cell(const typename parallel::distributed::Triangulation<dim>::active_cell_iterator &cell);

        /**
         * Remove a particle pointed to by the iterator.
         */
        void
        remove_particle(const particle_iterator &particle);

        /**
         * Return the total number of particles that were managed by this class
         * the last time the update_n_global_particles() function was called.
         * The actual number of particles may have changed since then if
         * particles have been added or removed.
         *
         * @return Total number of particles in simulation.
         */
        types::particle_index n_global_particles() const;

        /**
         * Return the number of particles in the local part of the
         * triangulation.
         */
        types::particle_index n_locally_owned_particles() const;

        /**
         * Return the number of particles in the given cell.
         */
        unsigned int n_particles_in_cell(const typename Triangulation<dim>::active_cell_iterator &cell) const;

        /**
         * Serialize the contents of this class.
         */
        template <class Archive>
        void serialize (Archive &ar, const unsigned int version);

      protected:
        /**
         * Access to particles in this handler.
         * TODO: This function eventually needs to go, to not expose internal structure.
         * This can only be done once World no longer uses this function.
         */
        std::multimap<types::LevelInd, Particle<dim,spacedim> > &
        get_particles();

        /**
         * Const access to particles in this world.
         * TODO: This function needs to go to not expose internal structure.
         * This can only be done once World no longer uses this function.
         */
        const std::multimap<types::LevelInd, Particle<dim,spacedim> > &
        get_particles() const;

      private:
        /**
         * Address of the triangulation to work on.
         */
        SmartPointer<const parallel::distributed::Triangulation<dim,spacedim>,ParticleHandler<dim,spacedim> > triangulation;

        /**
         * MPI communicator.
         */
        MPI_Comm mpi_communicator;

        /**
         * Set of particles currently in the local domain, organized by
         * the level/index of the cell they are in.
         */
        std::multimap<types::LevelInd, Particle<dim,spacedim> > particles;

        /**
         * This variable stores how many particles are stored globally. It is
         * calculated by update_n_global_particles().
         */
        types::particle_index global_number_of_particles;

        /**
         * The maximum number of particles per cell in the global domain. This
         * variable is important to store and load particle data during
         * repartition and serialization of the solution. Note that the
         * variable is only updated when it is needed, e.g. after particle
         * movement, before/after mesh refinement, before creating a
         * checkpoint and after resuming from a checkpoint.
         */
        unsigned int global_max_particles_per_cell;

        /**
         * This variable stores the next free particle index that is available
         * globally in case new particles need to be generated.
         */
        types::particle_index next_free_particle_index;

        /**
         * Calculates the number of particles in the global model domain.
         */
        void
        update_n_global_particles();

        /**
         * Calculates and stores the number of particles in the cell that
         * contains the most particles in the global model (stored in the
         * member variable global_max_particles_per_cell). This variable is a
         * state variable, because it is needed to serialize and deserialize
         * the particle data correctly in parallel (it determines the size of
         * the data chunks per cell that are stored and read). Before accessing
         * the variable this function has to be called, unless the state was
         * read from another source (e.g. after resuming from a checkpoint).
         */
        void
        update_global_max_particles_per_cell();

        /**
         * Calculates the next free particle index in the global model domain.
         * This equals one plus the highest particle index currently active.
         */
        void
        update_next_free_particle_index();

        /**
         * Make World a friend to access private functions while we transition
         * functionality from World to ParticleHandler.
         * TODO: remove when done moving functionality.
         */
        template <int> friend class World;
    };

    /* -------------------------- inline and template functions ---------------------- */

    template <int dim, int spacedim>
    template <class Archive>
    void ParticleHandler<dim,spacedim>::serialize (Archive &ar, const unsigned int)
    {
      // Note that we do not serialize the particle data itself. Instead we
      // use the serialization functionality of the triangulation class, because
      // this guarantees that data is immediately shipped to new processes if
      // the domain is distributed differently after resuming from a checkpoint.
      ar //&particles
      &global_number_of_particles
      &global_max_particles_per_cell
      &next_free_particle_index;
    }
  }
}

#endif
