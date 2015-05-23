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

#ifndef __aspect__particle_world_h
#define __aspect__particle_world_h

#include <aspect/particle/base_particle.h>
#include <aspect/particle/definitions.h>

#include <aspect/particle/integrator/interface.h>
#include <aspect/particle/property/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace Particle
  {
    using namespace dealii;

    // TODO: in the future, upgrade multimap to ParticleMap typedef
    // with C++11 standard "using" syntax

    template <int dim>
    class World : public SimulatorAccess<dim>
    {
      private:
        /// Integration scheme for moving particles in this world
        Integrator::Interface<dim>     *integrator;

        /// Integration scheme for moving particles in this world
        Property::Manager<dim>         *property_manager;

        /// Whether the triangulation was changed (e.g. through refinement), in which
        /// case we must treat all recorded particle level/index values as invalid
        bool                            triangulation_changed;

        /// Set of particles currently in the local domain, organized by
        /// the level/index of the cell they are in
        std::multimap<LevelInd, BaseParticle<dim> >      particles;

        // Total number of particles in simulation
        unsigned int                    global_num_particles;

        // MPI related variables
        /// MPI registered datatype encapsulating the MPI_Particle struct
        MPI_Datatype                    particle_type;

        /// Size and rank in the MPI communication world
        unsigned int                    world_size;
        unsigned int                    self_rank;

        /*
         * Information about the data of particles
         */
        std::vector<MPIDataInfo>        data_info;

        /**
         * Recursively determines which cell the given particle belongs to.
         * Returns true if the particle is in the specified cell and sets the
         * particle cell information appropriately, false otherwise.
         *
         * @param [in,out] particle The particle for which a cell is being
         * searched for. The particle will be marked to indicate whether it is
         * in the local subdomain or not.
         * @param [in] cur_cell The current cell level and index being
         * investigated as potentially containing the particle.
         * @return The level and index of the cell the particle was determined
         * to be in.  If no cell was found this returns (-1, -1).
         */
        LevelInd recursive_find_cell(BaseParticle<dim> &particle,
                                     const LevelInd cur_cell);

        /**
         * Called by listener functions to indicate that the mesh of this
         * subdomain has changed.
         */
        void mesh_changed();

      public:
        /**
         * Default World constructor.
         */
        World();

        /**
         * Default World destructor, deallocates all relevant arrays and
         * structures.
         */
        ~World();

        /**
         * Set the particle Integrator scheme for this particle world.
         *
         * @param [in] new_integrator The new Integrator scheme for this
         * world.
         */
        void set_integrator(Integrator::Interface<dim> *new_integrator);

        /**
         * Set the particle property manager for this particle world.
         *
         * @param [in] new_manager The new property manager for this
         * world.
         */
        void set_manager(Property::Manager<dim> *new_manager);

        /**
         * All processes must call this function when finished adding
         * particles to the world. This function will determine the total
         * number of particles.
         */
        void finished_adding_particles();

        /**
         * Add a particle to this world. If the specified cell does not exist
         * in the local subdomain an exception will be thrown.
         */
        void add_particle(const BaseParticle<dim> &particle, const LevelInd &cell);

        void initialize_particles();

        void update_particles();

        /**
         * Access to particles in this world.
         */
        std::multimap<LevelInd, BaseParticle<dim> > &get_particles();

        /**
         * Const access to particles in this world.
         */
        const std::multimap<LevelInd, BaseParticle<dim> > &get_particles() const;

        /**
         * Const access to particles in this world.
         */
        const std::vector<MPIDataInfo> &get_mpi_datainfo() const;

        /**
         * Initialize the particle world by creating appropriate MPI data
         * types for transferring particles, and allocating memory for MPI
         * related functions.
         */
        void init();

        /**
         * Calculate the cells containing each particle for all particles.
         */
        void find_all_cells();

        /**
         * Advance particles by the specified timestep using the current
         * integration scheme.
         *
         * @param [in] timestep Length of timestep to integrate particle
         * movement
         * @param [in] solution Current Aspect solution vector
         */
        void advance_timestep(const double timestep, const LinearAlgebra::BlockVector &solution);

        void move_particles_back_in_mesh();

        /**
         * Mark all particles to be checked for velocity at their current
         * position
         */
        void mark_particles_for_check();

        /**
         * Finds the cell the particle is contained in and returns the
         * appropriate cell level/index.
         *
         * @param [in,out] particle The particle to find the cell for. This
         * particle will be updated to indicate whether it is in the local
         * subdomain or not.
         * @param [in] cur_cell The current cell (level and index) being
         * checked.
         * @return The level and index of the active cell the particle is in.
         * If no cell was found to contain the particle, return the
         * level/index (-1, -1)
         */
        LevelInd find_cell(BaseParticle<dim> &particle, const LevelInd &cur_cell);

        /**
         * Transfer particles that have crossed subdomain boundaries to other
         * processors. Because subdomains can change drastically during mesh
         * refinement, particle transfer occurs as follows: - Each subdomain
         * finds the particles it owns which have fallen outside it - For each
         * particle outside the subdomain, send the particle to all subdomains
         * and let them determine which one owns it. This assumes there is no
         * overlap between subdomains. - Each process determines which of the
         * received particles is in its subdomain, keeps these and deletes the
         * others - TODO: handle particles outside any domain - TODO: if we
         * know the domain of a particle (e.g. bordering domains), send it
         * only to that domain
         */
        void send_recv_particles();

        /**
         * Calculates the velocities for each particle at its location given
         * the input solution velocity field. The calculated velocities are
         * stored in the Particle objects for this world.
         *
         * @param [in] solution The current solution vector for this
         * simulation.
         */
        void get_particle_velocities(const LinearAlgebra::BlockVector &solution);

        /**
         * Calculates the global sum of particles over all processes. This is
         * done to ensure no particles have fallen out of the simulation
         * domain.
         *
         * @return Total number of particles in simulation.
         */
        unsigned int get_global_particle_count() const;

        /**
         * Checks that the number of particles in the simulation has not
         * unexpectedly changed. If the particle count changes then the
         * simulation will be aborted.
         */
        void check_particle_count();

        /**
         * Read or write the data of this object for serialization
         */
        template <class Archive>
        void serialize(Archive &ar, const unsigned int version);
    };
  }
}

#endif
