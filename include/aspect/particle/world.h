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

#ifndef __aspect__particle_world_h
#define __aspect__particle_world_h

#include <aspect/global.h>
#include <aspect/particle/particle.h>
#include <aspect/particle/definitions.h>

#include <aspect/particle/integrator/interface.h>
#include <aspect/particle/interpolator/interface.h>
#include <aspect/particle/property/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/simulator_signals.h>

namespace aspect
{
  namespace Particle
  {
    using namespace dealii;

    template <int dim>
    class World : public SimulatorAccess<dim>
    {
      public:
        /**
         * Default World constructor.
         */
        World();

        /**
         * Default World destructor.
         */
        ~World();

        enum ParticleLoadBalancing
        {
          no_balancing,
          remove_particles,
          remove_and_add_particles,
          repartition
        };

        /**
         * Initialize the particle world.
         *
         * @param [in] new_integrator The new Integrator scheme for this
         * world.
         * @param [in] new_manager The new property manager for this
         * world.
         * @param [in] max_part_per_cell Threshold for removing
         * particles from cells. Because the
         * particle world currently does not declare own input parameters this
         * is read by the Tracer postprocessor and set in this function.
         */
        void initialize(Integrator::Interface<dim> *new_integrator,
                        Interpolator::Interface<dim> *new_interpolator,
                        Property::Manager<dim> *new_manager,
                        const ParticleLoadBalancing &load_balancing,
                        const unsigned int max_part_per_cell,
                        const unsigned int weight);

        /**
         * Get the particle property manager for this particle world.
         *
         * @param return manager The property manager for this
         * world.
         */
        const Property::Manager<dim> &
        get_manager() const;

        /**
         * Add a particle to this world. If the specified cell does not exist
         * in the local subdomain an exception will be thrown.
         */
        void add_particle(const Particle<dim> &particle,
                          const LevelInd &cell);

        /**
         * Initialize the particle properties.
         */
        void initialize_particles();

        /**
         * Update the particle properties if necessary.
         */
        void update_particles();

        /**
         * Advect the particle positions by one integration step. Needs to be
         * called until integrator->continue() returns false.
         */
        void advect_particles();

        /**
         * Access to particles in this world.
         */
        std::multimap<LevelInd, Particle<dim> > &
        get_particles();

        /**
         * Const access to particles in this world.
         */
        const std::multimap<LevelInd, Particle<dim> > &
        get_particles() const;

        /**
         * Advance particles by the old timestep using the current
         * integration scheme. This accounts for the fact that the tracers
         * are actually still at their old positions and the current timestep
         * length is already updated for the next step at the time this
         * function is called.
         */
        void advance_timestep();

        /**
         * Calculates the global sum of particles over all processes. This is
         * done to monitor the number of particles that have fallen out of the
         * domain.
         *
         * @return Total number of particles in simulation.
         */
        unsigned int get_global_particle_count() const;

        /**
         * This callback function is called by Simulator to allow us to connect
         * to the SimulatorSignals.
         */
        void
        connector_function(aspect::SimulatorSignals<dim> &signals);

        /**
         * Callback function that is called from Simulator before every
         * refinement. Allows registering store_tracers() in the triangulation.
         */
        void
        register_store_callback_function(typename parallel::distributed::Triangulation<dim> &triangulation);

        /**
         * Callback function that is called from Simulator after every
         * refinement. Allows registering load_tracers() in the triangulation.
         */
        void
        register_load_callback_function(typename parallel::distributed::Triangulation<dim> &triangulation);

        /**
         * Called by listener functions from Triangulation for every cell
         * before a refinement step. A weight is attached to every cell
         * depending on the number of contained tracers.
         */
        unsigned int
        cell_weight(const typename parallel::distributed::Triangulation<dim>::cell_iterator &cell,
                    const typename parallel::distributed::Triangulation<dim>::CellStatus status);

        /**
         * Called by listener functions from Triangulation for every cell
         * before a refinement step. All tracers have to be attached to their
         * element to be sent around to the new cell/processes.
         */
        void
        store_tracers(const typename parallel::distributed::Triangulation<dim>::cell_iterator &cell,
                      const typename parallel::distributed::Triangulation<dim>::CellStatus status,
                      void *data);

        /**
         * Called by listener functions after a refinement step. The local map
         * of particles has to be read from the triangulation user_pointer.
         */
        void
        load_tracers(const typename parallel::distributed::Triangulation<dim>::cell_iterator &cell,
                     const typename parallel::distributed::Triangulation<dim>::CellStatus status,
                     const void *data);

        /**
         * Serialize the contents of this class.
         */
        template <class Archive>
        void serialize (Archive &ar, const unsigned int version);

      private:
        /**
         * Integration scheme for moving particles in this world
         */
        Integrator::Interface<dim>     *integrator;

        /**
         * Integration scheme for moving particles in this world
         */
        Interpolator::Interface<dim>     *interpolator;

        /**
         * The property manager stores information about the additional
         * particle properties and handles the initialization and update of
         * these properties.
         */
        Property::Manager<dim>         *property_manager;

        /**
         * Set of particles currently in the local domain, organized by
         * the level/index of the cell they are in
         */
        std::multimap<LevelInd, Particle<dim> >      particles;

        /**
         * This variable is set by the register_store_callback_function()
         * function and used by the register_load_callback_function() function
         * to check where the tracer data was stored.
         */
        unsigned int data_offset;

        /**
         * Strategy for tracer load balancing.
         */
        ParticleLoadBalancing particle_load_balancing;

        /**
         * Limit for how many particles are allowed per cell. This limit is
         * useful for adaptive meshes to prevent coarse cells from slowing down
         * the whole model. It will only be checked and enforced during mesh
         * refinement and MPI transfer of particles.
         */
        unsigned int max_particles_per_cell;

        /**
         * Weight of a single particle.
         */
        unsigned int tracer_weight;


        /**
         * Returns the number of particles in the cell that contains the
         * most tracers in the global model.
         */
        unsigned int
        get_global_max_tracer_per_cell() const;

        /**
         * Calculate the cells containing each particle for all particles. If
         * particles moved out of the domain of this process they will be send
         * to their new process and inserted there. After this function call
         * every particle is either on its current process and in its current
         * cell, or deleted (if it could not find its new process or cell).
         */
        void find_all_cells();

        /**
         * TODO: This needs to be implemented in case some particles fall out of the
         * domain. In particular for periodic boundary conditions.
         */
        void move_particles_back_in_mesh();

        /**
         * Transfer particles that have crossed subdomain boundaries to other
         * processors. Because subdomains can change drastically during mesh
         * refinement, particle transfer occurs as follows: - Each subdomain
         * finds the particles it owns which have fallen outside it - For each
         * particle outside the subdomain, send the particle to all subdomains
         * and let them determine which one owns it. This assumes there is no
         * overlap between subdomains. - Each process determines which of the
         * received particles is in its subdomain, keeps these and deletes the
         * others
         *
         * @param [in,out] send_particles All particles that should be send
         * are in this vector.
         */
        void send_recv_particles(const std::multimap<types::subdomain_id,Particle <dim> > &send_particles);

        /**
         * Initialize the particle properties of one cell.
         */
        void
        local_initialize_particles(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                   const typename std::multimap<LevelInd, Particle<dim> >::iterator &begin_particle,
                                   const typename std::multimap<LevelInd, Particle<dim> >::iterator &end_particle);

        /**
         * Update the particle properties of one cell.
         */
        void
        local_update_particles(const typename DoFHandler<dim>::active_cell_iterator &cell,
                               const typename std::multimap<LevelInd, Particle<dim> >::iterator &begin_particle,
                               const typename std::multimap<LevelInd, Particle<dim> >::iterator &end_particle);

        /**
         * Advect the particles of one cell. Performs only one step for
         * multi-step integrators. Needs to be called until integrator->continue()
         * evaluates to false.
         */
        void
        local_advect_particles(const typename DoFHandler<dim>::active_cell_iterator &cell,
                               const typename std::multimap<LevelInd, Particle<dim> >::iterator &begin_particle,
                               const typename std::multimap<LevelInd, Particle<dim> >::iterator &end_particle);
    };

    /* -------------------------- inline and template functions ---------------------- */

    template <int dim>
    template <class Archive>
    void World<dim>::serialize (Archive &ar, const unsigned int)
    {
      ar &particles
      ;
    }
  }
}

#endif
