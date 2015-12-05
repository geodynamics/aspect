/*
 Copyright (C) 2012 - 2015 by the authors of the ASPECT code.

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

#include <aspect/particle/generator/interface.h>
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

    /**
     * This class manages the storage and handling of particles. It provides
     * interfaces to generate and store tracers, functions to initialize,
     * update and advect them, and ways to retrieve information about the
     * particles. The implementation of most of these methods is outsourced
     * to different plugin systems, this class is mostly concerned with
     * managing the interactions of the different systems with the code
     * outside the particle world.
     *
     * @ingroup Particle
     */
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
         * Initialize the particle world. Most of the arguments are read by the
         * Tracer postprocessor and handed to this function, because the
         * particle world currently does not declare own input parameters.
         *
         * @param [in] generator The particle generator for this world.
         * @param [in] integrator The Integrator scheme for this world.
         * @param [in] interpolator The Interpolator scheme for this
         * world. This object defines how particle properties are interpolated
         * to arbitrary positions within the domain. An example could be to
         * take the properties of the closest particle, but more complicated
         * schemes are possible.
         * @param [in] manager The property manager for this world.
         * @param [in] load_balancing The strategy used to balance the
         * computational load for particle advection across processes.
         * @param [in] max_part_per_cell Threshold for removing
         * particles from cells.
         * @param [in] weight The computational load that is associated with
         * integrating and updating one particle.
         *
         * @note World takes ownership of the @p generator, @p integrator, @p
         * interpolator, and @p manager object in this function.
         *
         */
        void initialize(Generator::Interface<dim> *generator,
                        Integrator::Interface<dim> *integrator,
                        Interpolator::Interface<dim> *interpolator,
                        Property::Manager<dim> *manager,
                        const ParticleLoadBalancing &load_balancing,
                        const unsigned int max_part_per_cell,
                        const unsigned int weight);

        /**
         * Get the particle property manager for this particle world.
         *
         * @return The property manager for this world.
         */
        const Property::Manager<dim> &
        get_property_manager() const;

        /**
         * Initialize the particle properties.
         */
        void generate_particles();
        /**
         * Initialize the particle properties.
         */
        void initialize_particles();

        /**
         * Access to particles in this world.
         */
        std::multimap<types::LevelInd, Particle<dim> > &
        get_particles();

        /**
         * Const access to particles in this world.
         */
        const std::multimap<types::LevelInd, Particle<dim> > &
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
         * Return the total number of particles in the simulation. This
         * function is useful for monitoring how many particles have been
         * lost by falling out of the domain. Not that this function does
         * not compute the number of particles, because that is an expensive
         * global MPI operation. Instead it returns the number, which is
         * updated internally every time it might change by a call to
         * update_n_global_particles().
         *
         * @return Total number of particles in simulation.
         */
        types::particle_index n_global_particles() const;

        /**
         * This callback function is registered within Simulator by the
         * constructor of this class and will be
         * called from Simulator during construction. It allows to attach slot
         * functions to the provided SimulatorSignals. This world will register
         * the register_store_callback_function() and
         * register_load_callback_function() to the
         * pre_refinement_store_user_data signal and the
         * post_refinement_load_user_data signal respectively.
         */
        void
        connect_to_signals(aspect::SimulatorSignals<dim> &signals);

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
         * Update the particle properties if necessary.
         */
        void update_particles();

        /**
         * Serialize the contents of this class.
         */
        template <class Archive>
        void serialize (Archive &ar, const unsigned int version);

      private:
        /**
         * Generation scheme for creating particles in this world
         */
        std_cxx11::unique_ptr<Generator::Interface<dim> > generator;

        /**
         * Integration scheme for moving particles in this world
         */
        std_cxx11::unique_ptr<Integrator::Interface<dim> > integrator;

        /**
         * Integration scheme for moving particles in this world
         */
        std_cxx11::unique_ptr<Interpolator::Interface<dim> > interpolator;

        /**
         * The property manager stores information about the additional
         * particle properties and handles the initialization and update of
         * these properties.
         */
        std_cxx11::unique_ptr<Property::Manager<dim> > property_manager;

        /**
         * Set of particles currently in the local domain, organized by
         * the level/index of the cell they are in.
         */
        std::multimap<types::LevelInd, Particle<dim> > particles;

        /**
         * This variable stores how many particles are stored globally. It is
         * calculated by update_n_global_particles().
         */
        types::particle_index global_number_of_particles;

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
         * The computational cost of a single particle. This is an input
         * parameter that is set during initialization and is only used if the
         * particle_load_balancing strategy 'repartition' is used. This value
         * determines how costly the computation of a single tracer is compared
         * to the computation of a whole cell, which is arbitrarily defined
         * to represent a cost of 1000.
         */
        unsigned int tracer_weight;

        /**
         * Calculates the number of particles in the global model domain.
         */
        void
        update_n_global_particles();

        /**
         * Returns the number of particles in the cell that contains the
         * most tracers in the global model.
         */
        unsigned int
        get_global_max_tracers_per_cell() const;

        /**
         * Finds the cells containing each particle for all particles. If
         * particles moved out of this subdomain they will be sent
         * to their new process and inserted there. After this function call
         * every particle is either on its current process and in its current
         * cell, or deleted (if it could not find its new process or cell).
         */
        void
        sort_particles_in_subdomains_and_cells();

        /**
         * TODO: Implement this for arbitrary meshes.
         * This function checks if the @p lost_particles moved across a
         * periodic boundary and tries to reinsert them into
         * @p moved_particles_cell or @p moved_particles_domain. All particles
         * that can not be found are discarded.
         */
        void
        move_particles_back_into_mesh(std::multimap<types::LevelInd, Particle<dim> >            &lost_particles,
                                      std::multimap<types::LevelInd, Particle<dim> >            &moved_particles_cell,
                                      std::multimap<types::subdomain_id, Particle<dim> >        &moved_particles_domain);

        /**
         * Transfer particles that have crossed subdomain boundaries to other
         * processors. The transfer occurs in two steps. As a first step all
         * processes notify their neighbor processes how many particles will
         * be sent to them. Because neighbor processes are defined as the owner
         * of ghost cells of the current process, this also handles
         * periodic boundaries correctly. Afterwards the transfer is done in the
         * same way as local communication between neighbor processes.
         * All received particles will immediately be inserted into the
         * particles member variable.
         *
         * @param [in] sent_particles All particles that should be sent and
         * their new subdomain_ids are in this map.
         */
        void
        send_recv_particles(const std::multimap<types::subdomain_id,Particle <dim> > &sent_particles);

        /**
         * Advect the particle positions by one integration step. Needs to be
         * called until integrator->continue() returns false.
         */
        void advect_particles();

        /**
         * Initialize the particle properties of one cell.
         */
        void
        local_initialize_particles(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                   const typename std::multimap<types::LevelInd, Particle<dim> >::iterator &begin_particle,
                                   const typename std::multimap<types::LevelInd, Particle<dim> >::iterator &end_particle);

        /**
         * Update the particle properties of one cell.
         */
        void
        local_update_particles(const typename DoFHandler<dim>::active_cell_iterator &cell,
                               const typename std::multimap<types::LevelInd, Particle<dim> >::iterator &begin_particle,
                               const typename std::multimap<types::LevelInd, Particle<dim> >::iterator &end_particle);

        /**
         * Advect the particles of one cell. Performs only one step for
         * multi-step integrators. Needs to be called until integrator->continue()
         * evaluates to false.
         */
        void
        local_advect_particles(const typename DoFHandler<dim>::active_cell_iterator &cell,
                               const typename std::multimap<types::LevelInd, Particle<dim> >::iterator &begin_particle,
                               const typename std::multimap<types::LevelInd, Particle<dim> >::iterator &end_particle);
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
