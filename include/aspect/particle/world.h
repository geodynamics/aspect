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
#include <aspect/particle/output/interface.h>

#include <aspect/simulator_access.h>
#include <aspect/simulator_signals.h>

#include <deal.II/base/timer.h>

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
         * Initialize the particle world.
         */
        void initialize();

        /**
         * Get the particle property manager for this particle world.
         *
         * @return The property manager for this world.
         */
        const Property::Manager<dim> &
        get_property_manager() const;

        /**
         * Do initial logic for handling pre-refinement steps
         */
        void setup_initial_state ();

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
         * Returns the number of particles in the cell that contains the
         * most tracers in the global model.
         */
        unsigned int
        get_global_max_tracers_per_cell() const;

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
         * Generate the selected particle output.
         */
        std::string
        generate_output() const;

        /**
         * Serialize the contents of this class.
         */
        template <class Archive>
        void serialize (Archive &ar, const unsigned int version);

        /**
         * Save the state of the object.
         */
        virtual
        void
        save (std::ostringstream &os) const;

        /**
         * Restore the state of the object.
         */
        virtual
        void
        load (std::istringstream &is);

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);

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
         * Pointer to an output object
         */
        std_cxx11::unique_ptr<Output::Interface<dim> > output;

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
         * This variable stores the next free particle index that is available
         * globally in case new particles need to be generated.
         */
        types::particle_index next_free_particle_index;

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
         * Lower limit for particle number per cell. This limit is
         * useful for adaptive meshes to prevent fine cells from being empty
         * of particles. It will be checked and enforced after mesh
         * refinement and after particle movement. If there are
         * n_number_of_particles < min_particles_per_cell
         * particles in one cell then
         * min_particles_per_cell - n_number_of_particles particles are
         * generated and randomly placed in this cell. If the particles carry
         * properties the individual property plugins control how the
         * properties of the new particles are initialized.
         */
        unsigned int min_particles_per_cell;

        /**
         * Upper limit for particle number per cell. This limit is
         * useful for adaptive meshes to prevent coarse cells from slowing down
         * the whole model. It will be checked and enforced after mesh
         * refinement, after MPI transfer of particles and after particle
         * movement. If there are
         * n_number_of_particles > max_particles_per_cell
         * particles in one cell then
         * n_number_of_particles - max_particles_per_cell
         * particles in this cell are randomly chosen and destroyed.
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
         * Calculates the next free particle index in the global model domain.
         * This equals one plus the highest particle index currently active.
         */
        void
        update_next_free_particle_index();

        /**
         * Returns whether a given particle is in the given cell.
         */
        bool
        particle_is_in_cell(const Particle<dim> &particle,
                            const typename parallel::distributed::Triangulation<dim>::active_cell_iterator &cell) const;

        /**
         * Returns a map of neighbor cells of the current cell. This map is
         * sorted according to the distance between the particle and the face
         * of cell that is shared with the neighbor cell. I.e. the first
         * entries of the map are the most likely ones to find the particle in.
         */
        std::multimap<double, typename parallel::distributed::Triangulation<dim>::active_cell_iterator>
        neighbor_cells_to_search(const Particle<dim> &particle,
                                 const typename parallel::distributed::Triangulation<dim>::active_cell_iterator &cell) const;

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
         * Apply the bounds for the maximum and minimum number of particles
         * per cell, if the appropriate @p particle_load_balancing strategy
         * has been selected.
         */
        void
        apply_particle_per_cell_bounds();

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
      &next_free_particle_index
      ;
    }
  }
}

#endif
