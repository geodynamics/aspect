/*
 Copyright (C) 2012 - 2016 by the authors of the ASPECT code.

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

#ifndef _aspect_particle_world_h
#define _aspect_particle_world_h

#include <aspect/global.h>
#include <aspect/particle/particle.h>
#include <aspect/particle/particle_accessor.h>
#include <aspect/particle/particle_iterator.h>
#include <aspect/particle/particle_handler.h>

#include <aspect/particle/generator/interface.h>
#include <aspect/particle/integrator/interface.h>
#include <aspect/particle/interpolator/interface.h>
#include <aspect/particle/property/interface.h>
#include <aspect/particle/property_pool.h>
#include <aspect/particle/output/interface.h>

#include <aspect/simulator_access.h>
#include <aspect/simulator_signals.h>

#include <deal.II/base/timer.h>
#include <deal.II/base/array_view.h>

#include <boost/serialization/unique_ptr.hpp>

namespace aspect
{
  namespace Particle
  {
    using namespace dealii;

    /**
     * This class manages the storage and handling of particles. It provides
     * interfaces to generate and store particles, functions to initialize,
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
         * A type that can be used to iterate over all particles in the domain.
         */
        typedef ParticleIterator<dim> particle_iterator;

        /**
         * Default World constructor.
         */
        World();

        /**
         * Default World destructor.
         */
        ~World();

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
         * Get the particle handler for this particle world.
         *
         * @return The particle handler for this world.
         */
        const ParticleHandler<dim> &
        get_particle_handler() const;

        /**
         * Do initial logic for handling pre-refinement steps
         */
        void setup_initial_state ();

        /**
         * Get the particle interpolator for this particle world.
         *
         * @return The interpolator for this world.
         */
        const Interpolator::Interface<dim> &
        get_interpolator() const;

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
         * Const access to ghost particles in this world.
         * Ghost particles are all particles that are owned by another process
         * and live in one of the ghost cells of the local subdomain.
         */
        const std::multimap<types::LevelInd, Particle<dim> > &
        get_ghost_particles() const;

        /**
         * Advance particles by the old timestep using the current
         * integration scheme. This accounts for the fact that the particles
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
         * refinement and when writing checkpoints.
         * Allows registering store_particles() in the triangulation.
         */
        void
        register_store_callback_function(const bool serialization,
                                         typename parallel::distributed::Triangulation<dim> &triangulation);

        /**
         * Callback function that is called from Simulator after every
         * refinement and after resuming from a checkpoint.
         * Allows registering load_particles() in the triangulation.
         */
        void
        register_load_callback_function(const bool serialization,
                                        typename parallel::distributed::Triangulation<dim> &triangulation);

        /**
         * Called by listener functions from Triangulation for every cell
         * before a refinement step. A weight is attached to every cell
         * depending on the number of contained particles.
         */
        unsigned int
        cell_weight(const typename parallel::distributed::Triangulation<dim>::cell_iterator &cell,
                    const typename parallel::distributed::Triangulation<dim>::CellStatus status);

        /**
         * Called by listener functions from Triangulation for every cell
         * before a refinement step. All particles have to be attached to their
         * element to be sent around to the new cell/processes.
         */
        void
        store_particles(const typename parallel::distributed::Triangulation<dim>::cell_iterator &cell,
                        const typename parallel::distributed::Triangulation<dim>::CellStatus status,
                        void *data);

        /**
         * Called by listener functions after a refinement step. The local map
         * of particles has to be read from the triangulation user_pointer.
         */
        void
        load_particles(const typename parallel::distributed::Triangulation<dim>::cell_iterator &cell,
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
        struct ParticleLoadBalancing
        {
          enum Kind
          {
            no_balancing = 0x0,
            remove_particles = 0x1,
            add_particles = 0x2,
            repartition = 0x4,
            remove_and_add_particles = remove_particles | add_particles
          };
        };

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
         * Particle handler object that is responsible for storing and
         * managing the internal particle structures.
         */
        std_cxx11::unique_ptr<ParticleHandler<dim> > particle_handler;

        /**
         * Set of particles currently in the ghost cells of the local domain,
         * organized by the level/index of the cell they are in. These
         * particles are marked read-only.
         */
        std::multimap<types::LevelInd, Particle<dim> > ghost_particles;

        /**
         * This variable is set by the register_store_callback_function()
         * function and used by the register_load_callback_function() function
         * to check where the particle data was stored.
         */
        unsigned int data_offset;

        /**
         * Strategy for particle load balancing.
         */
        typename ParticleLoadBalancing::Kind particle_load_balancing;

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
         * particle load balancing strategy 'repartition' is used. This value
         * determines how costly the computation of a single particle is compared
         * to the computation of a whole cell, which is arbitrarily defined
         * to represent a cost of 1000.
         */
        unsigned int particle_weight;

        /**
         * Some particle interpolation algorithms require knowledge
         * about particles in neighboring cells. To allow this,
         * particles in ghost cells need to be exchanged between the
         * processes neighboring this cell. This parameter determines
         * whether this transport is happening.
         */
        bool update_ghost_particles;

        /**
         * Get a map between subdomain id and the neighbor index. In other words
         * the returned map answers the question: Given a subdomain id, which
         * neighbor of the current processor's domain (in terms of a contiguous
         * number from 0 to n_neighbors) owns this subdomain?
         */
        std::map<types::subdomain_id, unsigned int>
        get_subdomain_id_to_neighbor_map() const;

        /**
         * Exchanges all particles that live in cells adjacent to ghost cells
         * (i.e. cells that are ghosts to other processes) with the neighboring
         * domains. Clears and re-populates the ghost_neighbors member variable.
         */
        void
        exchange_ghost_particles();

        /**
         * Returns a vector that contains a tensor for every vertex-cell
         * combination of the output of dealii::GridTools::vertex_to_cell_map()
         * (which is expected as input parameter for this function).
         * Each tensor represents a geometric vector from the vertex to the
         * respective cell center.
         */
        std::vector<std::vector<Tensor<1,dim> > >
        vertex_to_cell_centers_directions(const std::vector<std::set<typename parallel::distributed::Triangulation<dim>::active_cell_iterator> > &vertex_to_cells) const;

        /**
         * Finds the cells containing each particle for all particles in
         * @p particles_to_sort. If particles moved out of the local subdomain
         * they will be sent to their new process and inserted there.
         * After this function call every particle is either on its current
         * process and in its current cell, or deleted (if it could not find
         * its new process or cell).
         *
         * @param [in] particles_to_sort Vector containing all pairs of
         * particles and their old cells that will be sorted into the
         * 'particles' member variable in this function.
         */
        void
        sort_particles_in_subdomains_and_cells(const std::vector<std::pair<types::LevelInd, Particle<dim> > > &particles_to_sort);

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
        move_particles_back_into_mesh(const std::vector<std::pair<types::LevelInd, Particle<dim> > >         &lost_particles,
                                      std::vector<std::pair<types::LevelInd, Particle<dim> > >               &moved_particles_cell,
                                      std::vector<std::vector<std::pair<types::LevelInd, Particle<dim> > > > &moved_particles_domain);

        /**
         * Transfer particles that have crossed subdomain boundaries to other
         * processors. The transfer occurs in two steps. As a first step all
         * processes notify their neighbor processes how many particles will
         * be sent to them. Because neighbor processes are defined as the owner
         * of ghost cells of the current process, this also handles
         * periodic boundaries correctly. Afterwards the transfer is done in the
         * same way as local communication between neighbor processes.
         * All received particles and their new cells will be appended to the
         * @p received_particles vector.
         *
         * @param [in] sent_particles All particles that should be sent and
         * their new subdomain_ids are in this map.
         *
         * @param [in,out] received_particles Vector that stores all received
         * particles. Note that it is not required nor checked that the list
         * is empty, received particles are simply attached to the end of
         * the vector.
         */
        void
        send_recv_particles(const std::vector<std::vector<std::pair<types::LevelInd,Particle <dim> > > > &sent_particles,
                            std::vector<std::pair<types::LevelInd, Particle<dim> > >                     &received_particles);

        /**
         * Advect the particle positions by one integration step. Needs to be
         * called until integrator->continue() returns false.
         */
        void advect_particles();

        /**
         * Initialize the particle properties of one cell.
         */
        void
        local_initialize_particles(const typename std::multimap<types::LevelInd, Particle<dim> >::iterator &begin_particle,
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
         * evaluates to false. Particles that moved out of their old cell
         * during this advection step are removed from the local multimap and
         * stored in @p particles_out_of_cell for further treatment (sorting
         * them into the new cell).
         */
        void
        local_advect_particles(const typename DoFHandler<dim>::active_cell_iterator &cell,
                               const typename std::multimap<types::LevelInd, Particle<dim> >::iterator &begin_particle,
                               const typename std::multimap<types::LevelInd, Particle<dim> >::iterator &end_particle,
                               std::vector<std::pair<types::LevelInd, Particle <dim> > >               &particles_out_of_cell);
    };

    /* -------------------------- inline and template functions ---------------------- */

    template <int dim>
    template <class Archive>
    void World<dim>::serialize (Archive &ar, const unsigned int)
    {
      ar
      &particle_handler
      ;

      // Note that initialize is not necessary when saving the particle handler,
      // but it does not harm either. When loading the triangulation we need to
      // recreate the links to the triangulation and the MPI communicator.
      particle_handler->initialize(this->get_triangulation(),this->get_mpi_communicator());
    }
  }
}

#endif
