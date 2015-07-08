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

#include <aspect/particle/particle.h>
#include <aspect/particle/definitions.h>

#include <aspect/particle/integrator/interface.h>
#include <aspect/particle/property/interface.h>
#include <aspect/simulator_access.h>

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
         * Default World destructor, deallocates all relevant arrays and
         * structures.
         */
        ~World();

        /**
         * Initialize the particle world by connecting to be informed when
         * the triangulation changes.
         */
        void init();

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
         * Get the particle property manager for this particle world.
         *
         * @param return manager The new property manager for this
         * world.
         */
        const Property::Manager<dim> &
        get_manager() const;

        /**
         * Add a particle to this world. If the specified cell does not exist
         * in the local subdomain an exception will be thrown.
         */
        void add_particle(const Particle<dim> &particle, const LevelInd &cell);

        /**
         * Initialize the particle properties.
         */
        void initialize_particles();

        /**
         * Update the particle properties if necessary.
         */
        void update_particles();

        /**
         * Access to particles in this world.
         */
        std::multimap<LevelInd, Particle<dim> > &get_particles();

        /**
         * Const access to particles in this world.
         */
        const std::multimap<LevelInd, Particle<dim> > &get_particles() const;

        /**
         * Get the names and number of particle properties from the
         * property_manager.
         *
         * @param [inout] names Vector of property names attached to particles
         * @param [inout] length Number of doubles needed to represent properties
         */
        void
        get_data_info(std::vector<std::string> &names,
                      std::vector<unsigned int> &length) const;

        /**
         * Calculate the cells containing each particle for all particles. If
         * particles moved out of the domain of this process they will be send
         * to their new process and inserted there. After this function call
         * every particle is either on its current process and in its current
         * cell, or deleted (if it could not find its new process or cell).
         */
        void find_all_cells();

        /**
         * Advance particles by the old timestep using the current
         * integration scheme. This accounts for the fact that the tracers
         * are actually still at their old positions and the current timestep
         * length is already updated for the next step at the time this
         * function is called.
         */
        void advance_timestep();

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
         * others - TODO: handle particles outside any domain - TODO: if we
         * know the domain of a particle (e.g. bordering domains), send it
         * only to that domain
         *
         * @param [in,out] send_particles All particles that should be send
         * are in this vector.
         */
        void send_recv_particles(const std::vector<Particle <dim> > &send_particles);

        /**
         * Calculates the current and old velocities for each particle at its
         * current location. The velocities are stored in the vectors given
         * to the function, ordered in the same way an iterator iterates over
         * the map of particles. This means the ordering is only correct as
         * long as the map does not change.
         *
         * @param [inout] velocities The current solution vector for this
         * simulation.
         * @param [inout] old_velocities A vector of velocities evaluated at
         * the particle positions according to the given solution.
         */
        void
        get_particle_velocities(std::vector<Tensor<1,dim> > &velocities,
                                std::vector<Tensor<1,dim> > &old_velocities);

        /**
         * Calculates the global sum of particles over all processes. This is
         * done to monitor the number of particles that have fallen out of the
         * domain.
         *
         * @return Total number of particles in simulation.
         */
        unsigned int get_global_particle_count() const;

        /**
         * Read or write the data of this object for serialization
         */
        template <class Archive>
        void serialize(Archive &ar, const unsigned int version);

      private:
        /**
         * Integration scheme for moving particles in this world
         */
        Integrator::Interface<dim>     *integrator;

        /**
         * The property manager stores information about the additional
         * particle properties and handles the initialization and update of
         * these properties.
         */
        Property::Manager<dim>         *property_manager;

        /**
         * Whether the triangulation was changed (e.g. through refinement), in which
         * case we must treat all recorded particle level/index values as invalid
         */
        bool                            triangulation_changed;

        /**
         * Set of particles currently in the local domain, organized by
         * the level/index of the cell they are in
         */
        std::multimap<LevelInd, Particle<dim> >      particles;

        /**
         * Called by listener functions to indicate that the mesh of this
         * subdomain has changed.
         */
        void mesh_changed();
    };
  }
}

#endif
