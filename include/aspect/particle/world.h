/*
 Copyright (C) 2011, 2012, 2013 by the authors of the ASPECT code.

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

#include <deal.II/numerics/fe_field_function.h>
#include <aspect/particle/particle.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace Particle
  {
    // TODO: in the future, upgrade multimap to ParticleMap typedef
    // with C++11 standard "using" syntax

    /// MPI tag for particle transfers
    const int           PARTICLE_XFER_TAG = 382;

    namespace Integrator
    {
      template <int dim, class T>
      class Interface;
    }

    template <int dim, class T>
    class World
    {
      private:
        /// Mapping for the simulation this particle world exists in
        const Mapping<dim>              *mapping;

        /// Triangulation for the simulation this particle world exists in
        const parallel::distributed::Triangulation<dim>   *triangulation;

        /// DoFHandler for the simulation this particle world exists in
        const DoFHandler<dim>           *dof_handler;

        /// DoFHandler for the simulation this particle world exists in
        const LinearAlgebra::BlockVector           *solution;

        /// Integration scheme for moving particles in this world
        Integrator::Interface<dim, T>   *integrator;

        /// MPI communicator to be used for this world
        MPI_Comm                        communicator;

        /// Whether the triangulation was changed (e.g. through refinement), in which
        /// case we must treat all recorded particle level/index values as invalid
        bool                            triangulation_changed;

        /// Set of particles currently in the local domain, organized by
        /// the level/index of the cell they are in
        std::multimap<LevelInd, T>      particles;

        // Total number of particles in simulation
        unsigned int                    global_num_particles;

        // MPI related variables
        /// MPI registered datatype encapsulating the MPI_Particle struct
        MPI_Datatype                    particle_type;

        /// Size and rank in the MPI communication world
        unsigned int                    world_size;
        unsigned int                    self_rank;

        /// Buffers indicating how many particles to send/recv to each process
        int                             *num_send, *num_recv;
        /// Total number of particles to send/recv
        int                             total_send, total_recv;
        /// Send/recv offset into data buffer for each process
        int                             *send_offset, *recv_offset;
        /// MPI_Request object buffers to allow for non-blocking communication
        MPI_Request                     *send_reqs, *recv_reqs;


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
        LevelInd recursive_find_cell(T &particle,
                                     const LevelInd cur_cell)
        {
          typename parallel::distributed::Triangulation<dim>::cell_iterator it, found_cell, child_cell;
          unsigned int    child_num;
          LevelInd        res, child_li;

          // If the particle is in the specified cell
          found_cell = typename parallel::distributed::Triangulation<dim>::cell_iterator(triangulation, cur_cell.first, cur_cell.second);
          if (found_cell != triangulation->end() && found_cell->point_inside(particle.get_location()))
            {
              // If the cell is active, we're at the finest level of refinement and can finish
              if (found_cell->active())
                {
                  particle.set_local(found_cell->is_locally_owned());
                  return cur_cell;
                }
              else
                {
                  // Otherwise we need to search deeper
                  for (child_num=0; child_num<found_cell->n_children(); ++child_num)
                    {
                      child_cell = found_cell->child(child_num);
                      child_li = LevelInd(child_cell->level(), child_cell->index());
                      res = recursive_find_cell(particle, child_li);
                      if (res.first != -1 && res.second != -1) return res;
                    }
                }
            }

          // If we still can't find it, return false
          return LevelInd(-1, -1);
        };

        /**
         * Called by listener functions to indicate that the mesh of this
         * subdomain has changed.
         */
        void mesh_changed()
        {
          triangulation_changed = true;
        };

      public:
        /**
         * Default World constructor.
         */
        World()
        {
          triangulation_changed = true;
          total_send = total_recv = 0;
          world_size = self_rank = 0;
          num_send = num_recv = send_offset = recv_offset = NULL;
          send_reqs = recv_reqs = NULL;
          triangulation = NULL;
          mapping = NULL;
          dof_handler = NULL;
          solution = NULL;
          integrator = NULL;
        };

        /**
         * Default World destructor, deallocates all relevant arrays and
         * structures.
         */
        ~World()
        {
          if (world_size) MPI_Type_free(&particle_type);

          if (num_send) delete [] num_send;
          if (num_recv) delete [] num_recv;
          if (send_offset) delete [] send_offset;
          if (recv_offset) delete [] recv_offset;
          if (send_reqs) delete [] send_reqs;
          if (recv_reqs) delete [] recv_reqs;
        };

        /**
         * Get the deal.II Mapping associated with this particle world.
         *
         * @return The Mapping for this world.
         */
        const Mapping<dim> *get_mapping() const
        {
          return mapping;
        };

        /**
         * Set the deal.II Mapping associated with this particle world.
         *
         * @param [in] new_mapping The new Mapping for this world.
         */
        void set_mapping(const Mapping<dim> *new_mapping)
        {
          mapping = new_mapping;
        };

        /**
         * Set the deal.II Triangulation associated with this particle world
         * and connects relevant listener for mesh changes.
         *
         * @param [in] new_tria The new Triangulation for this world.
         */
        void set_triangulation(const parallel::distributed::Triangulation<dim> *new_tria)
        {
          //if (triangulation) triangulation->signals.post_refinement.disconnect(std_cxx1x::bind(&World::mesh_changed, std_cxx1x::ref(*this)));
          triangulation = new_tria;
          triangulation->signals.post_refinement.connect(std_cxx1x::bind(&World::mesh_changed, std_cxx1x::ref(*this)));
        };

        /**
         * Get the deal.II Triangulation associated with this particle world.
         *
         * @return const pointer to associated Triangulation
         */
        const parallel::distributed::Triangulation<dim> *get_triangulation()
        {
          return triangulation;
        };

        /**
         * Get the deal.II DoFHandler associated with this particle world.
         *
         * @return The DoFHandler for this world.
         */
        const DoFHandler<dim> *get_dof_handler() const
        {
          return dof_handler;
        };

        /**
         * Set the deal.II DoFHandler associated with this particle world.
         *
         * @param [in] new_dh The new DoFHandler for this world.
         */
        void set_dof_handler(const DoFHandler<dim> *new_dh)
        {
          dof_handler = new_dh;
        };

        /**
         * Get the deal.II BlockVector solution associated with this particle
         * world.
         *
         * @return The BlockVector solution for this world.
         */
        const LinearAlgebra::BlockVector *get_solution() const
        {
          return solution;
        };

        /**
         * Set the deal.II BlockVector solution associated with this particle
         * world.
         *
         * @param [in] new_solution The new LinearAlgebra::BlockVector
         * solution for this world.
         */
        void set_solution(const LinearAlgebra::BlockVector *new_solution)
        {
          solution = new_solution;
        };

        /**
         * Set the particle Integrator scheme for this particle world.
         *
         * @param [in] new_integrator The new Integrator scheme for this
         * world.
         */
        void set_integrator(Integrator::Interface<dim, T> *new_integrator)
        {
          integrator = new_integrator;
        };

        /**
         * Set the MPI communicator for this world.
         *
         * @param [in] new_comm_world The new MPI_Comm object for this world.
         */
        void set_mpi_comm(const MPI_Comm new_comm_world)
        {
          communicator = new_comm_world;

          // Determine the size of the MPI comm world
          world_size = Utilities::MPI::n_mpi_processes(communicator);
          self_rank  = Utilities::MPI::this_mpi_process(communicator);
        };

        /**
         * Get the MPI communicator associated with this particle world.
         *
         * @return associated MPI_Comm
         */
        MPI_Comm mpi_comm()
        {
          return communicator;
        };

        /**
         * All processes must call this function when finished adding
         * particles to the world. This function will determine the total
         * number of particles.
         */
        void finished_adding_particles()
        {
          unsigned int local_num_particles = particles.size();

          MPI_Allreduce(&local_num_particles, &global_num_particles, 1, MPI_UNSIGNED, MPI_SUM, communicator);
        }

        /**
         * Add a particle to this world. If the specified cell does not exist
         * in the local subdomain an exception will be thrown.
         */
        void add_particle(const T &particle, const LevelInd &cell)
        {
          const typename parallel::distributed::Triangulation<dim>::active_cell_iterator it
          (triangulation, cell.first, cell.second);
          AssertThrow(it != triangulation->end(),
                      ExcMessage("Particles may only be added to cells in local subdomain."));
          particles.insert(std::make_pair(cell, particle));
        }

        /**
         * Access to particles in this world.
         */
        std::multimap<LevelInd, T> &get_particles()
        {
          return particles;
        };

        /**
         * Const access to particles in this world.
         */
        const std::multimap<LevelInd, T> &get_particles() const
        {
          return particles;
        };

        /**
         * Initialize the particle world by creating appropriate MPI data
         * types for transferring particles, and allocating memory for MPI
         * related functions.
         */
        void init()
        {
          int                                   *block_lens;
          MPI_Aint                              *indices;
          MPI_Datatype                          *old_types;
          std::vector<MPIDataInfo>              data_info;
          std::vector<MPIDataInfo>::iterator    it;
          int                                   num_entries, res, i;

          // Assert that all necessary parameters have been set
          AssertThrow (triangulation != NULL, ExcMessage ("Particle world triangulation must be set before calling init()."));
          AssertThrow (mapping != NULL, ExcMessage ("Particle world mapping must be set before calling init()."));
          AssertThrow (dof_handler != NULL, ExcMessage ("Particle world dof_handler must be set before calling init()."));
          AssertThrow (integrator != NULL, ExcMessage ("Particle world integrator must be set before calling init()."));

          // Construct MPI data type for this particle
          T::add_mpi_types(data_info);

          // And data associated with the integration scheme
          integrator->add_mpi_types(data_info);

          // Set up the block lengths, indices and internal types
          num_entries = data_info.size();
          block_lens = new int[num_entries];
          indices = new MPI_Aint[num_entries];
          old_types = new MPI_Datatype[num_entries];
          for (i=0; i<num_entries; ++i)
            {
              block_lens[i] = data_info[i].n_elements;
              indices[i] = (i == 0 ? 0 : indices[i-1]+sizeof(double)*data_info[i-1].n_elements);
              old_types[i] = MPI_DOUBLE;
            }

          // Create and commit the MPI type
          res = MPI_Type_struct(num_entries, block_lens, indices, old_types, &particle_type);
          if (res != MPI_SUCCESS) exit(-1);

          res = MPI_Type_commit(&particle_type);
          if (res != MPI_SUCCESS) exit(-1);

          // Delete temporary arrays
          delete [] old_types;
          delete [] indices;
          delete [] block_lens;

          // Initialize send/recv structures appropriately
          num_send = new int[world_size];
          num_recv = new int[world_size];
          send_offset = new int[world_size];
          recv_offset = new int[world_size];
          send_reqs = new MPI_Request[world_size];
          recv_reqs = new MPI_Request[world_size];
        };

        /**
         * Calculate the cells containing each particle for all particles.
         */
        void find_all_cells()
        {
          typename std::multimap<LevelInd, T>::iterator   it;
          std::multimap<LevelInd, T>                      tmp_map;
          LevelInd                                        found_cell;

          // Find the cells that the particles moved to
          tmp_map.clear();
          for (it=particles.begin(); it!=particles.end();)
            {
              found_cell = find_cell(it->second, it->first);
              if (found_cell != it->first)
                {
                  tmp_map.insert(std::make_pair(found_cell, it->second));
                  particles.erase(it++);
                }
              else ++it;
            }
          particles.insert(tmp_map.begin(),tmp_map.end());
        };

        /**
         * Advance particles by the specified timestep using the current
         * integration scheme.
         *
         * @param [in] timestep Length of timestep to integrate particle
         * movement
         * @param [in] solution Current Aspect solution vector
         */
        void advance_timestep(const double timestep, const LinearAlgebra::BlockVector &solution)
        {
          bool        continue_integrator = true;

          // Find the cells that the particles moved to
          find_all_cells();

          // If the triangulation changed, we may need to move particles between processors
          if (triangulation_changed) send_recv_particles();

          // If particles fell out of the mesh, put them back in at the closest point in the mesh
          move_particles_back_in_mesh();

          // Mark all particles to be checked for velocity based on position
          mark_particles_for_check();

          // Keep calling the integrator until it indicates it is finished
          while (continue_integrator)
            {
              // Starting out, particles must know which cells they belong to
              // Using this we can quickly interpolate the velocities
              get_particle_velocities(solution);

              // Call the integrator
              continue_integrator = integrator->integrate_step(this, timestep);

              // Find the cells that the particles moved to
              find_all_cells();

              // If particles fell out of the mesh, put them back in at the closest point in the mesh
              move_particles_back_in_mesh();

              // Swap particles between processors if needed
              send_recv_particles();
            }

          // Ensure we didn't lose any particles
          check_particle_count();
        };

        void move_particles_back_in_mesh()
        {
          // TODO: fix this to work with arbitrary meshes
        };

        /**
         * Mark all particles to be checked for velocity at their current
         * position
         */
        void mark_particles_for_check()
        {
          typename std::multimap<LevelInd, T>::iterator  it;
          for (it=particles.begin(); it!=particles.end(); ++it) it->second.set_vel_check(true);
        }

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
        LevelInd find_cell(T &particle, const LevelInd &cur_cell)
        {
          typename parallel::distributed::Triangulation<dim>::cell_iterator         it, found_cell;
          typename parallel::distributed::Triangulation<dim>::active_cell_iterator  ait;
          LevelInd    res;

          // First check the last recorded cell since particles will generally stay in the same area
          if (!triangulation_changed)
            {
              found_cell = typename parallel::distributed::Triangulation<dim>::cell_iterator(triangulation, cur_cell.first, cur_cell.second);
              if (found_cell != triangulation->end() && found_cell->point_inside(particle.get_location()) && found_cell->active())
                {
                  // If the cell is active, we're at the finest level of refinement and can finish
                  particle.set_local(found_cell->is_locally_owned());
                  return cur_cell;
                }
            }

          // Check all the cells on level 0 and recurse down
          for (it=triangulation->begin(0); it!=triangulation->end(0); ++it)
            {
              res = recursive_find_cell(particle, std::make_pair(it->level(), it->index()));
              if (res.first != -1 && res.second != -1) return res;
            }

          // If we couldn't find it there, we need to check the active cells
          // since the particle could be on a curved boundary not included in the
          // coarse grid
          for (ait=triangulation->begin_active(); ait!=triangulation->end(); ++ait)
            {
              if (ait->point_inside(particle.get_location()))
                {
                  particle.set_local(ait->is_locally_owned());
                  return std::make_pair(ait->level(), ait->index());
                }
            }

          // If it failed all these tests, the particle is outside the mesh
          particle.set_local(false);
          return std::make_pair(-1, -1);
        };

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
        void send_recv_particles()
        {
          typename std::multimap<LevelInd, T>::iterator  it;
          typename parallel::distributed::Triangulation<dim>::cell_iterator found_cell;
          int                 i;
          unsigned int        rank;
          std::vector<T>      send_particles;
          typename std::vector<T>::const_iterator    sit;

          // Go through the particles and take out those which need to be moved to another processor
          for (it=particles.begin(); it!=particles.end();)
            {
              if (!it->second.local())
                {
                  send_particles.push_back(it->second);
                  particles.erase(it++);
                }
              else
                {
                  ++it;
                }
            }

          // Determine the total number of particles we will send to other processors
          total_send = send_particles.size();
          for (rank=0; rank<world_size; ++rank)
            {
              if (rank != self_rank) num_send[rank] = total_send;
              else num_send[rank] = 0;
              send_offset[rank] = 0;
            }

          // Notify other processors how many particles we will be sending
          MPI_Alltoall(num_send, 1, MPI_INT, num_recv, 1, MPI_INT, communicator);

          total_recv = 0;
          for (rank=0; rank<world_size; ++rank)
            {
              recv_offset[rank] = total_recv;
              total_recv += num_recv[rank];
            }

          // Allocate space for sending and receiving particle data
          unsigned int        integrator_data_len = integrator->data_len();
          unsigned int        particle_data_len = T::data_len();
          std::vector<double> send_data, recv_data;

          // Set up the space for the received particle data
          recv_data.resize(total_recv*(integrator_data_len+particle_data_len));

          // Copy the particle data into the send array
          for (i=0,sit=send_particles.begin(); sit!=send_particles.end(); ++sit,++i)
            {
              sit->write_data(send_data);
              integrator->write_data(send_data, sit->get_id());
            }

          // Exchange the particle data between domains
          double *recv_data_ptr = &(recv_data[0]);
          double *send_data_ptr = &(send_data[0]);
          MPI_Alltoallv(send_data_ptr, num_send, send_offset, particle_type,
                        recv_data_ptr, num_recv, recv_offset, particle_type,
                        communicator);

          int put_in_domain = 0;
          unsigned int  pos = 0;
          // Put the received particles into the domain if they are in the triangulation
          for (i=0; i<total_recv; ++i)
            {
              T                   recv_particle;
              LevelInd            found_cell;
              pos = recv_particle.read_data(recv_data, pos);
              pos = integrator->read_data(recv_data, pos, recv_particle.get_id());
              found_cell = find_cell(recv_particle, std::make_pair(-1,-1));
              if (recv_particle.local())
                {
                  put_in_domain++;
                  particles.insert(std::make_pair(found_cell, recv_particle));
                }
            }
        };

        /**
         * Calculates the velocities for each particle at its location given
         * the input solution velocity field. The calculated velocities are
         * stored in the Particle objects for this world.
         *
         * @param [in] solution The current solution vector for this
         * simulation.
         */
        void get_particle_velocities(const LinearAlgebra::BlockVector &solution)
        {
          Vector<double>                single_res(dim+2);
          std::vector<Vector<double> >  result;
          Point<dim>                    velocity;
          unsigned int                  i, num_cell_particles;
          LevelInd                      cur_cell;
          typename std::multimap<LevelInd, T>::iterator  it, sit;
          typename DoFHandler<dim>::active_cell_iterator  found_cell;
          std::vector<Point<dim> >      particle_points;

          // Prepare the field function
          Functions::FEFieldFunction<dim, DoFHandler<dim>, LinearAlgebra::BlockVector> fe_value(*dof_handler, solution, *mapping);

          // Get the velocity for each cell at a time so we can take advantage of knowing the active cell
          for (it=particles.begin(); it!=particles.end();)
            {
              // Save a pointer to the first particle in this cell
              sit = it;

              // Get the current cell
              cur_cell = it->first;

              // Resize the vectors to the number of particles in this cell
              num_cell_particles = particles.count(cur_cell);
              particle_points.resize(num_cell_particles);

              // Get a vector of the particle locations in this cell
              i=0;
              while (it != particles.end() && it->first == cur_cell)
                {
                  if (it->second.vel_check()) particle_points[i++] = it->second.get_location();
                  it++;
                }
              result.resize(i, single_res);
              particle_points.resize(i);

              // Get the cell the particle is in
              found_cell = typename DoFHandler<dim>::active_cell_iterator(triangulation, cur_cell.first, cur_cell.second, dof_handler);

              // Interpolate the velocity field for each of the particles
              fe_value.set_active_cell(found_cell);
              fe_value.vector_value_list(particle_points, result);

              // Copy the resulting velocities to the appropriate vector
              it = sit;
              i = 0;
              while (it != particles.end() && it->first == cur_cell)
                {
                  if (it->second.vel_check())
                    {
                      for (int d=0; d<dim; ++d) velocity(d) = result[i](d);
                      it->second.set_velocity(velocity);
                      i++;
                    }
                  it++;
                }
            }
        };

        /**
         * Calculates the global sum of particles over all processes. This is
         * done to ensure no particles have fallen out of the simulation
         * domain.
         *
         * @return Total number of particles in simulation.
         */
        unsigned int get_global_particle_count()
        {
          return Utilities::MPI::sum (particles.size(), communicator);
        };

        /**
         * Checks that the number of particles in the simulation has not
         * unexpectedly changed. If the particle count changes then the
         * simulation will be aborted.
         */
        void check_particle_count()
        {
          unsigned int    global_particles = get_global_particle_count();

          AssertThrow (global_particles==global_num_particles,
                       ExcMessage ("Particle count unexpectedly changed."));
        };

        /**
         * Read or write the data of this object for serialization
         */
        template <class Archive>
        void serialize(Archive &ar, const unsigned int version)
        {
          ar &particles
          &global_num_particles
          ;
        }
    };
  }
}

#endif
