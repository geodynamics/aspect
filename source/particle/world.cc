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

#include <aspect/particle/world.h>
#include <deal.II/numerics/fe_field_function.h>

namespace aspect
{
  namespace Particle
  {
    template <int dim>
    LevelInd
    World<dim>::recursive_find_cell(BaseParticle<dim> &particle,
                                    const LevelInd cur_cell)
        {
      typename parallel::distributed::Triangulation<dim>::cell_iterator it, found_cell, child_cell;
      unsigned int    child_num;
      LevelInd        res, child_li;

      // If the particle is in the specified cell
      found_cell = typename parallel::distributed::Triangulation<dim>::cell_iterator(&(this->get_triangulation()), cur_cell.first, cur_cell.second);
      if (found_cell != this->get_triangulation().end() && found_cell->point_inside(particle.get_location()))
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
        }

        /**
         * Called by listener functions to indicate that the mesh of this
         * subdomain has changed.
         */
    template <int dim>
    void
    World<dim>::mesh_changed()
    {
      triangulation_changed = true;
    }

    template <int dim>
    World<dim>::World()
    {
      triangulation_changed = true;
      total_send = total_recv = 0;
      world_size = self_rank = 0;
      num_send = num_recv = send_offset = recv_offset = NULL;
      send_reqs = recv_reqs = NULL;
      integrator = NULL;
    }

    template <int dim>
    World<dim>::~World()
    {
      if (world_size) MPI_Type_free(&particle_type);

      if (num_send) delete [] num_send;
      if (num_recv) delete [] num_recv;
      if (send_offset) delete [] send_offset;
      if (recv_offset) delete [] recv_offset;
      if (send_reqs) delete [] send_reqs;
      if (recv_reqs) delete [] recv_reqs;
    }

    template <int dim>
    void
    World<dim>::set_integrator(Integrator::Interface<dim> *new_integrator)
    {
      integrator = new_integrator;
    }

    template <int dim>
    void
    World<dim>::set_manager(Property::Manager<dim> *new_manager)
    {
      property_manager = new_manager;
    }

    template <int dim>
    void
    World<dim>::finished_adding_particles()
    {
      unsigned int local_num_particles = particles.size();

      MPI_Allreduce(&local_num_particles, &global_num_particles, 1, MPI_UNSIGNED, MPI_SUM, this->get_mpi_communicator());
    }

    template <int dim>
    void
    World<dim>::add_particle(const BaseParticle<dim> &particle, const LevelInd &cell)
    {
      const typename parallel::distributed::Triangulation<dim>::active_cell_iterator it
      (&(this->get_triangulation()), cell.first, cell.second);
      AssertThrow(it != this->get_triangulation().end(),
                  ExcMessage("Particles may only be added to cells in local subdomain."));
      particles.insert(std::make_pair(cell, particle));
    }

    template <int dim>
    void
    World<dim>::initialize_particles()
    {
      Vector<double>                single_res(dim+2+this->n_compositional_fields());
      std::vector<Vector<double> >  result(50,single_res);
      unsigned int                  i, num_cell_particles;
      LevelInd                      cur_cell;
      typename std::multimap<LevelInd, BaseParticle<dim> >::iterator  it, sit;
      typename DoFHandler<dim>::active_cell_iterator  found_cell;
      std::vector<Point<dim> >      particle_points(50);

      const DoFHandler<dim> *dof_handler = &(this->get_dof_handler());
      const Mapping<dim> *mapping = &(this->get_mapping());
      const LinearAlgebra::BlockVector *solution = &(this->get_solution());
      const typename parallel::distributed::Triangulation<dim> *triangulation = &(this->get_triangulation());

      // Prepare the field function
      Functions::FEFieldFunction<dim, DoFHandler<dim>, LinearAlgebra::BlockVector> fe_value(*dof_handler, *solution, *mapping);

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
              particle_points[i++] = it->second.get_location();
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
              property_manager->initialize_particle(it->second,
                  result[i]);
              i++;
              it++;
            }
        }
    }

    template <int dim>
    void
    World<dim>::update_particles()
    {
      Vector<double>                single_res(dim+2+this->n_compositional_fields());
      std::vector<Vector<double> >  result(50,single_res);
      unsigned int                  i, num_cell_particles;
      LevelInd                      cur_cell;
      typename std::multimap<LevelInd, BaseParticle<dim> >::iterator  it, sit;
      typename DoFHandler<dim>::active_cell_iterator  found_cell;
      std::vector<Point<dim> >      particle_points(50);

      const DoFHandler<dim> *dof_handler = &(this->get_dof_handler());
      const Mapping<dim> *mapping = &(this->get_mapping());
      const LinearAlgebra::BlockVector *solution = &(this->get_solution());
      const typename parallel::distributed::Triangulation<dim> *triangulation = &(this->get_triangulation());

      // Prepare the field function
      Functions::FEFieldFunction<dim, DoFHandler<dim>, LinearAlgebra::BlockVector> fe_value(*dof_handler, *solution, *mapping);

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
              particle_points[i++] = it->second.get_location();
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
              property_manager->update_particle(it->second,
                  result[i]);
              i++;
              it++;
            }
        }
    }

    template <int dim>
    std::multimap<LevelInd, BaseParticle<dim> > &
    World<dim>::get_particles()
    {
      return particles;
    }

    template <int dim>
    const std::multimap<LevelInd, BaseParticle<dim> > &
    World<dim>::get_particles() const
    {
      return particles;
    }

    template <int dim>
    const std::vector<MPIDataInfo> &
    World<dim>::get_mpi_datainfo() const
    {
      return data_info;
    }

    template <int dim>
    void
    World<dim>::init()
    {
      int                                   *block_lens;
      MPI_Aint                              *indices;
      MPI_Datatype                          *old_types;
      std::vector<MPIDataInfo>::iterator    it;
      int                                   num_entries, res, i;

      // Assert that all necessary parameters have been set
      AssertThrow (integrator != NULL, ExcMessage ("Particle world integrator must be set before calling init()."));
      AssertThrow (property_manager != NULL, ExcMessage ("Particle world integrator must be set before calling init()."));

      const typename std::multimap<LevelInd, BaseParticle<dim> >::iterator   pit (particles.begin());

      // Construct MPI data type for this particle
      property_manager->add_mpi_types(data_info);

      // And data associated with the integration scheme
      integrator->add_mpi_types(data_info);

      this->get_triangulation().signals.post_refinement.connect(std_cxx11::bind(&World::mesh_changed, std_cxx1x::ref(*this)));

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

      // Determine the size of the MPI comm world
      world_size = Utilities::MPI::n_mpi_processes(this->get_mpi_communicator());
      self_rank  = Utilities::MPI::this_mpi_process(this->get_mpi_communicator());

      // Initialize send/recv structures appropriately
      num_send = new int[world_size];
      num_recv = new int[world_size];
      send_offset = new int[world_size];
      recv_offset = new int[world_size];
      send_reqs = new MPI_Request[world_size];
      recv_reqs = new MPI_Request[world_size];
    }

    template <int dim>
    void
    World<dim>::find_all_cells()
    {
      typename std::multimap<LevelInd, BaseParticle<dim> >::iterator   it;
      std::multimap<LevelInd, BaseParticle<dim> >                      tmp_map;
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
    }

    template <int dim>
    void
    World<dim>::advance_timestep(const double timestep, const LinearAlgebra::BlockVector &solution)
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
          continue_integrator = integrator->integrate_step(particles, timestep);

          // Find the cells that the particles moved to
          find_all_cells();

          // If particles fell out of the mesh, put them back in at the closest point in the mesh
          move_particles_back_in_mesh();

          // Swap particles between processors if needed
          send_recv_particles();
        }

      // Update particle properties
      if (property_manager->need_update())
        update_particles();

      // Ensure we didn't lose any particles
      check_particle_count();
    }

    template <int dim>
    void
    World<dim>::move_particles_back_in_mesh()
    {
      // TODO: fix this to work with arbitrary meshes
    }

    template <int dim>
    void
    World<dim>::mark_particles_for_check()
    {
      typename std::multimap<LevelInd, BaseParticle<dim> >::iterator  it;
      for (it=particles.begin(); it!=particles.end(); ++it) it->second.set_vel_check(true);
    }

    template <int dim>
    LevelInd
    World<dim>::find_cell(BaseParticle<dim> &particle, const LevelInd &cur_cell)
    {
      typename parallel::distributed::Triangulation<dim>::cell_iterator         it, found_cell;
      typename parallel::distributed::Triangulation<dim>::active_cell_iterator  ait;
      LevelInd    res;

      const typename parallel::distributed::Triangulation<dim> *triangulation = &(this->get_triangulation());

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
    }

    template <int dim>
    void
    World<dim>::send_recv_particles()
    {
      typename std::multimap<LevelInd, BaseParticle<dim> >::iterator  it;
      typename parallel::distributed::Triangulation<dim>::cell_iterator found_cell;
      int                 i;
      unsigned int        rank;
      std::vector<BaseParticle<dim> >      send_particles;
      typename std::vector<BaseParticle<dim> >::const_iterator    sit;

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
      MPI_Alltoall(num_send, 1, MPI_INT, num_recv, 1, MPI_INT, this->get_mpi_communicator());

      total_recv = 0;
      for (rank=0; rank<world_size; ++rank)
        {
          recv_offset[rank] = total_recv;
          total_recv += num_recv[rank];
        }

      // Allocate space for sending and receiving particle data
      unsigned int        integrator_data_len = integrator->data_len();
      unsigned int        particle_data_len = property_manager->get_data_len();
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
                    this->get_mpi_communicator());

      int put_in_domain = 0;
      unsigned int  pos = 0;
      // Put the received particles into the domain if they are in the triangulation
      for (i=0; i<total_recv; ++i)
        {
          BaseParticle<dim>       recv_particle;
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
    }

    template <int dim>
    void
    World<dim>::get_particle_velocities(const LinearAlgebra::BlockVector &solution)
    {
      Vector<double>                single_res(dim);
      std::vector<Vector<double> >  result(50,single_res);
      Point<dim>                    velocity;
      unsigned int                  i, num_cell_particles;
      LevelInd                      cur_cell;
      typename std::multimap<LevelInd, BaseParticle<dim> >::iterator  it, sit;
      typename DoFHandler<dim>::active_cell_iterator  found_cell;
      std::vector<Point<dim> >      particle_points(50);

      const DoFHandler<dim> *dof_handler = &(this->get_dof_handler());
      const Mapping<dim> *mapping = &(this->get_mapping());
      const typename parallel::distributed::Triangulation<dim> *triangulation = &(this->get_triangulation());

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
    }

    template <int dim>
    unsigned int
    World<dim>::get_global_particle_count() const
    {
      return Utilities::MPI::sum (particles.size(), this->get_mpi_communicator());
    }

    template <int dim>
    void
    World<dim>::check_particle_count()
    {
      const unsigned int global_particles = get_global_particle_count();

      global_num_particles = global_particles;
    }

    template <int dim>
    template <class Archive>
    void
    World<dim>::serialize(Archive &ar, const unsigned int version)
    {
        ar &particles
        &global_num_particles
        ;
    }
  }
}



// explicit instantiation of the functions we implement in this file
namespace aspect
{
  namespace Particle
  {
#define INSTANTIATE(dim) \
    template class World<dim>;

    ASPECT_INSTANTIATE(INSTANTIATE)
  }
}
