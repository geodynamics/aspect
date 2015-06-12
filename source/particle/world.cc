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

#include <aspect/particle/world.h>
#include <deal.II/numerics/fe_field_function.h>

namespace aspect
{
  namespace Particle
  {
    template <int dim>
    LevelInd
    World<dim>::recursive_find_cell(Particle<dim> &particle,
                                    const LevelInd cur_cell)
    {
      // If the particle is in the specified cell
      typename parallel::distributed::Triangulation<dim>::cell_iterator found_cell(&(this->get_triangulation()), cur_cell.first, cur_cell.second);
      if (found_cell != this->get_triangulation().end() && found_cell->point_inside(particle.get_location()))
        {
          // If the cell is active, we're at the finest level of refinement and can finish
          if (found_cell->active())
            {
              return cur_cell;
            }
          else
            {
              // Otherwise we need to search deeper
              for (unsigned int child_num=0; child_num<found_cell->n_children(); ++child_num)
                {
                  const typename parallel::distributed::Triangulation<dim>::cell_iterator child_cell = found_cell->child(child_num);
                  const LevelInd child_li(child_cell->level(), child_cell->index());
                  const LevelInd res = recursive_find_cell(particle, child_li);
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
      integrator = NULL;
    }

    template <int dim>
    World<dim>::~World()
    {}

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
    const Property::Manager<dim> &
    World<dim>::get_manager() const
    {
      return *property_manager;
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
    World<dim>::add_particle(const Particle<dim> &particle, const LevelInd &cell)
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
      const unsigned int solution_components = this->introspection().n_components;
      Vector<double> value(solution_components);
      std::vector<Tensor<1,dim> > gradient (solution_components,Tensor<1,dim>());
      std::vector<Vector<double> >  values(50,value);
      std::vector<std::vector<Tensor<1,dim> > > gradients(50,gradient);
      std::vector<Point<dim> >      particle_points(50);

      const DoFHandler<dim> *dof_handler = &(this->get_dof_handler());
      const typename parallel::distributed::Triangulation<dim> *triangulation = &(this->get_triangulation());

      // Prepare the field function
      Functions::FEFieldFunction<dim, DoFHandler<dim>, LinearAlgebra::BlockVector> fe_value(*dof_handler, this->get_solution(), this->get_mapping());

      // Get the velocity for each cell at a time so we can take advantage of knowing the active cell
      typename std::multimap<LevelInd, Particle<dim> >::iterator  it, sit;
      for (it=particles.begin(); it!=particles.end();)
        {
          // Save a pointer to the first particle in this cell
          sit = it;

          // Get the current cell
          const LevelInd cur_cell = it->first;

          // Resize the vectors to the number of particles in this cell
          particle_points.resize(particles.count(cur_cell));

          // Get a vector of the particle locations in this cell
          unsigned int i=0;
          while (it != particles.end() && it->first == cur_cell)
            {
              particle_points[i++] = it->second.get_location();
              it++;
            }
          values.resize(i, value);
          gradients.resize(i,gradient);
          particle_points.resize(i);

          // Get the cell the particle is in
          const typename DoFHandler<dim>::active_cell_iterator found_cell (triangulation, cur_cell.first, cur_cell.second, dof_handler);

          // Interpolate the velocity field for each of the particles
          fe_value.set_active_cell(found_cell);
          fe_value.vector_value_list(particle_points, values);
          fe_value.vector_gradient_list(particle_points, gradients);

          // Copy the resulting velocities to the appropriate vector
          it = sit;
          i = 0;
          while (it != particles.end() && it->first == cur_cell)
            {
              property_manager->initialize_particle(it->second,
                                                    values[i],
                                                    gradients[i]);
              i++;
              it++;
            }
        }
    }

    template <int dim>
    void
    World<dim>::update_particles()
    {
      const unsigned int solution_components = this->introspection().n_components;
      Vector<double> value(solution_components);
      std::vector<Tensor<1,dim> > gradient (solution_components,Tensor<1,dim>());
      std::vector<Vector<double> >  values(50,value);
      std::vector<std::vector<Tensor<1,dim> > > gradients(50,gradient);
      std::vector<Point<dim> >      particle_points(50);

      const DoFHandler<dim> *dof_handler = &(this->get_dof_handler());
      const typename parallel::distributed::Triangulation<dim> *triangulation = &(this->get_triangulation());

      // Prepare the field function
      Functions::FEFieldFunction<dim, DoFHandler<dim>, LinearAlgebra::BlockVector> fe_value(*dof_handler, this->get_solution(), this->get_mapping());

      // Get the velocity for each cell at a time so we can take advantage of knowing the active cell
      typename std::multimap<LevelInd, Particle<dim> >::iterator  sit;
      for (typename std::multimap<LevelInd, Particle<dim> >::iterator
           it=particles.begin(); it!=particles.end();)
        {
          // Save a pointer to the first particle in this cell
          sit = it;

          // Get the current cell
          const LevelInd cur_cell = it->first;

          // Resize the vectors to the number of particles in this cell
          particle_points.resize(particles.count(cur_cell));

          // Get a vector of the particle locations in this cell
          unsigned int i=0;
          while (it != particles.end() && it->first == cur_cell)
            {
              particle_points[i++] = it->second.get_location();
              it++;
            }
          values.resize(i, value);
          gradients.resize(i,gradient);
          particle_points.resize(i);

          // Get the cell the particle is in
          const typename DoFHandler<dim>::active_cell_iterator found_cell =
            typename DoFHandler<dim>::active_cell_iterator(triangulation, cur_cell.first, cur_cell.second, dof_handler);

          // Interpolate the velocity field for each of the particles
          fe_value.set_active_cell(found_cell);
          fe_value.vector_value_list(particle_points, values);
          fe_value.vector_gradient_list(particle_points, gradients);

          // Copy the resulting velocities to the appropriate vector
          it = sit;
          i = 0;
          while (it != particles.end() && it->first == cur_cell)
            {
              property_manager->update_particle(it->second,
                                                values[i],
                                                gradients[i]);
              i++;
              it++;
            }
        }
    }

    template <int dim>
    std::multimap<LevelInd, Particle<dim> > &
    World<dim>::get_particles()
    {
      return particles;
    }

    template <int dim>
    const std::multimap<LevelInd, Particle<dim> > &
    World<dim>::get_particles() const
    {
      return particles;
    }

    template <int dim>
    void
    World<dim>::get_data_info(std::vector<std::string> &names,
                              std::vector<unsigned int> &length) const
    {
      property_manager->get_data_info(names,length);
    }

    template <int dim>
    void
    World<dim>::init()
    {
      this->get_triangulation().signals.post_refinement.connect(std_cxx11::bind(&World::mesh_changed, std_cxx1x::ref(*this)));
    }

    template <int dim>
    void
    World<dim>::find_all_cells(std::vector<Particle<dim> > &lost_particles)
    {
      std::multimap<LevelInd, Particle<dim> > moved_particles;

      // Find the cells that the particles moved to.
      // Note that the iterator in the following loop is increased in a
      // very particular way, because it is changed, if elements
      // get erased. A change can result in invalid memory access.
      moved_particles.clear();
      typename std::multimap<LevelInd, Particle<dim> >::iterator   it;
      for (it=particles.begin(); it!=particles.end();)
        {
          const LevelInd found_cell = find_cell(it->second, it->first);
          if (found_cell != it->first)
            {
              // Reinsert the particle into our domain if we found its cell
              // Mark it for MPI transfer otherwise
              typename parallel::distributed::Triangulation<dim>::cell_iterator cell(&this->get_triangulation(), found_cell.first, found_cell.second);
              if (cell->is_locally_owned())
                moved_particles.insert(std::make_pair(found_cell, it->second));
              else
                lost_particles.push_back(it->second);

              particles.erase(it++);
            }
          else
            ++it;
        }
      particles.insert(moved_particles.begin(),moved_particles.end());
    }

    template <int dim>
    void
    World<dim>::advance_timestep()
    {
      bool        continue_integrator = true;

      // If the triangulation changed, we may need to move particles between processors
      if (triangulation_changed)
        {
          // Find the cells that the particles moved to
          std::vector<Particle<dim> > lost_particles;
          find_all_cells(lost_particles);
          send_recv_particles(lost_particles);
        }

      // If particles fell out of the mesh, put them back in at the closest point in the mesh
      move_particles_back_in_mesh();

      // Keep calling the integrator until it indicates it is finished
      while (continue_integrator)
        {
          // Starting out, particles must know which cells they belong to
          // Using this we can quickly interpolate the velocities
          std::vector<Tensor<1,dim> > old_velocities(particles.size());
          std::vector<Tensor<1,dim> > velocities(particles.size());

          get_particle_velocities(old_velocities,velocities);

          // Call the integrator
          continue_integrator = integrator->integrate_step(particles,
                                                           old_velocities,
                                                           velocities,
                                                           this->get_old_timestep());

          std::vector<Particle<dim> > lost_particles;
          // Find the cells that the particles moved to
          find_all_cells(lost_particles);

          // If particles fell out of the mesh, put them back in at the closest point in the mesh
          move_particles_back_in_mesh();

          // Swap particles between processors if needed
          send_recv_particles(lost_particles);
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
    LevelInd
    World<dim>::find_cell(Particle<dim> &particle, const LevelInd &cur_cell)
    {
      const typename parallel::distributed::Triangulation<dim> *triangulation = &(this->get_triangulation());

      // First check the last recorded cell since particles will generally stay in the same area
      if (!triangulation_changed)
        {
          typename parallel::distributed::Triangulation<dim>::cell_iterator found_cell(triangulation, cur_cell.first, cur_cell.second);
          if (found_cell != triangulation->end() && found_cell->point_inside(particle.get_location()) && found_cell->active())
            {
              // If the cell is active, we're at the finest level of refinement and can finish
              return cur_cell;
            }
        }

      // Check all the cells on level 0 and recurse down
      for (typename parallel::distributed::Triangulation<dim>::cell_iterator
           it=triangulation->begin(0); it!=triangulation->end(0); ++it)
        {
          const LevelInd res = recursive_find_cell(particle, std::make_pair(it->level(), it->index()));
          if (res.first != -1 && res.second != -1) return res;
        }

      // If we couldn't find it there, we need to check the active cells
      // since the particle could be on a curved boundary not included in the
      // coarse grid
      for (typename parallel::distributed::Triangulation<dim>::active_cell_iterator
           ait=triangulation->begin_active(); ait!=triangulation->end(); ++ait)
        {
          if (ait->point_inside(particle.get_location()))
            {
              return std::make_pair(ait->level(), ait->index());
            }
        }

      // If it failed all these tests, the particle is outside the mesh
      return std::make_pair(-1, -1);
    }

    template <int dim>
    void
    World<dim>::send_recv_particles(const std::vector<Particle <dim> > &send_particles)
    {
      // Determine the size of the MPI comm world
      const unsigned int world_size = Utilities::MPI::n_mpi_processes(this->get_mpi_communicator());
      const unsigned int self_rank  = Utilities::MPI::this_mpi_process(this->get_mpi_communicator());

      // Determine the total number of particles we will send to other processors
      int total_send_particles = send_particles.size();
      const int total_recv_particles = Utilities::MPI::sum (total_send_particles, this->get_mpi_communicator());

      // Allocate space for sending and receiving particle data
      std::vector<double> send_data(total_send_particles * (property_manager->get_data_len() + integrator->data_length()));
      std::vector<double>::iterator data = send_data.begin();

      // Copy the particle data into the send array
      for (typename std::vector<Particle<dim> >::const_iterator particle = send_particles.begin(); particle!=send_particles.end(); ++particle)
        {
          particle->write_data(data);
          integrator->write_data(data, particle->get_id());
        }

      int num_send_data = send_data.size();

      AssertThrow(data == send_data.end(),
                  ExcMessage("The amount of data written into the array that is send to other processes "
                             "is inconsistent with the number and size of particles."));

      // Notify other processors how many particles we will be sending
      std::vector<int> num_recv_data(world_size,0);
      std::vector<int> recv_offset(world_size,0);

      MPI_Allgather(&num_send_data, 1, MPI_INT, &(num_recv_data[0]), 1, MPI_INT, this->get_mpi_communicator());

      int total_recv_data = 0;
      for (unsigned int rank=0; rank<world_size; ++rank)
        {
          recv_offset[rank] = total_recv_data;
          total_recv_data += num_recv_data[rank];
        }

      // Set up the space for the received particle data
      std::vector<double> recv_data(total_recv_data);

      // Exchange the particle data between domains
      MPI_Allgatherv (&(send_data[0]), num_send_data, MPI_DOUBLE,
                      &(recv_data[0]), &(num_recv_data[0]), &(recv_offset[0]), MPI_DOUBLE,
                      this->get_mpi_communicator());

      // Put the received particles into the domain if they are in the triangulation
      std::vector<double>::const_iterator recv_data_it = recv_data.begin();

      for (int i=0; i<total_recv_particles; ++i)
        {
          Particle<dim> recv_particle(recv_data_it,property_manager->get_data_len());
          integrator->read_data(recv_data_it, recv_particle.get_id());

          const LevelInd found_cell = find_cell(recv_particle, std::make_pair(-1,-1));
          const typename parallel::distributed::Triangulation<dim>::cell_iterator cell(&this->get_triangulation(), found_cell.first, found_cell.second);

          if (cell->is_locally_owned())
            {
              particles.insert(std::make_pair(found_cell, recv_particle));
            }
        }

      AssertThrow(recv_data_it == recv_data.end(),
                  ExcMessage("The amount of data that was read into new particles "
                             "does not match the amount of data sent around."));
    }

    template <int dim>
    void
    World<dim>::get_particle_velocities(std::vector<Tensor<1,dim> > &velocities,
                                        std::vector<Tensor<1,dim> > &old_velocities)
    {
      Vector<double>                single_res(dim);
      std::vector<Vector<double> >  result(50,single_res);
      std::vector<Vector<double> >  old_result(50,single_res);
      std::vector<Point<dim> >      particle_points(50);

      const DoFHandler<dim> *dof_handler = &(this->get_dof_handler());
      const typename parallel::distributed::Triangulation<dim> *triangulation = &(this->get_triangulation());

      // Prepare the field function
      Functions::FEFieldFunction<dim, DoFHandler<dim>, LinearAlgebra::BlockVector> fe_value(*dof_handler, this->get_solution(), this->get_mapping());
      Functions::FEFieldFunction<dim, DoFHandler<dim>, LinearAlgebra::BlockVector> old_fe_value(*dof_handler, this->get_old_solution(), this->get_mapping());

      unsigned int particle_idx = 0;
      // Get the velocity for each cell at a time so we can take advantage of knowing the active cell
      for (typename std::multimap<LevelInd, Particle<dim> >::iterator
           it=particles.begin(); it!=particles.end();)
        {
          // Get the current cell
          const LevelInd cur_cell = it->first;

          // Resize the vectors to the number of particles in this cell
          const unsigned int num_cell_particles = particles.count(cur_cell);
          particle_points.resize(num_cell_particles);

          // Get a vector of the particle locations in this cell
          unsigned int n_particles_in_cell = 0;
          while (it != particles.end() && it->first == cur_cell)
            {
              particle_points[n_particles_in_cell++] = it->second.get_location();
              it++;
            }

          result.resize(n_particles_in_cell, single_res);
          old_result.resize(n_particles_in_cell, single_res);
          particle_points.resize(n_particles_in_cell);

          // Get the cell the particle is in
          const typename DoFHandler<dim>::active_cell_iterator found_cell (triangulation, cur_cell.first, cur_cell.second, dof_handler);

          // Interpolate the velocity field for each of the particles
          old_fe_value.set_active_cell(found_cell);
          old_fe_value.vector_value_list(particle_points, old_result);

          fe_value.set_active_cell(found_cell);
          fe_value.vector_value_list(particle_points, result);

          // Copy the resulting velocities to the appropriate vector#
          for (unsigned int id = 0; id < n_particles_in_cell; ++id)
            {
              Tensor<1,dim> velocity, old_velocity;
              for (unsigned int d=0; d<dim; ++d)
                {
                  velocity[d] = result[id][d];
                  old_velocity[d] = old_result[id][d];
                }
              velocities[particle_idx] = velocity;
              old_velocities[particle_idx] = old_velocity;
              particle_idx++;
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
