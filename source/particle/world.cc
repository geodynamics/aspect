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
#include <aspect/global.h>

#include <deal.II/numerics/fe_field_function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/grid_tools.h>
#include <boost/serialization/map.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

namespace aspect
{
  namespace Particle
  {
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
    World<dim>::init()
    {
      this->get_triangulation().signals.post_refinement.connect(std_cxx11::bind(&World<dim>::mesh_changed, std_cxx1x::ref(*this)));
    }

    template <int dim>
    unsigned int
    World<dim>::get_max_tracer_per_cell() const
    {
      unsigned int local_max_tracer_per_cell(0);

      for (typename parallel::distributed::Triangulation<dim>::active_cell_iterator cell=this->get_triangulation().begin_active();
           cell!=this->get_triangulation().end(); ++cell)
        {
          if (cell->is_locally_owned())
            {
              const LevelInd found_cell = std::make_pair<int, int> (cell->level(),cell->index());
              const std::pair<typename std::multimap<LevelInd, Particle<dim> >::const_iterator, typename std::multimap<LevelInd, Particle<dim> >::const_iterator> particles_in_cell
              = particles.equal_range(found_cell);
              const unsigned int tracers_in_cell = std::distance(particles_in_cell.first,particles_in_cell.second);
              local_max_tracer_per_cell = std::max(local_max_tracer_per_cell,
                                                   tracers_in_cell);
            }
        }

      //std::cout << "Tracers per cell: " << local_max_tracer_per_cell << std::endl;
      const unsigned int global_max_tracer_per_cell = Utilities::MPI::max(local_max_tracer_per_cell,this->get_mpi_communicator());
      const size_t size = std::pow(2,dim) * global_max_tracer_per_cell;
      //std::cout << "Size per cell: " << size << std::endl;
      return size;
    }

    template <int dim>
    void
    World<dim>::store_tracers(const typename parallel::distributed::Triangulation<dim>::cell_iterator &cell,
                              const typename parallel::distributed::Triangulation<dim>::CellStatus status,
                              void *data)
    {
      unsigned int n_particles_in_cell = 0;

      if (status == parallel::distributed::Triangulation<dim>::CellStatus::CELL_PERSIST
          || status == parallel::distributed::Triangulation<dim>::CellStatus::CELL_REFINE)
        {
          const LevelInd found_cell = std::make_pair<int, int> (cell->level(),cell->index());
          const std::pair<typename std::multimap<LevelInd, Particle<dim> >::iterator, typename std::multimap<LevelInd, Particle<dim> >::iterator> particles_in_cell
          = particles.equal_range(found_cell);
          n_particles_in_cell = std::distance(particles_in_cell.first,particles_in_cell.second);

          unsigned int* ndata = static_cast<unsigned int*> (data);
          *ndata++ = n_particles_in_cell;
          data = static_cast<void*> (ndata);

          for (typename std::multimap<LevelInd, Particle<dim> >::iterator particle = particles_in_cell.first;
              particle != particles_in_cell.second;++particle)
            {
              particle->second.write_data(data);
            }
          particles.erase(particles_in_cell.first,particles_in_cell.second);
        }
      else if (status == parallel::distributed::Triangulation<dim>::CellStatus::CELL_COARSEN)
        {
          for (unsigned int child_index = 0; child_index < cell->number_of_children(); ++child_index)
            {
              const typename parallel::distributed::Triangulation<dim>::cell_iterator child = cell->child(child_index);
              const LevelInd found_cell = std::make_pair<int, int> (child->level(),child->index());
              const std::pair<typename std::multimap<LevelInd, Particle<dim> >::iterator, typename std::multimap<LevelInd, Particle<dim> >::iterator> particles_in_cell
              = particles.equal_range(found_cell);
              n_particles_in_cell += std::distance(particles_in_cell.first,particles_in_cell.second);
            }

          unsigned int* ndata = static_cast<unsigned int*> (data);
          *ndata++ = n_particles_in_cell;

          data = static_cast<void*> (ndata);

          for (unsigned int child_index = 0; child_index < cell->number_of_children(); ++child_index)
            {
              const typename parallel::distributed::Triangulation<dim>::cell_iterator child = cell->child(child_index);
              const LevelInd found_cell = std::make_pair<int, int> (child->level(),child->index());
              const std::pair<typename std::multimap<LevelInd, Particle<dim> >::iterator, typename std::multimap<LevelInd, Particle<dim> >::iterator>
              particles_in_cell = particles.equal_range(found_cell);

              for (typename std::multimap<LevelInd, Particle<dim> >::iterator particle = particles_in_cell.first;
                  particle != particles_in_cell.second;++particle)
                {
                  particle->second.write_data(data);
                }
              particles.erase(particles_in_cell.first,particles_in_cell.second);
            }
        }

      //std::cout << "Store Cell index: " << cell->index() << ". Particles: " << n_particles_in_cell << std::endl;
    }

    template <int dim>
    void
    World<dim>::load_tracers(const typename parallel::distributed::Triangulation<dim>::cell_iterator &cell,
                             const typename parallel::distributed::Triangulation<dim>::CellStatus status,
                             const void *data)
    {
      const unsigned int* n_particles_in_cell = static_cast<const unsigned int*> (data);
      const unsigned int particles_in_cell = *n_particles_in_cell++;
      //std::cout << "Loading Cell index: " << cell->index() << ". Particles: " << particles_in_cell <<  std::endl;

      void *pdata = (void *) n_particles_in_cell;

      for (unsigned int i = 0; i < particles_in_cell;++i)
        {
          Particle<dim> p(pdata,property_manager->get_data_len());

          if (status == parallel::distributed::Triangulation<dim>::CellStatus::CELL_COARSEN
              || status == parallel::distributed::Triangulation<dim>::CellStatus::CELL_PERSIST)
            particles.insert(std::make_pair(std::make_pair(cell->level(),cell->index()),p));
          else if (status == parallel::distributed::Triangulation<dim>::CellStatus::CELL_REFINE)
            {
              for (unsigned int child_index = 0; child_index < cell->number_of_children(); ++child_index)
                {
                  const typename parallel::distributed::Triangulation<dim>::cell_iterator child = cell->child(child_index);

                  try
                  {
                      const Point<dim> p_unit = this->get_mapping().transform_real_to_unit_cell(child, p.get_location());
                      if (GeometryInfo<dim>::is_inside_unit_cell(p_unit))
                        {
                          particles.insert(std::make_pair(std::make_pair(child->level(),child->index()),p));
                        }
                  }
                  catch (...)
                  {}
                }
            }
        }
    }

    template <int dim>
    void
    World<dim>::mesh_changed()
    {
      triangulation_changed = true;
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
    const Property::Manager<dim> &
    World<dim>::get_manager() const
    {
      return *property_manager;
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
    World<dim>::find_all_cells()
    {
      std::multimap<types::subdomain_id,Particle<dim> > lost_particles;
      std::multimap<LevelInd, Particle<dim> > moved_particles;

      // Find the cells that the particles moved to.
      // Note that the iterator in the following loop is increased in a
      // very particular way, because it is changed, if elements
      // get erased. A change can result in invalid memory access.
      moved_particles.clear();
      typename std::multimap<LevelInd, Particle<dim> >::iterator   it;
      for (it=particles.begin(); it!=particles.end();)
        {
          if ((it->first != std::make_pair(-1,-1)) && !triangulation_changed)
            {
              const typename parallel::distributed::Triangulation<dim>::active_cell_iterator
              old_cell (&(this->get_triangulation()), it->first.first, it->first.second);

              try
                {
                  const Point<dim> p_unit = this->get_mapping().transform_real_to_unit_cell(old_cell, it->second.get_location());
                  if (GeometryInfo<dim>::is_inside_unit_cell(p_unit))
                    {
                      ++it;
                      continue;
                    }
                }
              catch (...)
                {}
            }

          typename parallel::distributed::Triangulation<dim>::active_cell_iterator cell;
          try
            {
              cell = (GridTools::find_active_cell_around_point<> (this->get_mapping(), this->get_triangulation(), it->second.get_location())).first;
            }
          catch (GridTools::ExcPointNotFound<dim> &)
            {
              // If we can find no cell for this particle it has left the domain due
              // to an integration error or an open boundary. Simply remove the
              // tracer in this case.
              particles.erase(it++);
              continue;
            }


          // Reinsert the particle into our domain if we found its cell
          // Mark it for MPI transfer otherwise
          if (cell->is_locally_owned())
            {
              const LevelInd found_cell = std::make_pair(cell->level(),cell->index());
              moved_particles.insert(std::make_pair(found_cell, it->second));
            }
          else
            lost_particles.insert(std::make_pair(cell->subdomain_id(),it->second));

          particles.erase(it++);
        }
      particles.insert(moved_particles.begin(),moved_particles.end());


      // If particles fell out of the mesh, put them back in at the closest point in the mesh
      move_particles_back_in_mesh();

      // Swap lost particles between processors if any process lost some
      const unsigned int global_lost_particles = Utilities::MPI::sum(lost_particles.size(),
                                                                     this->get_mpi_communicator());
      if (global_lost_particles > 0)
        send_recv_particles(lost_particles);
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
          //find_all_cells();
          triangulation_changed = false;
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

          // Find the cells that the particles moved to
          find_all_cells();
        }

      // Update particle properties
      if (property_manager->need_update() == Property::update_time_step)
        update_particles();
    }

    template <int dim>
    void
    World<dim>::move_particles_back_in_mesh()
    {
      // TODO: fix this to work with arbitrary meshes
    }

    template <int dim>
    void
    World<dim>::send_recv_particles(const std::multimap<types::subdomain_id,Particle <dim> > &send_particles)
    {
      // Determine the size of the MPI comm world
      const unsigned int world_size = Utilities::MPI::n_mpi_processes(this->get_mpi_communicator());
      const unsigned int self_rank  = Utilities::MPI::this_mpi_process(this->get_mpi_communicator());
      const unsigned int particle_size = property_manager->get_particle_size() + integrator->data_length() * sizeof(double);

      // Determine the total number of particles we will send to other processors
      std::vector<int> num_send_particles(world_size);
      std::vector<int> num_recv_particles(world_size);

      std::vector<int> num_send_data(world_size);
      std::vector<int> num_recv_data(world_size);

      std::vector<int> send_offsets(world_size,0);
      std::vector<int> recv_offsets(world_size,0);

      int total_send_data = 0;
      for (types::subdomain_id rank = 0; rank < world_size; ++rank)
        {
          send_offsets[rank] = total_send_data;
          std::pair< const typename std::multimap<types::subdomain_id,Particle <dim> >::const_iterator,
                     const typename std::multimap<types::subdomain_id,Particle <dim> >::const_iterator>
          send_particle_range = send_particles.equal_range(rank);
          num_send_particles[rank] = std::distance(send_particle_range.first,send_particle_range.second);
          num_send_data[rank] = num_send_particles[rank] * particle_size;
          total_send_data += num_send_particles[rank] * particle_size;
        }

      Utilities::MPI::sum (num_send_particles, this->get_mpi_communicator(), num_recv_particles);

      // Allocate space for sending and receiving particle data
      std::vector<char> send_data(send_particles.size() * particle_size);
      void* data = static_cast<void *> (&send_data.front());

      // Copy the particle data into the send array
      for (typename std::multimap<types::subdomain_id,Particle<dim> >::const_iterator particle = send_particles.begin(); particle!=send_particles.end(); ++particle)
        {
          particle->second.write_data(data);
          integrator->write_data(data, particle->second.get_id());
        }

      AssertThrow(data == &(send_data.back())+1,
                  ExcMessage("The amount of data written into the array that is send to other processes "
                             "is inconsistent with the number and size of particles."));

      // Notify other processors how many particles we will be sending
      std::vector<int> recv_offset(world_size,0);

      MPI_Alltoall(&(num_send_data[0]), 1, MPI_INT, &(num_recv_data[0]), 1, MPI_INT, this->get_mpi_communicator());

      int total_recv_data = 0;
      for (unsigned int rank=0; rank<world_size; ++rank)
        {
          recv_offset[rank] = total_recv_data;
          total_recv_data += num_recv_data[rank];
        }

      // Set up the space for the received particle data
      std::vector<char> recv_data(total_recv_data);

      // Exchange the particle data between domains
      MPI_Alltoallv (&(send_data[0]), &(num_send_data[0]), &(send_offsets[0]), MPI_CHAR,
                      &(recv_data[0]), &(num_recv_data[0]), &(recv_offset[0]), MPI_CHAR,
                      this->get_mpi_communicator());

      // Put the received particles into the domain if they are in the triangulation
      void *recv_data_it = static_cast<void*> (&recv_data.front());

      for (int i=0; i<num_recv_particles[self_rank]; ++i)
        {
          Particle<dim> recv_particle(recv_data_it,property_manager->get_data_len());
          integrator->read_data(recv_data_it, recv_particle.get_id());

          typename parallel::distributed::Triangulation<dim>::active_cell_iterator cell;
          try
            {
              cell = (GridTools::find_active_cell_around_point<> (this->get_mapping(), this->get_triangulation(), recv_particle.get_location())).first;
            }
          catch (GridTools::ExcPointNotFound<dim> &)
            {
              // If we can find no cell for this particle it has left the domain due
              // to an integration error or open boundary. Simply ignore the
              // tracer in this case.
              continue;
            }


          if (cell->is_locally_owned())
            {
              const LevelInd found_cell = std::make_pair(cell->level(),cell->index());
              particles.insert(std::make_pair(found_cell, recv_particle));
            }
        }

      AssertThrow(recv_data_it == &recv_data.back()+1,
                  ExcMessage("The amount of data that was read into new particles "
                             "does not match the amount of data sent around."));
    }

    template <int dim>
    void
    World<dim>::get_particle_velocities(std::vector<Tensor<1,dim> > &old_velocities,
                                        std::vector<Tensor<1,dim> > &velocities)
    {
      std::vector<Tensor<1,dim> >  result;
      std::vector<Tensor<1,dim> >  old_result;
      std::vector<Point<dim> >     particle_points;

      const DoFHandler<dim> *dof_handler = &(this->get_dof_handler());
      const typename parallel::distributed::Triangulation<dim> *triangulation = &(this->get_triangulation());

      unsigned int particle_idx = 0;
      // Get the velocity for each cell at a time so we can take advantage of knowing the active cell
      for (typename std::multimap<LevelInd, Particle<dim> >::iterator
           it=particles.begin(); it!=particles.end();)
        {
          // Get the current cell
          const LevelInd cur_cell = it->first;

          // Get the cell the particle is in
          const typename DoFHandler<dim>::active_cell_iterator cell (triangulation, cur_cell.first, cur_cell.second, dof_handler);

          // Resize the vectors to the number of particles in this cell
          const unsigned int num_cell_particles = particles.count(cur_cell);
          particle_points.resize(num_cell_particles);
          result.resize(num_cell_particles);
          old_result.resize(num_cell_particles);

          // Get a vector of the particle locations in this cell
          unsigned int n_particles_in_cell = 0;
          while (it != particles.end() && it->first == cur_cell)
            {
              const Point<dim> position = it->second.get_location();
              particle_points[n_particles_in_cell++] = this->get_mapping().transform_real_to_unit_cell(cell, position);
              it++;
            }

          const std::vector< double > ww(num_cell_particles, 1./((double) num_cell_particles));
          const Quadrature<dim> quadrature_formula(particle_points, ww);
          FEValues<dim> fe_value (this->get_mapping(),
                                  this->get_fe(),
                                  quadrature_formula,
                                  update_values);

          fe_value.reinit (cell);
          fe_value[this->introspection().extractors.velocities].get_function_values (this->get_solution(),
                                                                                    result);
          fe_value[this->introspection().extractors.velocities].get_function_values (this->get_old_solution(),
                                                                                    old_result);

          // Copy the resulting velocities to the appropriate vector#
          for (unsigned int id = 0; id < num_cell_particles; ++id)
            {
              velocities[particle_idx] = result[id];
              old_velocities[particle_idx] = old_result[id];
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
