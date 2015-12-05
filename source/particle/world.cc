/*
 Copyright (C) 2015 by the authors of the ASPECT code.

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
#include <aspect/utilities.h>
#include <aspect/compat.h>
#include <aspect/geometry_model/box.h>

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
    template <>
    World<2>::World()
      :
      global_number_of_particles(0),
      data_offset(numbers::invalid_unsigned_int)
    {

      aspect::internals::SimulatorSignals::register_connector_function_2d (std_cxx11::bind(&World<2>::connect_to_signals,
                                                                           std_cxx11::ref(*this),
                                                                           std_cxx11::_1));
    }

    template <>
    World<3>::World()
      :
      global_number_of_particles(0),
      data_offset(numbers::invalid_unsigned_int)
    {
      aspect::internals::SimulatorSignals::register_connector_function_3d (std_cxx11::bind(&World<3>::connect_to_signals,
                                                                           std_cxx11::ref(*this),
                                                                           std_cxx11::_1));
    }

    template <int dim>
    World<dim>::~World()
    {}

    template <int dim>
    void
    World<dim>::initialize(Generator::Interface<dim> *particle_generator,
                           Integrator::Interface<dim> *particle_integrator,
                           Interpolator::Interface<dim> *property_interpolator,
                           Property::Manager<dim> *manager,
                           const ParticleLoadBalancing &load_balancing,
                           const unsigned int max_part_per_cell,
                           const unsigned int weight)
    {
      generator.reset(particle_generator);
      integrator.reset(particle_integrator);
      interpolator.reset(property_interpolator);
      property_manager.reset(manager);
      particle_load_balancing = load_balancing;
      max_particles_per_cell = max_part_per_cell;
      tracer_weight = weight;

#if DEAL_II_VERSION_GTE(8,4,0)
      if (particle_load_balancing == repartition)
        this->get_triangulation().signals.cell_weight.connect(std_cxx11::bind(&aspect::Particle::World<dim>::cell_weight,
                                                                              std_cxx11::ref(*this),
                                                                              std_cxx11::_1,
                                                                              std_cxx11::_2));
#else
      AssertThrow(particle_load_balancing != repartition,
                  ExcMessage("You tried to select the load balancing strategy 'repartition', "
                             "which is only available for deal.II 8.4 and newer but the installed version "
                             "seems to be older. Please update your deal.II or choose a different strategy."));
#endif
    }

    template <int dim>
    const Property::Manager<dim> &
    World<dim>::get_property_manager() const
    {
      return *property_manager;
    }

    template <int dim>
    std::multimap<types::LevelInd, Particle<dim> > &
    World<dim>::get_particles()
    {
      return particles;
    }

    template <int dim>
    const std::multimap<types::LevelInd, Particle<dim> > &
    World<dim>::get_particles() const
    {
      return particles;
    }

    template <int dim>
    types::particle_index
    World<dim>::n_global_particles() const
    {
      return global_number_of_particles;
    }

    template <int dim>
    void
    World<dim>::update_n_global_particles()
    {
      global_number_of_particles = dealii::Utilities::MPI::sum (particles.size(), this->get_mpi_communicator());
    }

    template <int dim>
    unsigned int
    World<dim>::get_global_max_tracers_per_cell() const
    {
      unsigned int local_max_tracer_per_cell(0);
      typename parallel::distributed::Triangulation<dim>::active_cell_iterator cell = this->get_triangulation().begin_active();
      for (; cell!=this->get_triangulation().end(); ++cell)
        if (cell->is_locally_owned())
          {
            const types::LevelInd found_cell = std::make_pair<int, int> (cell->level(),cell->index());
            const std::pair<typename std::multimap<types::LevelInd, Particle<dim> >::const_iterator, typename std::multimap<types::LevelInd, Particle<dim> >::const_iterator> particles_in_cell
              = particles.equal_range(found_cell);
            const unsigned int tracers_in_cell = std::distance(particles_in_cell.first,particles_in_cell.second);
            local_max_tracer_per_cell = std::max(local_max_tracer_per_cell,
                                                 tracers_in_cell);
          }

      return dealii::Utilities::MPI::max(local_max_tracer_per_cell,this->get_mpi_communicator());
    }

    template <int dim>
    void
    World<dim>::connect_to_signals(aspect::SimulatorSignals<dim> &signals)
    {
      signals.pre_refinement_store_user_data.connect(std_cxx11::bind(&World<dim>::register_store_callback_function,
                                                                     std_cxx11::ref(*this),
                                                                     std_cxx11::_1));
      signals.post_refinement_load_user_data.connect(std_cxx11::bind(&World<dim>::register_load_callback_function,
                                                                     std_cxx11::ref(*this),
                                                                     std_cxx11::_1));
    }

    template <int dim>
    void
    World<dim>::register_store_callback_function(typename parallel::distributed::Triangulation<dim> &triangulation)
    {
      // Only save and load tracers if there are any, we might get here for
      // example before the tracer generation in timestep 0, or if somebody
      // selected the tracer postprocessor but generated 0 tracers
      const unsigned int max_tracers_per_cell = get_global_max_tracers_per_cell();

      if (max_tracers_per_cell > 0)
        {
          const std_cxx11::function<void(const typename parallel::distributed::Triangulation<dim>::cell_iterator &,
                                         const typename parallel::distributed::Triangulation<dim>::CellStatus, void *) > callback_function
            = std_cxx11::bind(&aspect::Particle::World<dim>::store_tracers,
                              std_cxx11::ref(*this),
                              std_cxx11::_1,
                              std_cxx11::_2,
                              std_cxx11::_3);

          // We need to transfer the number of tracers for this cell and
          // the tracer data itself and we need to provide 2^dim times the
          // space for the data in case a cell is coarsened
          const std::size_t transfer_size_per_cell = sizeof (unsigned int) +
                                                     (property_manager->get_particle_size() * max_tracers_per_cell)
                                                     *  std::pow(2,dim);
          data_offset = triangulation.register_data_attach(transfer_size_per_cell,callback_function);
        }
    }

    template <int dim>
    void
    World<dim>::register_load_callback_function(typename parallel::distributed::Triangulation<dim> &triangulation)
    {
      // All particles have been stored, when we reach this point. Empty the
      // map and fill with new particles.
      particles.clear();

      // Check if something was stored and load it
      if (data_offset != numbers::invalid_unsigned_int)
        {
          const std_cxx11::function<void(const typename parallel::distributed::Triangulation<dim>::cell_iterator &,
                                         const typename parallel::distributed::Triangulation<dim>::CellStatus,
                                         const void *) > callback_function
            = std_cxx11::bind(&aspect::Particle::World<dim>::load_tracers,
                              std_cxx11::ref(*this),
                              std_cxx11::_1,
                              std_cxx11::_2,
                              std_cxx11::_3);

          triangulation.notify_ready_to_unpack(data_offset,callback_function);

          // Reset offset and update global number of particles. The number
          // can change because of discarded or newly generated particles
          data_offset = numbers::invalid_unsigned_int;
          update_n_global_particles();
        }
    }

    template <int dim>
    unsigned int
    World<dim>::cell_weight(const typename parallel::distributed::Triangulation<dim>::cell_iterator &cell,
                            const typename parallel::distributed::Triangulation<dim>::CellStatus status)
    {
      if (cell->active() && !cell->is_locally_owned())
        return 0;

      if (status == parallel::distributed::Triangulation<dim>::CELL_PERSIST
          || status == parallel::distributed::Triangulation<dim>::CELL_REFINE)
        {
          const types::LevelInd found_cell = std::make_pair<int, int> (cell->level(),cell->index());
          const unsigned int n_particles_in_cell = particles.count(found_cell);
          return n_particles_in_cell * tracer_weight;
        }
      else if (status == parallel::distributed::Triangulation<dim>::CELL_COARSEN)
        {
          unsigned int n_particles_in_cell = 0;

          for (unsigned int child_index = 0; child_index < GeometryInfo<dim>::max_children_per_cell; ++child_index)
            {
              const typename parallel::distributed::Triangulation<dim>::cell_iterator child = cell->child(child_index);
              const types::LevelInd found_cell = std::make_pair<int, int> (child->level(),child->index());
              n_particles_in_cell += particles.count(found_cell);
            }
          return n_particles_in_cell * tracer_weight;
        }

      Assert (false, ExcInternalError());
      return 0;
    }

    template <int dim>
    void
    World<dim>::store_tracers(const typename parallel::distributed::Triangulation<dim>::cell_iterator &cell,
                              const typename parallel::distributed::Triangulation<dim>::CellStatus status,
                              void *data)
    {
      unsigned int n_particles_in_cell(0);

      // If the cell persist or is refined store all tracers of the current cell.
      if (status == parallel::distributed::Triangulation<dim>::CELL_PERSIST
          || status == parallel::distributed::Triangulation<dim>::CELL_REFINE)
        {
          const types::LevelInd found_cell = std::make_pair<int, int> (cell->level(),cell->index());
          const std::pair<typename std::multimap<types::LevelInd, Particle<dim> >::iterator, typename std::multimap<types::LevelInd, Particle<dim> >::iterator> particles_in_cell
            = particles.equal_range(found_cell);
          n_particles_in_cell = std::distance(particles_in_cell.first,particles_in_cell.second);

          unsigned int *ndata = static_cast<unsigned int *> (data);
          *ndata = n_particles_in_cell;
          data = static_cast<void *> (ndata + 1);

          for (typename std::multimap<types::LevelInd, Particle<dim> >::iterator particle = particles_in_cell.first;
               particle != particles_in_cell.second; ++particle)
            {
              particle->second.write_data(data);
            }
        }
      // If this cell is the parent of children that will be coarsened, collect
      // the tracers of all children.
      // First check if the maximum number of particles per cell is exceeded for
      // the new cell, and if that is the case, only store every 2^dim 'th
      // particle.
      else if (status == parallel::distributed::Triangulation<dim>::CELL_COARSEN)
        {
          for (unsigned int child_index = 0; child_index < GeometryInfo<dim>::max_children_per_cell; ++child_index)
            {
              const typename parallel::distributed::Triangulation<dim>::cell_iterator child = cell->child(child_index);
              const types::LevelInd found_cell = std::make_pair<int, int> (child->level(),child->index());
              const std::pair<typename std::multimap<types::LevelInd, Particle<dim> >::iterator, typename std::multimap<types::LevelInd, Particle<dim> >::iterator> particles_in_cell
                = particles.equal_range(found_cell);
              n_particles_in_cell += std::distance(particles_in_cell.first,particles_in_cell.second);
            }

          bool reduce_tracers = false;
          if (particle_load_balancing == remove_particles || particle_load_balancing == remove_and_add_particles)
            if ((max_particles_per_cell > 0) && (n_particles_in_cell > max_particles_per_cell))
              {
                n_particles_in_cell /= GeometryInfo<dim>::max_children_per_cell;
                reduce_tracers = true;
              }

          unsigned int *ndata = static_cast<unsigned int *> (data);
          *ndata = n_particles_in_cell;

          data = static_cast<void *> (ndata + 1);

          unsigned int particle_index = 0;

          for (unsigned int child_index = 0; child_index < GeometryInfo<dim>::max_children_per_cell; ++child_index)
            {
              const typename parallel::distributed::Triangulation<dim>::cell_iterator child = cell->child(child_index);
              const types::LevelInd found_cell = std::make_pair<int, int> (child->level(),child->index());
              const std::pair<typename std::multimap<types::LevelInd, Particle<dim> >::iterator, typename std::multimap<types::LevelInd, Particle<dim> >::iterator>
              particles_in_cell = particles.equal_range(found_cell);

              for (typename std::multimap<types::LevelInd, Particle<dim> >::iterator particle = particles_in_cell.first;
                   particle != particles_in_cell.second; ++particle, ++particle_index)
                {
                  if (!reduce_tracers || (particle_index % GeometryInfo<dim>::max_children_per_cell == 0))
                    particle->second.write_data(data);
                }
            }
        }
      else
        Assert (false, ExcInternalError());

    }

    template <int dim>
    void
    World<dim>::load_tracers(const typename parallel::distributed::Triangulation<dim>::cell_iterator &cell,
                             const typename parallel::distributed::Triangulation<dim>::CellStatus status,
                             const void *data)
    {
      const unsigned int *n_particles_in_cell_ptr = static_cast<const unsigned int *> (data);
      const void *pdata = reinterpret_cast<const void *> (n_particles_in_cell_ptr + 1);

      // Load all particles from the data stream and store them in the local
      // particle map.
      for (unsigned int i = 0; i < *n_particles_in_cell_ptr; ++i)
        {
          Particle<dim> p(pdata,property_manager->get_particle_size());

          if (status == parallel::distributed::Triangulation<dim>::CELL_COARSEN
              || status == parallel::distributed::Triangulation<dim>::CELL_PERSIST)
            particles.insert(std::make_pair(std::make_pair(cell->level(),cell->index()),p));
          else if (status == parallel::distributed::Triangulation<dim>::CELL_REFINE)
            {
              for (unsigned int child_index = 0; child_index < GeometryInfo<dim>::max_children_per_cell; ++child_index)
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

      if (particle_load_balancing == remove_and_add_particles)
        {
          // TODO: Add particles
        }
    }

    template <int dim>
    void
    World<dim>::sort_particles_in_subdomains_and_cells()
    {
      // TODO: The current algorithm only works for CFL numbers <= 1.0,
      // because it only knows the subdomain_id of ghost cells, but not
      // of artificial cells.

      // There are three reasons why a particle is not in its old cell:
      // It moved to another cell, to another domain or it left the mesh.
      // Sort the particles accordingly and deal with them
      std::multimap<types::LevelInd, Particle<dim> >     moved_particles_cell;
      std::multimap<types::subdomain_id, Particle<dim> > moved_particles_domain;
      std::multimap<types::LevelInd, Particle<dim> >     lost_particles;

      // Find the cells that the particles moved to.
      // Note that the iterator in the following loop is increased in a
      // very particular way, because it is changed, if elements
      // get erased. A change can result in invalid memory access.
      typename std::multimap<types::LevelInd, Particle<dim> >::iterator   it;
      for (it=particles.begin(); it!=particles.end();)
        {
          // If we know the particle's old cell, check if it is still inside
          if (it->first != std::make_pair(-1,-1))
            {
              const typename parallel::distributed::Triangulation<dim>::active_cell_iterator
              old_cell (&(this->get_triangulation()), it->first.first, it->first.second);

              try
                {
                  const Point<dim> p_unit = this->get_mapping().transform_real_to_unit_cell(old_cell, it->second.get_location());
                  if (GeometryInfo<dim>::is_inside_unit_cell(p_unit))
                    {
                      // The particle is still in the old cell, move on to next particle
                      ++it;
                      continue;
                    }
                }
              catch (...)
                {}
            }

          // The particle is not in its old cell. Look for the new cell.
          typename parallel::distributed::Triangulation<dim>::active_cell_iterator cell;
          try
            {
              cell = (GridTools::find_active_cell_around_point<> (this->get_mapping(), this->get_triangulation(), it->second.get_location())).first;
            }
          catch (GridTools::ExcPointNotFound<dim> &)
            {
              // If we can find no cell for this particle it has left the domain due
              // to an integration error or an open boundary.
              lost_particles.insert(*it);

              // Now remove the lost particle and continue with next particle.
              // Also make sure we do not invalidate the iterator we are increasing.
              const typename std::multimap<types::LevelInd, Particle<dim> >::iterator particle_to_delete = it;
              it++;
              particles.erase(particle_to_delete);
              continue;
            }


          // Reinsert the particle into our domain if we found its cell
          // Mark it for MPI transfer otherwise
          if (cell->is_locally_owned())
            {
              const types::LevelInd found_cell = std::make_pair(cell->level(),cell->index());
              moved_particles_cell.insert(std::make_pair(found_cell, it->second));
            }
          else
            moved_particles_domain.insert(std::make_pair(cell->subdomain_id(),it->second));

          // Now remove the lost particle and continue with next particle.
          // Also make sure we do not invalidate the iterator we are increasing.
          const typename std::multimap<types::LevelInd, Particle<dim> >::iterator particle_to_delete = it;
          it++;
          particles.erase(particle_to_delete);
        }

      // If particles fell out of the mesh, put them back in if they have crossed
      // a periodic boundary. If they have left the mesh otherwise, they will be
      // discarded by being deleted from lost_particles, and not inserted anywhere.
      move_particles_back_into_mesh(lost_particles,
                                    moved_particles_cell,
                                    moved_particles_domain);

      // Reinsert all local particles with their cells
      particles.insert(moved_particles_cell.begin(),moved_particles_cell.end());

      // Swap lost particles between processors if we have more than one process
      if (dealii::Utilities::MPI::n_mpi_processes(this->get_mpi_communicator()) > 1)
        send_recv_particles(moved_particles_domain);
    }

    template <int dim>
    void
    World<dim>::move_particles_back_into_mesh(std::multimap<types::LevelInd, Particle<dim> >          &lost_particles,
                                              std::multimap<types::LevelInd, Particle<dim> >            &moved_particles_cell,
                                              std::multimap<types::subdomain_id, Particle<dim> >        &moved_particles_domain)
    {
      // TODO: fix this to work with arbitrary meshes. Currently periodic boundaries only work for boxes.
      // If the geometry is not a box, we simply discard particles that have left the
      // model domain.

      const GeometryModel::Box<dim> *geometry
        = dynamic_cast<const GeometryModel::Box<dim>*> (&this->get_geometry_model());

      if (geometry != 0)
        {
          const Point<dim> origin = geometry->get_origin();
          const Point<dim> extent = geometry->get_extents();
          const std::set< std::pair< std::pair<types::boundary_id, types::boundary_id>, unsigned int> > periodic_boundaries =
            geometry->get_periodic_boundary_pairs();

          std::vector<bool> periodic(dim,false);
          std::set< std::pair< std::pair<types::boundary_id, types::boundary_id>, unsigned int> >::const_iterator boundary =
            periodic_boundaries.begin();
          for (; boundary != periodic_boundaries.end(); ++boundary)
            periodic[boundary->second] = true;

          typename std::multimap<types::LevelInd, Particle<dim> >::iterator lost_particle = lost_particles.begin();
          for (; lost_particle != lost_particles.end(); ++lost_particle)
            {
              // modify the particle position if it crossed a periodic boundary
              Point<dim> particle_position = lost_particle->second.get_location();
              for (unsigned int i = 0; i < dim; ++i)
                {
                  if (periodic[i])
                    {
                      if (particle_position[i] < origin[i])
                        particle_position[i] += extent[i];
                      else if (particle_position[i] > origin[i] + extent[i])
                        particle_position[i] -= extent[i];
                    }
                }
              lost_particle->second.set_location(particle_position);

              // Try again looking for the new cell with the updated position
              typename parallel::distributed::Triangulation<dim>::active_cell_iterator cell;
              try
                {
                  cell = (GridTools::find_active_cell_around_point<> (this->get_mapping(), this->get_triangulation(), lost_particle->second.get_location())).first;
                }
              catch (GridTools::ExcPointNotFound<dim> &)
                {
                  // If we can find no cell for this particle there is no hope left
                  // to find its cell. Simply delete the particle.
                  continue;
                }

              // Reinsert the particle into our domain if we found its cell
              // Mark it for MPI transfer otherwise
              if (cell->is_locally_owned())
                {
                  const types::LevelInd found_cell = std::make_pair(cell->level(),cell->index());
                  moved_particles_cell.insert(std::make_pair(found_cell, lost_particle->second));
                }
              else
                moved_particles_domain.insert(std::make_pair(cell->subdomain_id(),lost_particle->second));
            }
        }
    }

    template <int dim>
    void
    World<dim>::send_recv_particles(const std::multimap<types::subdomain_id,Particle <dim> > &send_particles)
    {
      // Determine the communication pattern
#if DEAL_II_VERSION_GTE(8,4,0)
      const std::vector<types::subdomain_id> neighbors (this->get_triangulation().ghost_owners().begin(),
                                                        this->get_triangulation().ghost_owners().end());
#else
      const std::set<types::subdomain_id> ghost_owner = ghost_owners(this->get_triangulation());
      const std::vector<types::subdomain_id> neighbors (ghost_owner.begin(),
                                                        ghost_owner.end());
#endif
      const unsigned int n_neighbors = neighbors.size();
      const unsigned int particle_size = property_manager->get_particle_size() + integrator->get_data_size();

      // Determine the amount of data we will send to other processors
      std::vector<int> n_send_data(n_neighbors);
      std::vector<int> n_recv_data(n_neighbors);

      std::vector<int> send_offsets(n_neighbors);
      std::vector<int> recv_offsets(n_neighbors);

      // Allocate space for sending and receiving particle data
      std::vector<char> send_data(send_particles.size() * particle_size);
      void *data = static_cast<void *> (&send_data.front());

      int total_send_data = 0;
      for (types::subdomain_id neighbor_id = 0; neighbor_id < n_neighbors; ++neighbor_id)
        {
          send_offsets[neighbor_id] = total_send_data;

          std::pair< const typename std::multimap<types::subdomain_id,Particle <dim> >::const_iterator,
              const typename std::multimap<types::subdomain_id,Particle <dim> >::const_iterator>
              send_particle_range = send_particles.equal_range(neighbors[neighbor_id]);

          int n_send_particles = std::distance(send_particle_range.first,send_particle_range.second);
          n_send_data[neighbor_id] = n_send_particles * particle_size;
          total_send_data += n_send_particles * particle_size;

          // Copy the particle data into the send array
          typename std::multimap<types::subdomain_id,Particle<dim> >::const_iterator particle = send_particle_range.first;
          for (; particle != send_particle_range.second; ++particle)
            {
              particle->second.write_data(data);
              data = integrator->write_data(data, particle->second.get_id());
            }
        }

      AssertThrow(data == &(send_data.back())+1,
                  ExcMessage("The amount of data written into the array that is send to other processes "
                             "is inconsistent with the number and size of particles."));

      // Notify other processors how many particles we will send
      std::vector<MPI_Request> n_requests(2*n_neighbors);
      for (unsigned int i=0; i<n_neighbors; ++i)
        MPI_Irecv(&(n_recv_data[i]), 1, MPI_INT, neighbors[i], 0, this->get_mpi_communicator(), &(n_requests[2*i]));
      for (unsigned int i=0; i<n_neighbors; ++i)
        MPI_Isend(&(n_send_data[i]), 1, MPI_INT, neighbors[i], 0, this->get_mpi_communicator(), &(n_requests[2*i+1]));
      MPI_Waitall(2*n_neighbors,&n_requests[0],MPI_STATUSES_IGNORE);

      // Determine how many particles and data we will receive
      int total_recv_data = 0;
      for (unsigned int neighbor_id=0; neighbor_id<n_neighbors; ++neighbor_id)
        {
          recv_offsets[neighbor_id] = total_recv_data;
          total_recv_data += n_recv_data[neighbor_id];
        }
      const int n_recv_particles = total_recv_data / particle_size;

      // Set up the space for the received particle data
      std::vector<char> recv_data(total_recv_data);

      // Exchange the particle data between domains
      std::vector<MPI_Request> requests(2*n_neighbors);
      unsigned int send_ops = 0;
      unsigned int recv_ops = 0;

      for (unsigned int i=0; i<n_neighbors; ++i)
        if (n_recv_data[i] > 0)
          {
            MPI_Irecv(&(recv_data[recv_offsets[i]]), n_recv_data[i], MPI_CHAR, neighbors[i], 1, this->get_mpi_communicator(),&(requests[send_ops]));
            send_ops++;
          }

      for (unsigned int i=0; i<n_neighbors; ++i)
        if (n_send_data[i] > 0)
          {
            MPI_Isend(&(send_data[send_offsets[i]]), n_send_data[i], MPI_CHAR, neighbors[i], 1, this->get_mpi_communicator(),&(requests[send_ops+recv_ops]));
            recv_ops++;
          }
      MPI_Waitall(send_ops+recv_ops,&requests[0],MPI_STATUSES_IGNORE);

      // Put the received particles into the domain if they are in the triangulation
      const void *recv_data_it = static_cast<const void *> (&recv_data.front());

      for (int i=0; i<n_recv_particles; ++i)
        {
          const Particle<dim> recv_particle(recv_data_it,property_manager->get_particle_size());
          recv_data_it = integrator->read_data(recv_data_it, recv_particle.get_id());

          const typename parallel::distributed::Triangulation<dim>::active_cell_iterator cell =
            (GridTools::find_active_cell_around_point<> (this->get_mapping(), this->get_triangulation(), recv_particle.get_location())).first;

          Assert(cell->is_locally_owned(),
                 ExcMessage("Another process sent us a particle, but the particle is not in our domain."));

          const types::LevelInd found_cell = std::make_pair(cell->level(),cell->index());

          if (particle_load_balancing == remove_particles || particle_load_balancing == remove_and_add_particles)
            {
              // Detect if we need to reduce the number of tracers in this cell,
              // we first reduce the incoming tracers, because they likely came from
              // a region, where the particle density is higher than in this cell
              // (otherwise this would not have been triggered).
              const bool reduce_particles = (max_particles_per_cell > 0) && (particles.count(found_cell) >= max_particles_per_cell);

              if ( !reduce_particles || (i % GeometryInfo<dim>::max_children_per_cell == 0))
                particles.insert(std::make_pair(found_cell, recv_particle));
            }
          else
            particles.insert(std::make_pair(found_cell, recv_particle));
        }

      AssertThrow(recv_data_it == &recv_data.back()+1,
                  ExcMessage("The amount of data that was read into new particles "
                             "does not match the amount of data sent around."));
    }

    template <int dim>
    void
    World<dim>::local_initialize_particles(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                           const typename std::multimap<types::LevelInd, Particle<dim> >::iterator &begin_particle,
                                           const typename std::multimap<types::LevelInd, Particle<dim> >::iterator &end_particle)
    {
      const unsigned int particles_in_cell = std::distance(begin_particle,end_particle);
      const unsigned int solution_components = this->introspection().n_components;

      Vector<double> value(solution_components);
      std::vector<Tensor<1,dim> > gradient (solution_components,Tensor<1,dim>());

      std::vector<Vector<double> >  values(particles_in_cell,value);
      std::vector<std::vector<Tensor<1,dim> > > gradients(particles_in_cell,gradient);

      std::vector<Point<dim> >     particle_points(particles_in_cell);

      typename std::multimap<types::LevelInd, Particle<dim> >::iterator it = begin_particle;
      for (unsigned int i = 0; it!=end_particle; ++it,++i)
        {
          const Point<dim> position = it->second.get_location();
          particle_points[i] = this->get_mapping().transform_real_to_unit_cell(cell, position);
        }

      const Quadrature<dim> quadrature_formula(particle_points);
      FEValues<dim> fe_value (this->get_mapping(),
                              this->get_fe(),
                              quadrature_formula,
                              update_values |
                              update_gradients);

      fe_value.reinit (cell);
      fe_value.get_function_values (this->get_solution(),
                                    values);
      fe_value.get_function_gradients (this->get_solution(),
                                       gradients);

      it = begin_particle;
      for (unsigned int i = 0; it!=end_particle; ++it,++i)
        {
          property_manager->initialize_one_particle(it->second,
                                                    values[i],
                                                    gradients[i]);
        }
    }

    template <int dim>
    void
    World<dim>::local_update_particles(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                       const typename std::multimap<types::LevelInd, Particle<dim> >::iterator &begin_particle,
                                       const typename std::multimap<types::LevelInd, Particle<dim> >::iterator &end_particle)
    {
      const unsigned int particles_in_cell = std::distance(begin_particle,end_particle);
      const unsigned int solution_components = this->introspection().n_components;

      Vector<double> value(solution_components);
      std::vector<Tensor<1,dim> > gradient (solution_components,Tensor<1,dim>());

      std::vector<Vector<double> >  values(particles_in_cell,value);
      std::vector<std::vector<Tensor<1,dim> > > gradients(particles_in_cell,gradient);

      std::vector<Point<dim> >     particle_points(particles_in_cell);

      typename std::multimap<types::LevelInd, Particle<dim> >::iterator it = begin_particle;
      for (unsigned int i = 0; it!=end_particle; ++it,++i)
        {
          const Point<dim> position = it->second.get_location();
          particle_points[i] = this->get_mapping().transform_real_to_unit_cell(cell, position);
        }

      const Quadrature<dim> quadrature_formula(particle_points);
      FEValues<dim> fe_value (this->get_mapping(),
                              this->get_fe(),
                              quadrature_formula,
                              update_values |
                              update_gradients);

      fe_value.reinit (cell);
      fe_value.get_function_values (this->get_solution(),
                                    values);
      fe_value.get_function_gradients (this->get_solution(),
                                       gradients);

      it = begin_particle;
      for (unsigned int i = 0; it!=end_particle; ++it,++i)
        {
          property_manager->update_one_particle(it->second,
                                                values[i],
                                                gradients[i]);
        }
    }

    template <int dim>
    void
    World<dim>::local_advect_particles(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                       const typename std::multimap<types::LevelInd, Particle<dim> >::iterator &begin_particle,
                                       const typename std::multimap<types::LevelInd, Particle<dim> >::iterator &end_particle)
    {
      const unsigned int particles_in_cell = std::distance(begin_particle,end_particle);

      std::vector<Tensor<1,dim> >  result(particles_in_cell);
      std::vector<Tensor<1,dim> >  old_result(particles_in_cell);
      std::vector<Point<dim> >     particle_points(particles_in_cell);

      typename std::multimap<types::LevelInd, Particle<dim> >::iterator it = begin_particle;
      for (unsigned int i = 0; it!=end_particle; ++it,++i)
        {
          const Point<dim> position = it->second.get_location();
          particle_points[i] = this->get_mapping().transform_real_to_unit_cell(cell, position);
        }

      const Quadrature<dim> quadrature_formula(particle_points);
      FEValues<dim> fe_value (this->get_mapping(),
                              this->get_fe(),
                              quadrature_formula,
                              update_values);

      fe_value.reinit (cell);
      fe_value[this->introspection().extractors.velocities].get_function_values (this->get_solution(),
                                                                                 result);
      fe_value[this->introspection().extractors.velocities].get_function_values (this->get_old_solution(),
                                                                                 old_result);

      integrator->local_integrate_step(begin_particle,
                                       end_particle,
                                       old_result,
                                       result,
                                       this->get_old_timestep());
    }

    template <int dim>
    void
    World<dim>::generate_particles()
    {
      generator->generate_particles(particles);
      update_n_global_particles();
    }

    template <int dim>
    void
    World<dim>::initialize_particles()
    {
      // TODO: Change this loop over all cells to use the WorkStream interface

      if (property_manager->get_n_property_components() > 0)
        {
          // Loop over all cells and initialize the particles cell-wise
          typename DoFHandler<dim>::active_cell_iterator
          cell = this->get_dof_handler().begin_active(),
          endc = this->get_dof_handler().end();

          for (; cell!=endc; ++cell)
            if (cell->is_locally_owned())
              {
                std::pair< const typename std::multimap<types::LevelInd,Particle <dim> >::iterator,
                    const typename std::multimap<types::LevelInd,Particle <dim> >::iterator>
                    particle_range_in_cell = particles.equal_range(std::make_pair(cell->level(),cell->index()));

                // Only initialize particles, if there are any in this cell
                if (particle_range_in_cell.first != particle_range_in_cell.second)
                  local_initialize_particles(cell,
                                             particle_range_in_cell.first,
                                             particle_range_in_cell.second);
              }
        }
    }

    template <int dim>
    void
    World<dim>::update_particles()
    {
      // TODO: Change this loop over all cells to use the WorkStream interface

      if (property_manager->get_n_property_components() > 0)
        {
          // Loop over all cells and update the particles cell-wise
          typename DoFHandler<dim>::active_cell_iterator
          cell = this->get_dof_handler().begin_active(),
          endc = this->get_dof_handler().end();

          for (; cell!=endc; ++cell)
            if (cell->is_locally_owned())
              {
                std::pair< const typename std::multimap<types::LevelInd,Particle <dim> >::iterator,
                    const typename std::multimap<types::LevelInd,Particle <dim> >::iterator>
                    particle_range_in_cell = particles.equal_range(std::make_pair(cell->level(),cell->index()));

                // Only update particles, if there are any in this cell
                if (particle_range_in_cell.first != particle_range_in_cell.second)
                  local_update_particles(cell,
                                         particle_range_in_cell.first,
                                         particle_range_in_cell.second);
              }
        }
    }

    template <int dim>
    void
    World<dim>::advect_particles()
    {
      // TODO: Change this loop over all cells to use the WorkStream interface

      // Loop over all cells and advect the particles cell-wise
      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();

      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
          {
            std::pair< const typename std::multimap<types::LevelInd,Particle <dim> >::iterator,
                const typename std::multimap<types::LevelInd,Particle <dim> >::iterator>
                particle_range_in_cell = particles.equal_range(std::make_pair(cell->level(),cell->index()));

            // Only advect particles, if there are any in this cell
            if (particle_range_in_cell.first != particle_range_in_cell.second)
              local_advect_particles(cell,
                                     particle_range_in_cell.first,
                                     particle_range_in_cell.second);
          }
    }

    template <int dim>
    void
    World<dim>::advance_timestep()
    {
      do
        {
          advect_particles();

          // Find the cells that the particles moved to
          sort_particles_in_subdomains_and_cells();
        }
      // Keep calling the integrator until it indicates it is finished
      while (integrator->new_integration_step());

      // Update particle properties
      if (property_manager->need_update() == Property::update_time_step)
        update_particles();

      // Update the number of global particles if some have left the domain
      update_n_global_particles();
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
