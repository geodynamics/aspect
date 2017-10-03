/*
  Copyright (C) 2015 - 2017 by the authors of the ASPECT code.

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
    template <int dim>
    World<dim>::World()
    {}

    template <int dim>
    World<dim>::~World()
    {}

    template <int dim>
    void
    World<dim>::initialize()
    {
      if (particle_load_balancing & ParticleLoadBalancing::repartition)
        this->get_triangulation().signals.cell_weight.connect(std_cxx11::bind(&aspect::Particle::World<dim>::cell_weight,
                                                                              std_cxx11::ref(*this),
                                                                              std_cxx11::_1,
                                                                              std_cxx11::_2));

      // Create a particle handler that stores the future particles.
      // If we restarted from a checkpoint we will fill this particle handler
      // later with its serialized variables and stored particles
      particle_handler.reset(new ParticleHandler<dim>(this->get_triangulation(),
                                                      this->get_mapping(),
                                                      this->get_mpi_communicator(),
                                                      property_manager->get_n_property_components()));

      const std_cxx11::function<std::size_t ()> size_callback_function
        = std_cxx11::bind(&aspect::Particle::Integrator::Interface<dim>::get_data_size,
                          std_cxx11::ref(*integrator));

      const std_cxx11::function<void *(const typename ParticleHandler<dim>::particle_iterator &,
                                       void *)> store_callback_function
        = std_cxx11::bind(&aspect::Particle::Integrator::Interface<dim>::write_data,
                          std_cxx11::ref(*integrator),
                          std_cxx11::_1,
                          std_cxx11::_2);

      const std_cxx11::function<const void *(const typename ParticleHandler<dim>::particle_iterator &,
                                             const void *)> load_callback_function
        = std_cxx11::bind(&aspect::Particle::Integrator::Interface<dim>::read_data,
                          std_cxx11::ref(*integrator),
                          std_cxx11::_1,
                          std_cxx11::_2);

      particle_handler->register_additional_store_load_functions(size_callback_function,
                                                                 store_callback_function,
                                                                 load_callback_function);

      connect_to_signals(this->get_signals());
    }

    template <int dim>
    const Property::Manager<dim> &
    World<dim>::get_property_manager() const
    {
      return *property_manager;
    }

    template <int dim>
    const ParticleHandler<dim> &
    World<dim>::get_particle_handler() const
    {
      return *particle_handler.get();
    }

    template <int dim>
    const Interpolator::Interface<dim> &
    World<dim>::get_interpolator() const
    {
      return *interpolator;
    }

    template <int dim>
    std::string
    World<dim>::generate_output() const
    {
      // If we do not write output
      // return early with the number of particles that were advected
      if (!output)
        return "";

      TimerOutput::Scope timer_section(this->get_computing_timer(), "Particles: Output");
      const double output_time = (this->convert_output_to_years() ?
                                  this->get_time() / year_in_seconds :
                                  this->get_time());

      const std::string filename = output->output_particle_data(*particle_handler,
                                                                property_manager->get_data_info(),
                                                                output_time);

      return filename;
    }

    template <int dim>
    types::particle_index
    World<dim>::n_global_particles() const
    {
      return particle_handler->n_global_particles();
    }


    template <int dim>
    void
    World<dim>::connect_to_signals(aspect::SimulatorSignals<dim> &signals)
    {
      signals.post_set_initial_state.connect(std_cxx11::bind(&World<dim>::setup_initial_state,
                                                             std_cxx11::ref(*this)));

      signals.pre_refinement_store_user_data.connect(std_cxx11::bind(&ParticleHandler<dim>::register_store_callback_function,
                                                                     std_cxx11::ref(*particle_handler),
                                                                     false));

      signals.post_refinement_load_user_data.connect(std_cxx11::bind(&ParticleHandler<dim>::register_load_callback_function,
                                                                     std_cxx11::ref(*particle_handler),
                                                                     false));

      signals.pre_checkpoint_store_user_data.connect(std_cxx11::bind(&ParticleHandler<dim>::register_store_callback_function,
                                                                     std_cxx11::ref(*particle_handler),
                                                                     true));

      signals.post_resume_load_user_data.connect(std_cxx11::bind(&ParticleHandler<dim>::register_load_callback_function,
                                                                 std_cxx11::ref(*particle_handler),
                                                                 true));

      signals.post_refinement_load_user_data.connect(std_cxx11::bind(&World<dim>::apply_particle_per_cell_bounds,
                                                                     std_cxx11::ref(*this)));
      signals.post_resume_load_user_data.connect(std_cxx11::bind(&World<dim>::apply_particle_per_cell_bounds,
                                                                 std_cxx11::ref(*this)));

      if (update_ghost_particles &&
          dealii::Utilities::MPI::n_mpi_processes(this->get_mpi_communicator()) > 1)
        {
          signals.post_refinement_load_user_data.connect(std_cxx11::bind(&ParticleHandler<dim>::exchange_ghost_particles,
                                                                         std_cxx11::ref(*particle_handler)));
          signals.post_resume_load_user_data.connect(std_cxx11::bind(&ParticleHandler<dim>::exchange_ghost_particles,
                                                                     std_cxx11::ref(*particle_handler)));
        }
    }



    template <int dim>
    void
    World<dim>::apply_particle_per_cell_bounds()
    {
      // If any load balancing technique is selected that creates/destroys particles
      if (particle_load_balancing & ParticleLoadBalancing::remove_and_add_particles)
        {
          // First do some preparation for particle generation in poorly
          // populated areas. For this we need to know which particle ids to
          // generate so that they are globally unique.
          // Ensure this by communicating the number of particles that every
          // process is going to generate.
          types::particle_index local_next_particle_index = particle_handler->next_free_particle_index;
          if (particle_load_balancing & ParticleLoadBalancing::add_particles)
            {
              types::particle_index particles_to_add_locally = 0;

              // Loop over all cells and determine the number of particles to generate
              typename DoFHandler<dim>::active_cell_iterator
              cell = this->get_dof_handler().begin_active(),
              endc = this->get_dof_handler().end();

              for (; cell!=endc; ++cell)
                if (cell->is_locally_owned())
                  {
                    const unsigned int particles_in_cell = particle_handler->n_particles_in_cell(cell);

                    if (particles_in_cell < min_particles_per_cell)
                      particles_to_add_locally += static_cast<types::particle_index> (min_particles_per_cell - particles_in_cell);
                  }

              // Determine the starting particle index of this process, which
              // is the highest currently existing particle index plus the sum
              // of the number of newly generated particles of all
              // processes with a lower rank.

              types::particle_index local_start_index = 0.0;
              MPI_Scan(&particles_to_add_locally, &local_start_index, 1, ASPECT_PARTICLE_INDEX_MPI_TYPE, MPI_SUM, this->get_mpi_communicator());
              local_start_index -= particles_to_add_locally;
              local_next_particle_index += local_start_index;

              const types::particle_index globally_generated_particles =
                dealii::Utilities::MPI::sum(particles_to_add_locally,this->get_mpi_communicator());

              AssertThrow (particle_handler->next_free_particle_index
                           <= std::numeric_limits<types::particle_index>::max() - globally_generated_particles,
                           ExcMessage("There is no free particle index left to generate a new particle id. Please check if your "
                                      "model generates unusually many new particles (by repeatedly deleting and regenerating particles), or "
                                      "recompile deal.II with the DEAL_II_WITH_64BIT_INDICES option enabled, to use 64-bit integers for "
                                      "particle ids."));

              particle_handler->next_free_particle_index += globally_generated_particles;
            }

          boost::mt19937 random_number_generator;

          // Loop over all cells and generate or remove the particles cell-wise
          typename DoFHandler<dim>::active_cell_iterator
          cell = this->get_dof_handler().begin_active(),
          endc = this->get_dof_handler().end();

          for (; cell!=endc; ++cell)
            if (cell->is_locally_owned())
              {
                const types::LevelInd found_cell(cell->level(),cell->index());
                const unsigned int n_particles_in_cell = particle_handler->n_particles_in_cell(cell);

                // Add particles if necessary
                if ((particle_load_balancing & ParticleLoadBalancing::add_particles) &&
                    (n_particles_in_cell < min_particles_per_cell))
                  {
                    for (unsigned int i = n_particles_in_cell; i < min_particles_per_cell; ++i,++local_next_particle_index)
                      {
                        std::pair<aspect::Particle::types::LevelInd,Particle<dim> > new_particle = generator->generate_particle(cell,local_next_particle_index);

                        const std::vector<double> particle_properties =
                          property_manager->initialize_late_particle(new_particle.second.get_location(),
                                                                     *particle_handler,
                                                                     *interpolator,
                                                                     cell);

                        typename ParticleHandler<dim>::particle_iterator particle = particle_handler->insert_particle(new_particle.second,
                                                                                    typename parallel::distributed::Triangulation<dim>::cell_iterator (&this->get_triangulation(),
                                                                                        new_particle.first.first,
                                                                                        new_particle.first.second));
                        particle->set_properties(particle_properties);
                      }
                  }

                // Remove particles if necessary
                else if ((particle_load_balancing & ParticleLoadBalancing::remove_particles) &&
                         (n_particles_in_cell > max_particles_per_cell))
                  {
                    const boost::iterator_range<typename ParticleHandler<dim>::particle_iterator> particles_in_cell
                      = particle_handler->particles_in_cell(cell);

                    const unsigned int n_particles_to_remove = n_particles_in_cell - max_particles_per_cell;

                    std::set<unsigned int> particle_ids_to_remove;
                    while (particle_ids_to_remove.size() < n_particles_to_remove)
                      particle_ids_to_remove.insert(random_number_generator() % n_particles_in_cell);

                    std::list<typename ParticleHandler<dim>::particle_iterator> particles_to_remove;

                    for (std::set<unsigned int>::const_iterator id = particle_ids_to_remove.begin();
                         id != particle_ids_to_remove.end(); ++id)
                      {
                        typename ParticleHandler<dim>::particle_iterator particle_to_remove = particles_in_cell.begin();
                        std::advance(particle_to_remove,*id);

                        particles_to_remove.push_back(particle_to_remove);
                      }

                    for (typename std::list<typename ParticleHandler<dim>::particle_iterator>::iterator particle = particles_to_remove.begin();
                         particle != particles_to_remove.end(); ++particle)
                      {
                        particle_handler->remove_particle(*particle);
                      }
                  }
              }

          particle_handler->update_n_global_particles();
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
          const unsigned int n_particles_in_cell = particle_handler->n_particles_in_cell(cell);
          return n_particles_in_cell * particle_weight;
        }
      else if (status == parallel::distributed::Triangulation<dim>::CELL_COARSEN)
        {
          unsigned int n_particles_in_cell = 0;

          for (unsigned int child_index = 0; child_index < GeometryInfo<dim>::max_children_per_cell; ++child_index)
            n_particles_in_cell += particle_handler->n_particles_in_cell(cell->child(child_index));

          return n_particles_in_cell * particle_weight;
        }

      Assert (false, ExcInternalError());
      return 0;
    }


    template <int dim>
    std::map<types::subdomain_id, unsigned int>
    World<dim>::get_subdomain_id_to_neighbor_map() const
    {
      std::map<types::subdomain_id, unsigned int> subdomain_id_to_neighbor_map;
      const std::set<types::subdomain_id> ghost_owners = this->get_triangulation().ghost_owners();
      std::set<types::subdomain_id>::const_iterator ghost_owner = ghost_owners.begin();

      for (unsigned int neighbor_id=0; neighbor_id<ghost_owners.size(); ++neighbor_id,++ghost_owner)
        {
          subdomain_id_to_neighbor_map.insert(std::make_pair(*ghost_owner,neighbor_id));
        }
      return subdomain_id_to_neighbor_map;
    }

    template <int dim>
    void
    World<dim>::move_particles_back_into_mesh()
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

          const std::map<types::subdomain_id, unsigned int> subdomain_to_neighbor_map(get_subdomain_id_to_neighbor_map());

          std::vector<bool> periodic(dim,false);
          std::set< std::pair< std::pair<types::boundary_id, types::boundary_id>, unsigned int> >::const_iterator boundary =
            periodic_boundaries.begin();
          for (; boundary != periodic_boundaries.end(); ++boundary)
            periodic[boundary->second] = true;

          typename ParticleHandler<dim>::particle_iterator particle = particle_handler->begin();
          for (; particle != particle_handler->end(); ++particle)
            {
              // modify the particle position if it crossed a periodic boundary
              Point<dim> particle_position = particle->get_location();
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
              particle->set_location(particle_position);
            }
        }
    }

    template <int dim>
    void
    World<dim>::local_initialize_particles(const typename ParticleHandler<dim>::particle_iterator &begin_particle,
                                           const typename ParticleHandler<dim>::particle_iterator &end_particle)
    {
      for (typename ParticleHandler<dim>::particle_iterator it = begin_particle; it!=end_particle; ++it)
        property_manager->initialize_one_particle(it);
    }

    template <int dim>
    void
    World<dim>::local_update_particles(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                       const typename ParticleHandler<dim>::particle_iterator &begin_particle,
                                       const typename ParticleHandler<dim>::particle_iterator &end_particle)
    {
      const unsigned int particles_in_cell = std::distance(begin_particle,end_particle);
      const unsigned int solution_components = this->introspection().n_components;

      Vector<double>              value (solution_components);
      std::vector<Tensor<1,dim> > gradient (solution_components,Tensor<1,dim>());

      std::vector<Vector<double> >              values(particles_in_cell,value);
      std::vector<std::vector<Tensor<1,dim> > > gradients(particles_in_cell,gradient);
      std::vector<Point<dim> >                  positions(particles_in_cell);

      typename ParticleHandler<dim>::particle_iterator it = begin_particle;
      for (unsigned int i = 0; it!=end_particle; ++it,++i)
        {
          positions[i] = it->get_reference_location();
        }

      const Quadrature<dim> quadrature_formula(positions);
      const UpdateFlags update_flags = property_manager->get_needed_update_flags();
      FEValues<dim> fe_value (this->get_mapping(),
                              this->get_fe(),
                              quadrature_formula,
                              update_flags);

      fe_value.reinit (cell);
      if (update_flags & update_values)
        fe_value.get_function_values (this->get_solution(),
                                      values);
      if (update_flags & update_gradients)
        fe_value.get_function_gradients (this->get_solution(),
                                         gradients);

      it = begin_particle;
      for (unsigned int i = 0; it!=end_particle; ++it,++i)
        {
          property_manager->update_one_particle(it,
                                                values[i],
                                                gradients[i]);
        }
    }

    template <int dim>
    void
    World<dim>::local_advect_particles(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                       const typename ParticleHandler<dim>::particle_iterator &begin_particle,
                                       const typename ParticleHandler<dim>::particle_iterator &end_particle)
    {
      const unsigned int particles_in_cell = std::distance(begin_particle,end_particle);

      std::vector<Tensor<1,dim> >  velocity(particles_in_cell);
      std::vector<Tensor<1,dim> >  old_velocity(particles_in_cell);

      // Below we manually evaluate the solution at all support points of the
      // current cell, and then use the shape functions to interpolate the
      // solution to the particle points. All of this can be done with less
      // code using an FEValues object, but since this object initializes a lot
      // of memory for other purposes and we can not reuse the FEValues object
      // for other cells, it is much faster to do the work manually. Also this
      // function is quite performance critical.

      std::vector<types::global_dof_index> cell_dof_indices (this->get_fe().dofs_per_cell);
      cell->get_dof_indices (cell_dof_indices);

      const FiniteElement<dim> &velocity_fe = this->get_fe().base_element(this->introspection()
                                                                          .base_elements.velocities);

      const bool compute_fluid_velocity = this->include_melt_transport() &&
                                          property_manager->get_data_info().fieldname_exists("melt_presence");

      const unsigned int fluid_component_index = (compute_fluid_velocity ?
                                                  this->introspection().variable("fluid velocity").first_component_index
                                                  :
                                                  numbers::invalid_unsigned_int);

      std::vector<bool> use_fluid_velocity((compute_fluid_velocity ?
                                            particles_in_cell
                                            :
                                            0), false);

      if (compute_fluid_velocity)
        {
          const unsigned int melt_property_index = property_manager->get_data_info()
                                                   .get_position_by_field_name("melt_presence");

          typename ParticleHandler<dim>::particle_iterator it = begin_particle;
          for (unsigned int particle_index = 0; it!=end_particle; ++it,++particle_index)
            if (it->get_properties()[melt_property_index] == 1.0)
              use_fluid_velocity[particle_index] = true;
        }

      for (unsigned int j=0; j<velocity_fe.dofs_per_cell; ++j)
        {
          Tensor<1,dim> velocity_at_support_point;
          Tensor<1,dim> old_velocity_at_support_point;

          for (unsigned int dir=0; dir<dim; ++dir)
            {
              const unsigned int support_point_index
                = this->get_fe().component_to_system_index(this->introspection()
                                                           .component_indices.velocities[dir],j);

              velocity_at_support_point[dir] = this->get_solution()[cell_dof_indices[support_point_index]];
              old_velocity_at_support_point[dir] = this->get_old_solution()[cell_dof_indices[support_point_index]];
            }

          Tensor<1,dim> fluid_velocity_at_support_point;
          Tensor<1,dim> old_fluid_velocity_at_support_point;

          if (compute_fluid_velocity)
            for (unsigned int dir=0; dir<dim; ++dir)
              {
                const unsigned int support_point_index
                  = this->get_fe().component_to_system_index(fluid_component_index + dir,j);

                fluid_velocity_at_support_point[dir] = this->get_solution()[cell_dof_indices[support_point_index]];
                old_fluid_velocity_at_support_point[dir] = this->get_old_solution()[cell_dof_indices[support_point_index]];
              }

          typename ParticleHandler<dim>::particle_iterator it = begin_particle;
          for (unsigned int particle_index = 0; it!=end_particle; ++it,++particle_index)
            {
              // melt FE uses the same FE so the shape value is the same
              const double shape_value = velocity_fe.shape_value(j,it->get_reference_location());

              if (compute_fluid_velocity && use_fluid_velocity[particle_index])
                {
                  velocity[particle_index] += fluid_velocity_at_support_point * shape_value;
                  old_velocity[particle_index] += old_fluid_velocity_at_support_point * shape_value;
                }
              else
                {
                  velocity[particle_index] += velocity_at_support_point * shape_value;
                  old_velocity[particle_index] += old_velocity_at_support_point * shape_value;
                }
            }
        }

      integrator->local_integrate_step(begin_particle,
                                       end_particle,
                                       old_velocity,
                                       velocity,
                                       this->get_timestep());
    }

    template <int dim>
    void
    World<dim>::setup_initial_state ()
    {
      // If we are in the first adaptive refinement cycle generate particles
      if (this->get_pre_refinement_step() == 0)
        generate_particles();

      // And initialize the particle properties according to the initial
      // conditions on the current mesh
      initialize_particles();
    }

    template <int dim>
    void
    World<dim>::generate_particles()
    {
      TimerOutput::Scope timer_section(this->get_computing_timer(), "Particles: Generate");

      std::multimap<types::LevelInd, Particle<dim> > particles;
      generator->generate_particles(particles);
      particle_handler->insert_particles(particles);

      particle_handler->update_n_global_particles();
      particle_handler->update_next_free_particle_index();
    }

    template <int dim>
    void
    World<dim>::initialize_particles()
    {
      // Initialize the particle's access to the property_pool. This is necessary
      // even if the Particle do not carry properties, because they need a
      // way to determine the number of properties they carry.
      for (ParticleIterator<dim> particle = particle_handler->begin(); particle!=particle_handler->end(); ++particle)
        particle->set_property_pool(particle_handler->get_property_pool());


      // TODO: Change this loop over all cells to use the WorkStream interface
      if (property_manager->get_n_property_components() > 0)
        {
          TimerOutput::Scope timer_section(this->get_computing_timer(), "Particles: Initialize properties");

          particle_handler->get_property_pool().reserve(2 * particle_handler->n_locally_owned_particles());

          // Loop over all cells and initialize the particles cell-wise
          typename DoFHandler<dim>::active_cell_iterator
          cell = this->get_dof_handler().begin_active(),
          endc = this->get_dof_handler().end();

          for (; cell!=endc; ++cell)
            if (cell->is_locally_owned())
              {
                typename ParticleHandler<dim>::particle_iterator_range
                particles_in_cell = particle_handler->particles_in_cell(cell);

                // Only initialize particles, if there are any in this cell
                if (particles_in_cell.begin() != particles_in_cell.end())
                  local_initialize_particles(particles_in_cell.begin(),
                                             particles_in_cell.end());
              }
          if (update_ghost_particles &&
              dealii::Utilities::MPI::n_mpi_processes(this->get_mpi_communicator()) > 1)
            {
              TimerOutput::Scope timer_section(this->get_computing_timer(), "Particles: Exchange ghosts");
              particle_handler->exchange_ghost_particles();
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
          TimerOutput::Scope timer_section(this->get_computing_timer(), "Particles: Update properties");

          // Loop over all cells and update the particles cell-wise
          typename DoFHandler<dim>::active_cell_iterator
          cell = this->get_dof_handler().begin_active(),
          endc = this->get_dof_handler().end();

          for (; cell!=endc; ++cell)
            if (cell->is_locally_owned())
              {
                typename ParticleHandler<dim>::particle_iterator_range
                particles_in_cell = particle_handler->particles_in_cell(cell);

                // Only update particles, if there are any in this cell
                if (particles_in_cell.begin() != particles_in_cell.end())
                  local_update_particles(cell,
                                         particles_in_cell.begin(),
                                         particles_in_cell.end());
              }
        }
    }

    template <int dim>
    void
    World<dim>::advect_particles()
    {
      {
        // TODO: Change this loop over all cells to use the WorkStream interface
        TimerOutput::Scope timer_section(this->get_computing_timer(), "Particles: Advect");

        // Loop over all cells and advect the particles cell-wise
        typename DoFHandler<dim>::active_cell_iterator
        cell = this->get_dof_handler().begin_active(),
        endc = this->get_dof_handler().end();

        for (; cell!=endc; ++cell)
          if (cell->is_locally_owned())
            {
              const typename ParticleHandler<dim>::particle_iterator_range
              particles_in_cell = particle_handler->particles_in_cell(cell);

              // Only advect particles, if there are any in this cell
              if (particles_in_cell.begin() != particles_in_cell.end())
                local_advect_particles(cell,
                                       particles_in_cell.begin(),
                                       particles_in_cell.end());
            }

        // If particles fell out of the mesh, put them back in if they have crossed
        // a periodic boundary. If they have left the mesh otherwise, they will be
        // discarded during the next call to
        // particle_handler->sort_particles_into_subdomains_and_cells()
        move_particles_back_into_mesh();
      }

      {
        TimerOutput::Scope timer_section(this->get_computing_timer(), "Particles: Sort");
        // Find the cells that the particles moved to
        particle_handler->sort_particles_into_subdomains_and_cells();
      }
    }

    template <int dim>
    void
    World<dim>::advance_timestep()
    {
      do
        {
          advect_particles();
        }
      // Keep calling the integrator until it indicates it is finished
      while (integrator->new_integration_step());

      apply_particle_per_cell_bounds();

      // Update particle properties
      if (property_manager->need_update() == Property::update_time_step)
        update_particles();

      // Update the number of global particles if some have left the domain
      particle_handler->update_n_global_particles();

      // Now that all particle information was updated, exchange the new
      // ghost particles.
      if (update_ghost_particles &&
          dealii::Utilities::MPI::n_mpi_processes(this->get_mpi_communicator()) > 1)
        {
          TimerOutput::Scope timer_section(this->get_computing_timer(), "Particles: Exchange ghosts");
          particle_handler->exchange_ghost_particles();
        }
    }

    template <int dim>
    void
    World<dim>::save (std::ostringstream &os) const
    {
      aspect::oarchive oa (os);
      oa << (*this);
#if !DEAL_II_VERSION_GTE(9,0,0)
      output->save(os);
#endif
    }

    template <int dim>
    void
    World<dim>::load (std::istringstream &is)
    {
      aspect::iarchive ia (is);
      ia >> (*this);
#if !DEAL_II_VERSION_GTE(9,0,0)
      output->load(is);
#endif
    }

    template <int dim>
    void
    World<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Particles");
        {
          prm.declare_entry ("Load balancing strategy", "repartition",
                             Patterns::MultipleSelection ("none|remove particles|add particles|"
                                                          "remove and add particles|repartition"),
                             "Strategy that is used to balance the computational "
                             "load across processors for adaptive meshes.");
          prm.declare_entry ("Minimum particles per cell", "0",
                             Patterns::Integer (0),
                             "Lower limit for particle number per cell. This limit is "
                             "useful for adaptive meshes to prevent fine cells from being empty "
                             "of particles. It will be checked and enforced after mesh "
                             "refinement and after particle movement. "
                             "If there are "
                             "\\texttt{n\\_number\\_of\\_particles} $<$ \\texttt{min\\_particles\\_per\\_cell} "
                             "particles in one cell then "
                             "\\texttt{min\\_particles\\_per\\_cell} - \\texttt{n\\_number\\_of\\_particles} "
                             "particles are generated and randomly placed in "
                             "this cell. If the particles carry properties the "
                             "individual property plugins control how the "
                             "properties of the new particles are initialized.");
          prm.declare_entry ("Maximum particles per cell", "100",
                             Patterns::Integer (0),
                             "Upper limit for particle number per cell. This limit is "
                             "useful for adaptive meshes to prevent coarse cells from slowing down "
                             "the whole model. It will be checked and enforced after mesh "
                             "refinement, after MPI transfer of particles and after particle "
                             "movement. If there are "
                             "\\texttt{n\\_number\\_of\\_particles} $>$ \\texttt{max\\_particles\\_per\\_cell} "
                             "particles in one cell then "
                             "\\texttt{n\\_number\\_of\\_particles} - \\texttt{max\\_particles\\_per\\_cell} "
                             "particles in this cell are randomly chosen and destroyed.");
          prm.declare_entry ("Particle weight", "10",
                             Patterns::Integer (0),
                             "Weight that is associated with the computational load of "
                             "a single particle. The sum of particle weights will be added "
                             "to the sum of cell weights to determine the partitioning of "
                             "the mesh if the `repartition' particle load balancing strategy "
                             "is selected. The optimal weight depends on the used "
                             "integrator and particle properties. In general for a more "
                             "expensive integrator and more expensive properties a larger "
                             "particle weight is recommended. Before adding the weights "
                             "of particles, each cell already carries a weight of 1000 to "
                             "account for the cost of field-based computations.");
          prm.declare_entry ("Update ghost particles", "false",
                             Patterns::Bool (),
                             "Some particle interpolation algorithms require knowledge "
                             "about particles in neighboring cells. To allow this, "
                             "particles in ghost cells need to be exchanged between the "
                             "processes neighboring this cell. This parameter determines "
                             "whether this transport is happening.");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();

      Generator::declare_parameters<dim>(prm);

      // Output of particle related data has been moved to deal.II from version 9.0 on
      // The relevant code now lives in postprocess/particles.cc
#if !DEAL_II_VERSION_GTE(9,0,0)
      Output::declare_parameters<dim>(prm);
#endif
      Integrator::declare_parameters<dim>(prm);
      Interpolator::declare_parameters<dim>(prm);
      Property::Manager<dim>::declare_parameters(prm);
    }


    template <int dim>
    void
    World<dim>::parse_parameters (ParameterHandler &prm)
    {
      // First do some error checking. The current algorithm does not find
      // the cells around particles, if the particles moved more than one
      // cell in one timestep and we are running in parallel, because they
      // skip the layer of ghost cells around our local domain. Assert this
      // is not possible.
      const double CFL_number = prm.get_double ("CFL number");
      const unsigned int n_processes = Utilities::MPI::n_mpi_processes(this->get_mpi_communicator());

      AssertThrow((n_processes == 1) || (CFL_number <= 1.0),
                  ExcMessage("The current particle algorithm does not work in "
                             "parallel if the CFL number is larger than 1.0, because "
                             "in this case particles can move more than one cell "
                             "diameter in one time step and therefore skip the layer "
                             "of ghost cells around the local subdomain."));

      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Particles");
        {
          min_particles_per_cell = prm.get_integer("Minimum particles per cell");
          max_particles_per_cell = prm.get_integer("Maximum particles per cell");

          AssertThrow(min_particles_per_cell <= max_particles_per_cell,
                      ExcMessage("Please select a 'Minimum particles per cell' parameter "
                                 "that is smaller than or equal to the 'Maximum particles per cell' parameter."));

          particle_weight = prm.get_integer("Particle weight");

          update_ghost_particles = prm.get_bool("Update ghost particles");

          const std::vector<std::string> strategies = Utilities::split_string_list(prm.get ("Load balancing strategy"));
          AssertThrow(Utilities::has_unique_entries(strategies),
                      ExcMessage("The list of strings for the parameter "
                                 "'Postprocess/Particles/Load balancing strategy' contains entries more than once. "
                                 "This is not allowed. Please check your parameter file."));

          particle_load_balancing = ParticleLoadBalancing::no_balancing;

          for (std::vector<std::string>::const_iterator strategy = strategies.begin(); strategy != strategies.end(); ++strategy)
            {
              if (*strategy == "remove particles")
                particle_load_balancing = typename ParticleLoadBalancing::Kind(particle_load_balancing | ParticleLoadBalancing::remove_particles);
              else if (*strategy == "add particles")
                particle_load_balancing = typename ParticleLoadBalancing::Kind(particle_load_balancing | ParticleLoadBalancing::add_particles);
              else if (*strategy == "remove and add particles")
                particle_load_balancing = typename ParticleLoadBalancing::Kind(particle_load_balancing | ParticleLoadBalancing::remove_and_add_particles);
              else if (*strategy == "repartition")
                particle_load_balancing = typename ParticleLoadBalancing::Kind(particle_load_balancing | ParticleLoadBalancing::repartition);
              else if (*strategy == "none")
                {
                  particle_load_balancing = ParticleLoadBalancing::no_balancing;
                  AssertThrow(strategies.size() == 1,
                              ExcMessage("The particle load balancing strategy `none' is not compatible "
                                         "with any other strategy, yet it seems another is selected as well. "
                                         "Please check the parameter file."));
                }
              else
                AssertThrow(false,
                            ExcMessage("The 'Load balancing strategy' parameter contains an unknown value: <" + *strategy
                                       + ">. This value does not correspond to any known load balancing strategy. Possible values "
                                       "are listed in the corresponding manual subsection."));
            }

        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();

      TimerOutput::Scope timer_section(this->get_computing_timer(), "Particles: Initialization");

      // Create a generator object depending on what the parameters specify
      generator.reset(Generator::create_particle_generator<dim> (prm));
      if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(generator.get()))
        sim->initialize_simulator (this->get_simulator());
      generator->parse_parameters(prm);
      generator->initialize();

      // Output of particle related data has been moved to deal.II from version 9.0 on
      // The relevant code now lives in postprocess/particles.cc
#if !DEAL_II_VERSION_GTE(9,0,0)
      // Create an output object depending on what the parameters specify
      output.reset(Output::create_particle_output<dim>
                   (prm));

      // We allow to not generate any output plugin, in which case output is
      // a null pointer. Only initialize output if it was created.
      if (output)
        {
          if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(output.get()))
            sim->initialize_simulator (this->get_simulator());
          output->parse_parameters(prm);
          output->initialize();
        }
#endif

      // Create an integrator object depending on the specified parameter
      integrator.reset(Integrator::create_particle_integrator<dim> (prm));
      if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(integrator.get()))
        sim->initialize_simulator (this->get_simulator());
      integrator->parse_parameters(prm);

      // Create an property_manager object and initialize its properties
      property_manager.reset(new Property::Manager<dim> ());
      SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(property_manager.get());
      sim->initialize_simulator (this->get_simulator());
      property_manager->parse_parameters(prm);
      property_manager->initialize();

      // Create an interpolator object depending on the specified parameter
      interpolator.reset(Interpolator::create_particle_interpolator<dim> (prm));
      if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(interpolator.get()))
        sim->initialize_simulator (this->get_simulator());
      interpolator->parse_parameters(prm);
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
