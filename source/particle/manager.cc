/*
  Copyright (C) 2015 - 2024 by the authors of the ASPECT code.

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

#include <aspect/particle/manager.h>
#include <aspect/global.h>
#include <aspect/utilities.h>
#include <aspect/simulator.h>
#include <aspect/melt.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/fe/mapping_cartesian.h>

#include <boost/serialization/map.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

namespace aspect
{
  namespace Particle
  {
    template <int dim>
    Manager<dim>::Manager()
      = default;

    template <int dim>
    Manager<dim>::~Manager()
      = default;

    template <int dim>
    Manager<dim>::Manager(Manager &&other) noexcept
  :
    generator(std::move(other.generator)),
              integrator(std::move(other.integrator)),
              interpolator(std::move(other.interpolator)),
              particle_handler(std::move(other.particle_handler)),
              particle_handler_backup(), // can not move
              property_manager(std::move(other.property_manager)),
              particle_load_balancing(other.particle_load_balancing),
              min_particles_per_cell(other.min_particles_per_cell),
              max_particles_per_cell(other.max_particles_per_cell),
              particle_weight(other.particle_weight)
    {}



    template <int dim>
    void
    Manager<dim>::initialize()
    {
      CitationInfo::add("particles");

      // Create a particle handler that stores the future particles.
      // If we restarted from a checkpoint we will fill this particle handler
      // later with its serialized variables and stored particles
      particle_handler = std::make_unique<ParticleHandler<dim>>(this->get_triangulation(),
                                                                 this->get_mapping(),
                                                                 property_manager->get_n_property_components());

      particle_handler_backup.initialize(this->get_triangulation(),
                                         this->get_mapping(),
                                         property_manager->get_n_property_components());

      connect_to_signals(this->get_signals());
    }



    template <int dim>
    void
    Manager<dim>::update()
    {
      generator->update();
      integrator->update();
      interpolator->update();
      property_manager->update();
    }



    template <int dim>
    const Property::Manager<dim> &
    Manager<dim>::get_property_manager() const
    {
      return *property_manager;
    }



    template <int dim>
    const Particles::ParticleHandler<dim> &
    Manager<dim>::get_particle_handler() const
    {
      return *particle_handler.get();
    }



    template <int dim>
    Particles::ParticleHandler<dim> &
    Manager<dim>::get_particle_handler()
    {
      return *particle_handler.get();
    }



    template <int dim>
    void
    Manager<dim>::copy_particle_handler (const Particles::ParticleHandler<dim> &from_particle_handler,
                                         Particles::ParticleHandler<dim> &to_particle_handler) const
    {
      {
        TimerOutput::Scope timer_section(this->get_computing_timer(), "Particles: Copy");

        to_particle_handler.copy_from(from_particle_handler);
      }
    }



    template <int dim>
    void
    Manager<dim>::backup_particles ()
    {
      copy_particle_handler (*particle_handler.get(), particle_handler_backup);
    }



    template <int dim>
    void
    Manager<dim>::restore_particles ()
    {
      copy_particle_handler (particle_handler_backup, *particle_handler.get());
    }



    template <int dim>
    const Interpolator::Interface<dim> &
    Manager<dim>::get_interpolator() const
    {
      return *interpolator;
    }



    template <int dim>
    types::particle_index
    Manager<dim>::n_global_particles() const
    {
      return particle_handler->n_global_particles();
    }



    template <int dim>
    void
    Manager<dim>::connect_to_signals(aspect::SimulatorSignals<dim> &signals)
    {
      signals.post_set_initial_state.connect(
        [&] (const SimulatorAccess<dim> &)
      {
        this->setup_initial_state();
      });

      connect_particle_handler_signals(signals,*particle_handler);
      // Particle handler backup will not be stored for checkpointing
      connect_particle_handler_signals(signals, particle_handler_backup, false);

      signals.post_refinement_load_user_data.connect(
        [&] (typename parallel::distributed::Triangulation<dim> &)
      {
        this->apply_particle_per_cell_bounds();
      });

      signals.post_resume_load_user_data.connect(
        [&] (typename parallel::distributed::Triangulation<dim> &)
      {
        this->apply_particle_per_cell_bounds();
      });
    }



    template <int dim>
    void
    Manager<dim>::connect_particle_handler_signals(aspect::SimulatorSignals<dim> &signals,
                                                   ParticleHandler<dim> &particle_handler_,
                                                   const bool connect_to_checkpoint_signals) const
    {
      signals.pre_refinement_store_user_data.connect(
        [&] (typename parallel::distributed::Triangulation<dim> &)
      {
        particle_handler_.prepare_for_coarsening_and_refinement();
      });

      signals.post_refinement_load_user_data.connect(
        [&] (typename parallel::distributed::Triangulation<dim> &)
      {
        particle_handler_.unpack_after_coarsening_and_refinement();
      });

      // Only connect to checkpoint signals if requested
      if (connect_to_checkpoint_signals)
        {
          signals.pre_checkpoint_store_user_data.connect(
            [&] (typename parallel::distributed::Triangulation<dim> &)
          {
            particle_handler_.prepare_for_serialization();
          });

          signals.post_resume_load_user_data.connect(
            [&] (typename parallel::distributed::Triangulation<dim> &)
          {
            particle_handler_.deserialize();
          });
        }

      if (dealii::Utilities::MPI::n_mpi_processes(this->get_mpi_communicator()) > 1)
        {
          auto do_ghost_exchange = [&] (typename parallel::distributed::Triangulation<dim> &)
          {
            particle_handler_.exchange_ghost_particles();
          };
          signals.post_refinement_load_user_data.connect(do_ghost_exchange);
          signals.post_resume_load_user_data.connect(do_ghost_exchange);
        }

      signals.post_mesh_deformation.connect(
        [&] (const SimulatorAccess<dim> &)
      {
        particle_handler->sort_particles_into_subdomains_and_cells();
      },
      boost::signals2::at_front);
    }



    template <int dim>
    void
    Manager<dim>::apply_particle_per_cell_bounds()
    {
      // If any load balancing technique is selected that creates/destroys particles
      if (particle_load_balancing & ParticleLoadBalancing::remove_and_add_particles)
        {
          // First do some preparation for particle generation in poorly
          // populated areas. For this we need to know which particle ids to
          // generate so that they are globally unique.
          // Ensure this by communicating the number of particles that every
          // process is going to generate.
          particle_handler->update_cached_numbers();
          types::particle_index local_next_particle_index = particle_handler->get_next_free_particle_index();
          if (particle_load_balancing & ParticleLoadBalancing::add_particles)
            {
              types::particle_index particles_to_add_locally = 0;

              // Loop over all cells and determine the number of particles to generate
              for (const auto &cell : this->get_dof_handler().active_cell_iterators())
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

              const int ierr = MPI_Scan(&particles_to_add_locally, &local_start_index, 1, DEAL_II_PARTICLE_INDEX_MPI_TYPE, MPI_SUM, this->get_mpi_communicator());
              AssertThrowMPI(ierr);

              local_start_index -= particles_to_add_locally;
              local_next_particle_index += local_start_index;

              const types::particle_index globally_generated_particles =
                dealii::Utilities::MPI::sum(particles_to_add_locally,this->get_mpi_communicator());

              AssertThrow (particle_handler->get_next_free_particle_index()
                           <= std::numeric_limits<types::particle_index>::max() - globally_generated_particles,
                           ExcMessage("There is no free particle index left to generate a new particle id. Please check if your "
                                      "model generates unusually many new particles (by repeatedly deleting and regenerating particles), or "
                                      "recompile deal.II with the DEAL_II_WITH_64BIT_INDICES option enabled, to use 64-bit integers for "
                                      "particle ids."));
            }

          std::mt19937 random_number_generator;

          // Loop over all cells and generate or remove the particles cell-wise
          for (const auto &cell : this->get_dof_handler().active_cell_iterators())
            if (cell->is_locally_owned())
              {
                const unsigned int n_particles_in_cell = particle_handler->n_particles_in_cell(cell);

                // Add particles if necessary
                if ((particle_load_balancing & ParticleLoadBalancing::add_particles) &&
                    (n_particles_in_cell < min_particles_per_cell))
                  {
                    for (unsigned int i = n_particles_in_cell; i < min_particles_per_cell; ++i,++local_next_particle_index)
                      {
                        std::pair<Particles::internal::LevelInd,Particles::Particle<dim>> new_particle = generator->generate_particle(cell,local_next_particle_index);

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
                    const unsigned int n_particles_to_remove = n_particles_in_cell - max_particles_per_cell;
                    for (unsigned int i=0; i < n_particles_to_remove; ++i)
                      {
                        const unsigned int current_n_particles_in_cell = particle_handler->n_particles_in_cell(cell);
                        const unsigned int index_to_remove = std::uniform_int_distribution<unsigned int>
                                                             (0,current_n_particles_in_cell-1)(random_number_generator);

                        auto particle_to_remove = particle_handler->particles_in_cell(cell).begin();
                        std::advance(particle_to_remove, index_to_remove);
                        particle_handler->remove_particle(particle_to_remove);
                      }

                  }
              }

          particle_handler->update_cached_numbers();
        }
    }

    template <int dim>
    unsigned int
    Manager<dim>::cell_weight(const typename parallel::distributed::Triangulation<dim>::cell_iterator &cell,
#if DEAL_II_VERSION_GTE(9,6,0)
                              const CellStatus status
#else
                              const typename parallel::distributed::Triangulation<dim>::CellStatus status
#endif
                             )
    {
      if (cell->is_active() && !cell->is_locally_owned())
        return 0;

      unsigned int n_particles_in_cell = 0;
      switch (status)
        {
#if DEAL_II_VERSION_GTE(9,6,0)
          case CellStatus::cell_will_persist:
          case CellStatus::cell_will_be_refined:
            n_particles_in_cell = particle_handler->n_particles_in_cell(cell);
            break;

          case CellStatus::cell_invalid:
            break;

          case CellStatus::children_will_be_coarsened:
            for (const auto &child : cell->child_iterators())
              n_particles_in_cell += particle_handler->n_particles_in_cell(child);
            break;
#else
          case parallel::distributed::Triangulation<dim>::CELL_PERSIST:
          case parallel::distributed::Triangulation<dim>::CELL_REFINE:
            n_particles_in_cell = particle_handler->n_particles_in_cell(cell);
            break;

          case parallel::distributed::Triangulation<dim>::CELL_INVALID:
            break;

          case parallel::distributed::Triangulation<dim>::CELL_COARSEN:
            for (const auto &child : cell->child_iterators())
              n_particles_in_cell += particle_handler->n_particles_in_cell(child);
            break;
#endif
          default:
            Assert(false, ExcInternalError());
            break;
        }
      return n_particles_in_cell * particle_weight;
    }


    template <int dim>
    std::map<types::subdomain_id, unsigned int>
    Manager<dim>::get_subdomain_id_to_neighbor_map() const
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
    Manager<dim>::local_initialize_particles(const typename ParticleHandler<dim>::particle_iterator &begin_particle,
                                             const typename ParticleHandler<dim>::particle_iterator &end_particle)
    {
      for (typename ParticleHandler<dim>::particle_iterator it = begin_particle; it!=end_particle; ++it)
        property_manager->initialize_one_particle(it);
    }



    template <int dim>
    void
    Manager<dim>::local_update_particles(Property::ParticleUpdateInputs<dim> &inputs,
                                         small_vector<Point<dim>> &positions,
                                         const std::vector<EvaluationFlags::EvaluationFlags> &evaluation_flags,
                                         SolutionEvaluator<dim> &evaluator)
    {
      const unsigned int n_particles = particle_handler->n_particles_in_cell(inputs.current_cell);

      typename ParticleHandler<dim>::particle_iterator_range particles = particle_handler->particles_in_cell(inputs.current_cell);

      positions.resize(n_particles);
      unsigned int p = 0;
      for (const auto &particle : particles)
        {
          positions[p] = particle.get_reference_location();
          ++p;
        }

      small_vector<double> solution_values(this->get_fe().dofs_per_cell);

      inputs.current_cell->get_dof_values(this->get_solution(),
                                          solution_values.begin(),
                                          solution_values.end());

      EvaluationFlags::EvaluationFlags evaluation_flags_union = EvaluationFlags::nothing;
      for (const auto &flag : evaluation_flags)
        evaluation_flags_union |= flag;

      if (evaluation_flags_union & (EvaluationFlags::values | EvaluationFlags::gradients))
        {
          // Reinitialize and evaluate the requested solution values and gradients
          evaluator.reinit(inputs.current_cell,
          {positions.data(), positions.size()});

          evaluator.evaluate({solution_values.data(),solution_values.size()},
                             evaluation_flags);
        }

      if (evaluation_flags_union & EvaluationFlags::values)
        inputs.solution.resize(n_particles,small_vector<double,50>(evaluator.n_components(), numbers::signaling_nan<double>()));

      if (evaluation_flags_union & EvaluationFlags::gradients)
        inputs.gradients.resize(n_particles,small_vector<Tensor<1,dim>,50>(evaluator.n_components(), numbers::signaling_nan<Tensor<1,dim>>()));

      for (unsigned int i = 0; i<n_particles; ++i)
        {
          // Evaluate the solution, but only if it is requested in the update_flags
          if (evaluation_flags_union & EvaluationFlags::values)
            evaluator.get_solution(i, {&inputs.solution[i][0],inputs.solution[i].size()}, evaluation_flags);

          // Evaluate the gradients, but only if they are requested in the update_flags
          if (evaluation_flags_union & EvaluationFlags::gradients)
            evaluator.get_gradients(i, {&inputs.gradients[i][0],inputs.gradients[i].size()}, evaluation_flags);
        }

      property_manager->update_particles(inputs,particles);
    }



    template <int dim>
    void
    Manager<dim>::local_advect_particles(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                         const typename ParticleHandler<dim>::particle_iterator &begin_particle,
                                         const typename ParticleHandler<dim>::particle_iterator &end_particle,
                                         SolutionEvaluator<dim> &evaluator)
    {
      const unsigned int n_particles_in_cell = particle_handler->n_particles_in_cell(cell);

      small_vector<Point<dim>> positions;
      positions.reserve(n_particles_in_cell);
      for (auto particle = begin_particle; particle!=end_particle; ++particle)
        positions.push_back(particle->get_reference_location());

      const std::array<bool, 3> required_solution_vectors = integrator->required_solution_vectors();

      AssertThrow (required_solution_vectors[0] == false,
                   ExcMessage("The integrator requires the old old solution vector, but it is not available."));


      const bool use_fluid_velocity = this->include_melt_transport() &&
                                      property_manager->get_data_info().fieldname_exists("melt_presence");

      auto &velocity_evaluator = evaluator.get_velocity_or_fluid_velocity_evaluator(use_fluid_velocity);
      auto &mapping_info = evaluator.get_mapping_info();
      mapping_info.reinit(cell, {positions.data(),positions.size()});

      std::vector<Tensor<1,dim>> velocities;
      std::vector<Tensor<1,dim>> old_velocities;

      if (required_solution_vectors[1] == true)
        {
          small_vector<double> old_solution_values(this->get_fe().dofs_per_cell);
          cell->get_dof_values(this->get_old_solution(),
                               old_solution_values.begin(),
                               old_solution_values.end());

          velocity_evaluator.evaluate({old_solution_values.data(),old_solution_values.size()},
                                      EvaluationFlags::values);

          old_velocities.resize(n_particles_in_cell);
          for (unsigned int i=0; i<n_particles_in_cell; ++i)
            old_velocities[i] = velocity_evaluator.get_value(i);
        }

      if (required_solution_vectors[2] == true)
        {
          small_vector<double> solution_values(this->get_fe().dofs_per_cell);
          cell->get_dof_values(this->get_current_linearization_point(),
                               solution_values.begin(),
                               solution_values.end());
          velocity_evaluator.evaluate({solution_values.data(),solution_values.size()},
                                      EvaluationFlags::values);

          velocities.resize(n_particles_in_cell);
          for (unsigned int i=0; i<n_particles_in_cell; ++i)
            velocities[i] = velocity_evaluator.get_value(i);
        }

      integrator->local_integrate_step(begin_particle,
                                       end_particle,
                                       old_velocities,
                                       velocities,
                                       this->get_timestep());
    }



    template <int dim>
    void
    Manager<dim>::setup_initial_state ()
    {
      // We want to generate a new set of particles in each adaptive refinement
      // cycle to get the right number of particles per cell and to accurately
      // initialize their properties. Delete existing particles beforehand.
      if (this->get_pre_refinement_step() > 0)
        particle_handler->clear();

      // Generate particles in each adaptive refinement cycle
      generate_particles();

      // And initialize the particle properties according to the initial
      // conditions on the current mesh
      initialize_particles();
    }



    template <int dim>
    void
    Manager<dim>::generate_particles()
    {
      TimerOutput::Scope timer_section(this->get_computing_timer(), "Particles: Generate");
      generator->generate_particles(*particle_handler);
    }



    template <int dim>
    void
    Manager<dim>::initialize_particles()
    {
      // TODO: Change this loop over all cells to use the WorkStream interface
      if (property_manager->get_n_property_components() > 0)
        {
          TimerOutput::Scope timer_section(this->get_computing_timer(), "Particles: Initialize properties");

          particle_handler->get_property_pool().reserve(2 * particle_handler->n_locally_owned_particles());


          if (particle_handler->n_locally_owned_particles() > 0)
            local_initialize_particles(particle_handler->begin(),
                                       particle_handler->end());

          if (dealii::Utilities::MPI::n_mpi_processes(this->get_mpi_communicator()) > 1)
            {
              TimerOutput::Scope timer_section(this->get_computing_timer(), "Particles: Exchange ghosts");
              particle_handler->exchange_ghost_particles();
            }
        }
    }



    template <int dim>
    void
    Manager<dim>::update_particles()
    {
      // TODO: Change this loop over all cells to use the WorkStream interface

      if (property_manager->get_n_property_components() > 0)
        {
          TimerOutput::Scope timer_section(this->get_computing_timer(), "Particles: Update properties");

          Assert(dealii::internal::FEPointEvaluation::is_fast_path_supported(this->get_mapping()) == true,
                 ExcMessage("The particle system was optimized for deal.II mappings that support the fast evaluation path "
                            "of the class FEPointEvaluation. The mapping currently in use does not support this path. "
                            "It is safe to uncomment this assertion, but you can expect a performance penalty."));

          const std::vector<UpdateFlags> update_flags = property_manager->get_update_flags();

          // combine all update flags to a single flag, which is the required information
          // for the mapping inside the solution evaluator
          UpdateFlags mapping_flags = update_flags[0];
          for (unsigned int i=1; i<update_flags.size(); ++i)
            mapping_flags |= update_flags[i];

          std::unique_ptr<SolutionEvaluator<dim>> evaluator = construct_solution_evaluator(*this,
                                                               mapping_flags);

          // FEPointEvaluation uses different evaluation flags than the common UpdateFlags.
          // Translate between the two.
          std::vector<EvaluationFlags::EvaluationFlags> evaluation_flags (update_flags.size(), EvaluationFlags::nothing);

          for (unsigned int i=0; i<update_flags.size(); ++i)
            {
              if (update_flags[i] & update_values)
                evaluation_flags[i] |= EvaluationFlags::values;

              if (update_flags[i] & update_gradients)
                evaluation_flags[i] |= EvaluationFlags::gradients;
            }

          Property::ParticleUpdateInputs<dim> inputs;
          small_vector<Point<dim>> positions;

          // Loop over all cells and update the particles cell-wise
          for (const auto &cell : this->get_dof_handler().active_cell_iterators())
            if (cell->is_locally_owned())
              {
                // Only update particles if there are any in this cell
                if (particle_handler->n_particles_in_cell(cell) > 0)
                  {
                    inputs.current_cell = cell;
                    local_update_particles(inputs,
                                           positions,
                                           evaluation_flags,
                                           *evaluator);
                  }

              }
        }
    }



    template <int dim>
    void
    Manager<dim>::advect_particles()
    {
      {
        // TODO: Change this loop over all cells to use the WorkStream interface
        TimerOutput::Scope timer_section(this->get_computing_timer(), "Particles: Advect");

        Assert(dealii::internal::FEPointEvaluation::is_fast_path_supported(this->get_mapping()) == true,
               ExcMessage("The particle system was optimized for deal.II mappings that support the fast evaluation path "
                          "of the class FEPointEvaluation. The mapping currently in use does not support this path. "
                          "It is safe to uncomment this assertion, but you can expect a performance penalty."));

        std::unique_ptr<SolutionEvaluator<dim>> evaluator = construct_solution_evaluator(*this,
                                                             update_values);

        // Loop over all cells and advect the particles cell-wise
        for (const auto &cell : this->get_dof_handler().active_cell_iterators())
          if (cell->is_locally_owned())
            {
              const typename ParticleHandler<dim>::particle_iterator_range
              particles_in_cell = particle_handler->particles_in_cell(cell);

              // Only advect particles, if there are any in this cell
              if (particles_in_cell.begin() != particles_in_cell.end())
                {
                  local_advect_particles(cell,
                                         particles_in_cell.begin(),
                                         particles_in_cell.end(),
                                         *evaluator);
                }
            }
      }

      {
        TimerOutput::Scope timer_section(this->get_computing_timer(), "Particles: Sort");
        // Find the cells that the particles moved to
        particle_handler->sort_particles_into_subdomains_and_cells();
      }
    }



    template <int dim>
    void
    Manager<dim>::advance_timestep()
    {
      this->get_pcout() << "   Advecting particles... " << std::flush;
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

      // Now that all particle information was updated, exchange the new
      // ghost particles.
      if (dealii::Utilities::MPI::n_mpi_processes(this->get_mpi_communicator()) > 1)
        {
          TimerOutput::Scope timer_section(this->get_computing_timer(), "Particles: Exchange ghosts");
          particle_handler->exchange_ghost_particles();
        }
      this->get_pcout() << " done." << std::endl;
    }



    template <int dim>
    void
    Manager<dim>::save (std::ostringstream &os) const
    {
      aspect::oarchive oa (os);
      oa << (*this);
    }



    template <int dim>
    void
    Manager<dim>::load (std::istringstream &is)
    {
      aspect::iarchive ia (is);
      ia >> (*this);
    }



    template <int dim>
    void
    Manager<dim>::declare_parameters (ParameterHandler &prm)
    {
      constexpr unsigned int number_of_particle_managers = ASPECT_MAX_NUM_PARTICLE_SYSTEMS;
      for (unsigned int particle_manager = 0; particle_manager < number_of_particle_managers; ++particle_manager)
        {
          if (particle_manager == 0)
            {
              prm.enter_subsection("Particles");
            }
          else
            {
              prm.enter_subsection("Particles " + std::to_string(particle_manager+1));
            }
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
            prm.declare_entry ("Update ghost particles", "true",
                               Patterns::Bool (),
                               "Some particle interpolation algorithms require knowledge "
                               "about particles in neighboring cells. To allow this, "
                               "particles in ghost cells need to be exchanged between the "
                               "processes neighboring this cell. This parameter determines "
                               "whether this transport is happening. This parameter is "
                               "deprecated and will be removed in the future. Ghost particle "
                               "updates are always performed. Please set the parameter to `true'.");

            Generator::declare_parameters<dim>(prm);
            Integrator::declare_parameters<dim>(prm);
            Interpolator::declare_parameters<dim>(prm);

            Property::Manager<dim>::declare_parameters(prm);
          }
          prm.leave_subsection ();
        }

    }



    template <int dim>
    void
    Manager<dim>::parse_parameters (ParameterHandler &prm, const unsigned int particle_manager)
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

      if (particle_manager == 0)
        {
          prm.enter_subsection("Particles");
        }
      else
        {
          prm.enter_subsection("Particles " + std::to_string(particle_manager+1));
        }
      {
        min_particles_per_cell = prm.get_integer("Minimum particles per cell");
        max_particles_per_cell = prm.get_integer("Maximum particles per cell");

        AssertThrow(min_particles_per_cell <= max_particles_per_cell,
                    ExcMessage("Please select a 'Minimum particles per cell' parameter "
                               "that is smaller than or equal to the 'Maximum particles per cell' parameter."));

        particle_weight = prm.get_integer("Particle weight");

        const bool update_ghost_particles = prm.get_bool("Update ghost particles");
        AssertThrow(update_ghost_particles == true,
                    ExcMessage("The 'Update ghost particles' parameter is deprecated and will be removed in the future. "
                               "Ghost particle updates are always performed. Please set the parameter to `true'."));

        const std::vector<std::string> strategies = Utilities::split_string_list(prm.get ("Load balancing strategy"));
        AssertThrow(Utilities::has_unique_entries(strategies),
                    ExcMessage("The list of strings for the parameter "
                               "'Particles/Load balancing strategy' contains entries more than once. "
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

        if (particle_load_balancing & ParticleLoadBalancing::repartition)
          this->get_triangulation().signals.weight.connect(
#if DEAL_II_VERSION_GTE(9,6,0)
            [ &, particle_manager] (const typename parallel::distributed::Triangulation<dim>::cell_iterator &cell,
                                    const CellStatus status)
            -> unsigned int
#else
            [ &, particle_manager] (const typename parallel::distributed::Triangulation<dim>::cell_iterator &cell,
                                    const typename parallel::distributed::Triangulation<dim>::CellStatus status)
            -> unsigned int
#endif
          {
            // Only add the base weight of cells in particle manager 0, because all weights will be summed
            // across all particle managers.
            return (particle_manager == 0) ? 1000 + this->cell_weight(cell, status) : this->cell_weight(cell, status);
          });


        TimerOutput::Scope timer_section(this->get_computing_timer(), "Particles: Initialization");

        // Create a generator object depending on what the parameters specify
        generator = Generator::create_particle_generator<dim> (prm);
        if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(generator.get()))
          sim->initialize_simulator (this->get_simulator());
        generator->set_particle_manager_index(particle_manager);
        generator->parse_parameters(prm);
        generator->initialize();

        // Create a property_manager object and initialize its properties
        property_manager = std::make_unique<Property::Manager<dim>> ();
        SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(property_manager.get());
        sim->initialize_simulator (this->get_simulator());
        property_manager->set_particle_manager_index(particle_manager);
        property_manager->parse_parameters(prm);
        property_manager->initialize();

        // Create an integrator object depending on the specified parameter
        integrator = Integrator::create_particle_integrator<dim> (prm);
        if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(integrator.get()))
          sim->initialize_simulator (this->get_simulator());
        integrator->set_particle_manager_index(particle_manager);
        integrator->parse_parameters(prm);
        integrator->initialize();

        // Create an interpolator object depending on the specified parameter
        interpolator = Interpolator::create_particle_interpolator<dim> (prm);
        if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(interpolator.get()))
          sim->initialize_simulator (this->get_simulator());
        interpolator->set_particle_manager_index(particle_manager);
        interpolator->parse_parameters(prm);
        interpolator->initialize();

      }
      prm.leave_subsection ();
    }
  }
}


// explicit instantiation of the functions we implement in this file
namespace aspect
{
  namespace Particle
  {
#define INSTANTIATE(dim) \
  template class Manager<dim>;

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
