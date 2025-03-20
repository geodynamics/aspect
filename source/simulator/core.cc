/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

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


#include <aspect/simulator.h>
#include <aspect/global.h>
#include <aspect/utilities.h>
#include <aspect/melt.h>
#include <aspect/volume_of_fluid/handler.h>
#include <aspect/newton.h>
#include <aspect/stokes_matrix_free.h>
#include <aspect/mesh_deformation/interface.h>
#include <aspect/postprocess/particles.h>

#ifdef ASPECT_WITH_WORLD_BUILDER
#include <world_builder/world.h>
#endif

#include <aspect/simulator/assemblers/interface.h>
#include <aspect/geometry_model/initial_topography_model/zero_topography.h>
#include <aspect/material_model/rheology/elasticity.h>
#include <aspect/time_stepping/repeat_on_nonlinear_fail.h>

#include <deal.II/base/index_set.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/signaling_nan.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/mapping_cartesian.h>
#include <deal.II/fe/mapping_q_cache.h>

#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/derivative_approximation.h>
#include <deal.II/numerics/vector_tools.h>

#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/distributed/grid_refinement.h>

#include <fstream>
#include <iostream>
#include <locale>
#include <string>




namespace aspect
{
  namespace
  {
    /**
     * Helper function to construct the final std::vector of FEVariable before
     * it is used to construct the Introspection object. Create the default
     * setup based on parameters followed by a signal to allow modifications.
     */
    template <int dim>
    std::vector<VariableDeclaration<dim>> construct_variables(const Parameters<dim> &parameters,
                                                               SimulatorSignals<dim> &signals,
                                                               std::unique_ptr<MeltHandler<dim>> &melt_handler)
    {
      std::vector<VariableDeclaration<dim>> variables
        = construct_default_variables (parameters);
      if (melt_handler)
        melt_handler->edit_finite_element_variables (parameters, variables);

      signals.edit_finite_element_variables(variables);
      return variables;
    }

    /**
     * Helper function to construct mapping for the model.
     * The mapping is given by a degree four MappingQ for the case
     * of a curved mesh, and a cartesian mapping for a rectangular mesh that is
     * not deformed. Use a MappingQ1 if the initial mesh is deformed.
     * If mesh deformation is enabled, each mapping is later swapped out for a
     * MappingQ1Eulerian, which allows for mesh deformation during the
     * computation.
     */
    template <int dim>
    std::unique_ptr<Mapping<dim>>
    construct_mapping(const GeometryModel::Interface<dim> &geometry_model,
                      const InitialTopographyModel::Interface<dim> &initial_topography_model)
    {
      if (geometry_model.has_curved_elements())
        return std::make_unique<MappingQCache<dim>>(4);

      if (Plugins::plugin_type_matches<const InitialTopographyModel::ZeroTopography<dim>>(initial_topography_model))
        return std::make_unique<MappingCartesian<dim>>();

      return std::make_unique<MappingQ1<dim>>();
    }
  }


  /**
   * Constructor of the IntermediaryConstructorAction class. Since the
   * class has no members, there is nothing to initialize -- all we
   * need to do is execute the 'action' argument.
   */
  template <int dim>
  Simulator<dim>::IntermediaryConstructorAction::
  IntermediaryConstructorAction (const std::function<void ()> &action)
  {
    action();
  }



  /**
   * Constructor. Initialize all member variables.
   **/
  template <int dim>
  Simulator<dim>::Simulator (const MPI_Comm mpi_communicator_,
                             ParameterHandler &prm)
    :
    simulator_is_past_initialization (false),
    assemblers (std::make_unique<Assemblers::Manager<dim>>()),
    parameters (prm, mpi_communicator_),
    melt_handler (parameters.include_melt_transport ?
                  std::make_unique<MeltHandler<dim>>(prm) :
                  nullptr),
    newton_handler (Parameters<dim>::is_defect_correction(parameters.nonlinear_solver) ?
                    std::make_unique<NewtonHandler<dim>>() :
                    nullptr),
    post_signal_creation(
      std::bind (&internals::SimulatorSignals::call_connector_functions<dim>,
                 std::ref(signals))),
    volume_of_fluid_handler (parameters.volume_of_fluid_tracking_enabled ?
                             std::make_unique<VolumeOfFluidHandler<dim>> (*this, prm) :
                             nullptr),
    introspection (construct_variables<dim>(parameters, signals, melt_handler), parameters),
    mpi_communicator (Utilities::MPI::duplicate_communicator (mpi_communicator_)),
    iostream_tee_device(std::cout, log_file_stream),
    iostream_tee_stream(iostream_tee_device),
    pcout (iostream_tee_stream,
           (Utilities::MPI::
            this_mpi_process(mpi_communicator)
            == 0)),

    statistics_last_write_size (0),
    statistics_last_hash (0),

    computing_timer (mpi_communicator,
                     pcout,
                     TimerOutput::never,
                     TimerOutput::wall_times),
    total_walltime_until_last_snapshot(0.),
    initial_topography_model(InitialTopographyModel::create_initial_topography_model<dim>(prm)),
    geometry_model (GeometryModel::create_geometry_model<dim>(prm)),
    // make sure the parameters object gets a chance to
    // parse those parameters that depend on symbolic names
    // for boundary components
    post_geometry_model_creation_action (std::bind (&Parameters<dim>::parse_geometry_dependent_parameters,
                                                    std::ref(parameters),
                                                    std::ref(prm),
                                                    std::cref(*geometry_model))),
    material_model (MaterialModel::create_material_model<dim>(prm)),
    gravity_model (GravityModel::create_gravity_model<dim>(prm)),
    prescribed_stokes_solution (PrescribedStokesSolution::create_prescribed_stokes_solution<dim>(prm)),
    adiabatic_conditions (AdiabaticConditions::create_adiabatic_conditions<dim>(prm)),
#ifdef ASPECT_WITH_WORLD_BUILDER
    world_builder (parameters.world_builder_file != "" ?
                   std::make_shared<WorldBuilder::World>(parameters.world_builder_file) :
                   nullptr),
#endif
    boundary_heat_flux (BoundaryHeatFlux::create_boundary_heat_flux<dim>(prm)),
    time (numbers::signaling_nan<double>()),
    time_step (numbers::signaling_nan<double>()),
    old_time_step (numbers::signaling_nan<double>()),
    timestep_number (numbers::invalid_unsigned_int),
    nonlinear_iteration (numbers::invalid_unsigned_int),
    nonlinear_solver_failures (0),

    // We need to disable eliminate_refined_boundary_islands as this leads to
    // a deadlock for deal.II <= 9.2.0 as described in
    // https://github.com/geodynamics/aspect/issues/3604 when an
    // refined_island is at a periodic boundary. This flag is not too
    // important as it does not improve accuracy. Otherwise, these flags
    // correspond to smoothing_on_refinement|smoothing_on_coarsening.
    triangulation (mpi_communicator,
                   static_cast<typename Triangulation<dim>::MeshSmoothing>
                   (
                     Triangulation<dim>::limit_level_difference_at_vertices |
                     (Triangulation<dim>::eliminate_unrefined_islands |
                      Triangulation<dim>::eliminate_refined_inner_islands |
                      // Triangulation<dim>::eliminate_refined_boundary_islands |
                      Triangulation<dim>::do_not_produce_unrefined_islands)
                   )
                   ,
                   (parameters.stokes_solver_type == Parameters<dim>::StokesSolverType::block_gmg ||
                    parameters.stokes_solver_type == Parameters<dim>::StokesSolverType::default_solver
                    ?
                    static_cast<typename parallel::distributed::Triangulation<dim>::Settings>
                    (parallel::distributed::Triangulation<dim>::mesh_reconstruction_after_repartitioning |
                     parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy)
                    :
                    parallel::distributed::Triangulation<dim>::mesh_reconstruction_after_repartitioning)),

    mapping(construct_mapping<dim>(*geometry_model,*initial_topography_model)),

    // define the finite element
    finite_element(introspection.get_fes(), introspection.get_multiplicities()),

    dof_handler (triangulation),

    last_pressure_normalization_adjustment (numbers::signaling_nan<double>()),
    pressure_scaling (numbers::signaling_nan<double>()),

    rebuild_stokes_matrix (true),
    assemble_newton_stokes_matrix (true),
    assemble_newton_stokes_system (Parameters<dim>::is_defect_correction(parameters.nonlinear_solver)
                                   ?
                                   true
                                   :
                                   false),
    rebuild_stokes_preconditioner (true)
  {
    wall_timer.start();

    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
      {
        // only open the log file on processor 0, the other processors won't be
        // writing into the stream anyway
        log_file_stream.open(parameters.output_directory + "log.txt",
                             parameters.resume_computation ? std::ios_base::app : std::ios_base::out);

        // we already printed the header to the screen, so here we just dump it
        // into the log file.
        print_aspect_header(log_file_stream);
      }

    // now that we have output set up, we can start timer sections
    TimerOutput::Scope timer (computing_timer, "Initialization");


    // if any plugin wants access to the Simulator by deriving from SimulatorAccess, initialize it and
    // call the initialize() functions immediately after.
    //
    // up front, we can not know whether a plugin derives from
    // SimulatorAccess. all we have is a pointer to the base class of
    // each plugin type (the 'Interface' class in the namespace
    // corresponding to each plugin type), but this base class is not
    // derived from SimulatorAccess. in order to find out whether a
    // concrete plugin derives from this base (interface) class AND
    // the SimulatorAccess class via multiple inheritance, we need to
    // do a sideways dynamic_cast to this putative sibling of the
    // interface class, and investigate if the dynamic_cast
    // succeeds. if it succeeds, the dynamic_cast returns a non-nullptr
    // result, and we can test this in an if-statement. there is a nice
    // idiom whereby we can write
    //    if (SiblingClass *ptr = dynamic_cast<SiblingClass*>(ptr_to_base))
    //      ptr->do_something()
    // where we declare a variable *inside* the 'if' condition, and only
    // enter the code block guarded by the 'if' in case the so-declared
    // variable evaluates to something non-zero, which here means that
    // the dynamic_cast succeeded and returned the address of the sibling
    // object.
    //
    // we also need to let all models parse their parameters. this is done *after* setting
    // up their SimulatorAccess base class so that they can query, for example, the
    // geometry model's description of symbolic names for boundary parts. note that
    // the geometry model is the only model whose run time parameters are already read
    // at the time it is created
    if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(initial_topography_model.get()))
      sim->initialize_simulator (*this);
    initial_topography_model->initialize ();

    if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(geometry_model.get()))
      sim->initialize_simulator (*this);
    geometry_model->initialize ();

    if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(material_model.get()))
      sim->initialize_simulator (*this);
    material_model->parse_parameters (prm);
    material_model->initialize ();

    // Make sure that the material model supports elasticity when it is requested,
    // by checking that the material model creates the necessary material model
    // outputs.
    if (parameters.enable_elasticity)
      {
        // Set up outputs for one point
        MaterialModel::MaterialModelOutputs<dim> out(1,
                                                     introspection.n_compositional_fields);

        material_model->create_additional_named_outputs(out);

        MaterialModel::ElasticAdditionalOutputs<dim> *elastic_outputs = out.template get_additional_output<MaterialModel::ElasticAdditionalOutputs<dim>>();

        // Throw if the elastic_outputs do not exist
        AssertThrow(elastic_outputs != nullptr,
                    ExcMessage("Elasticity is enabled, but not supported by the material model."));
      }

    heating_model_manager.initialize_simulator(*this);
    heating_model_manager.parse_parameters (prm);

    if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(gravity_model.get()))
      sim->initialize_simulator (*this);
    gravity_model->parse_parameters (prm);
    gravity_model->initialize ();

    // Create the initial temperature and condition plugins, and then
    // initialize them. Some of these objects store std::shared_ptrs
    // to the initial temperature and composition objects, so it is
    // important that we have the pointers of the current class
    // already set by the time we call the initialize() functions.
    initial_temperature_manager
      = std::make_shared<InitialTemperature::Manager<dim>>();
    initial_composition_manager
      = std::make_shared<InitialComposition::Manager<dim>>();

    initial_temperature_manager->initialize_simulator(*this);
    initial_temperature_manager->parse_parameters (prm);

    initial_composition_manager->initialize_simulator(*this);
    initial_composition_manager->parse_parameters (prm);

    // Create a boundary temperature manager
    boundary_temperature_manager.initialize_simulator (*this);
    boundary_temperature_manager.parse_parameters (prm);

    if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(boundary_heat_flux.get()))
      sim->initialize_simulator (*this);
    boundary_heat_flux->parse_parameters (prm);
    boundary_heat_flux->initialize ();

    // Create a boundary composition manager
    boundary_composition_manager.initialize_simulator (*this);
    boundary_composition_manager.parse_parameters (prm);

    boundary_velocity_manager.initialize_simulator (*this);
    boundary_velocity_manager.parse_parameters (prm);

    // Make sure we only have a prescribed Stokes plugin if needed
    if (parameters.nonlinear_solver == NonlinearSolver::single_Advection_no_Stokes)
      {
        AssertThrow(prescribed_stokes_solution.get()!=nullptr,
                    ExcMessage("For the 'single Advection, no Stokes' solver scheme you need to provide a Stokes plugin!")
                   );
      }
    else
      {
        AssertThrow(prescribed_stokes_solution.get()==nullptr,
                    ExcMessage("The prescribed stokes plugin you selected only works with the solver "
                               "scheme 'single Advection, no Stokes'.")
                   );
      }

    if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(prescribed_stokes_solution.get()))
      sim->initialize_simulator (*this);
    if (prescribed_stokes_solution.get())
      {
        prescribed_stokes_solution->parse_parameters (prm);
        prescribed_stokes_solution->initialize ();
      }


    if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(adiabatic_conditions.get()))
      sim->initialize_simulator (*this);
    adiabatic_conditions->parse_parameters (prm);
    adiabatic_conditions->initialize ();

    // Create a boundary traction manager
    boundary_traction_manager.initialize_simulator (*this);
    boundary_traction_manager.parse_parameters (prm);

    // Initialize the mesh deformation handler
    if (parameters.mesh_deformation_enabled)
      {
        // Models with deformed boundaries require
        // the full A block to converge consistently
        parameters.use_full_A_block_preconditioner = true;

        // Allocate the MeshDeformationHandler object
        mesh_deformation = std::make_unique<MeshDeformation::MeshDeformationHandler<dim>>(*this);
        mesh_deformation->initialize_simulator(*this);
        mesh_deformation->parse_parameters(prm);
        mesh_deformation->initialize();
      }

    // Initialize the melt handler
    if (parameters.include_melt_transport)
      {
        // Models with melt transport require
        // the full A block to converge consistently
        parameters.use_full_A_block_preconditioner = true;

        melt_handler->initialize_simulator (*this);
        melt_handler->initialize();
      }

    // If the solver type is a Newton or defect correction type of solver, we need to set make sure
    // assemble_newton_stokes_system set to true.
    if (Parameters<dim>::is_defect_correction(parameters.nonlinear_solver))
      {
        assemble_newton_stokes_system = true;
        newton_handler->initialize_simulator(*this);
        newton_handler->parameters.parse_parameters(prm);
      }

    // choose the default solver and averaging scheme
    select_default_solver_and_averaging();

    if (parameters.stokes_solver_type == Parameters<dim>::StokesSolverType::block_gmg)
      {
        switch (parameters.stokes_velocity_degree)
          {
            case 2:
              stokes_matrix_free = std::make_unique<StokesMatrixFreeHandlerImplementation<dim,2>>(*this, parameters);
              break;
            case 3:
              stokes_matrix_free = std::make_unique<StokesMatrixFreeHandlerImplementation<dim,3>>(*this, parameters);
              break;
            default:
              AssertThrow(false, ExcMessage("The finite element degree for the Stokes system you selected is not supported yet."));
          }

        stokes_matrix_free->initialize_simulator(*this);
        stokes_matrix_free->parse_parameters(prm);
        stokes_matrix_free->initialize();
      }

    postprocess_manager.initialize_simulator (*this);
    postprocess_manager.parse_parameters (prm);

    if (postprocess_manager.template has_matching_active_plugin<Postprocess::Particles<dim>>())
      {
        particle_managers.resize(parameters.n_particle_managers);

        AssertThrow(particle_managers.size() <= ASPECT_MAX_NUM_PARTICLE_SYSTEMS,
                    ExcMessage("You have selected " + std::to_string(particle_managers.size()) + " particle managers, but ASPECT "
                               "has been compiled with a maximum of " + std::to_string(ASPECT_MAX_NUM_PARTICLE_SYSTEMS) + ". "
                               "Please recompile ASPECT with a higher value for ASPECT_MAX_NUM_PARTICLE_SYSTEMS. You can set a higher number "
                               "specifying the CMake variable -DASPECT_MAX_NUM_PARTICLE_SYSTEMS=<number>"));

        for (unsigned int particle_manager_index = 0 ; particle_manager_index < particle_managers.size(); ++particle_manager_index)
          {
            if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(&particle_managers[particle_manager_index]))
              sim->initialize_simulator (*this);

            particle_managers[particle_manager_index].parse_parameters(prm,particle_manager_index);
            particle_managers[particle_manager_index].initialize();
          }
      }

    mesh_refinement_manager.initialize_simulator (*this);
    mesh_refinement_manager.parse_parameters (prm);

    // VoF Must be initialized after mesh_refinement_manager, due to needing to check
    // for a mesh refinement strategy
    if (parameters.volume_of_fluid_tracking_enabled)
      {
        volume_of_fluid_handler->initialize (prm);
      }

    time_stepping_manager.initialize_simulator (*this);
    time_stepping_manager.parse_parameters (prm);

    lateral_averaging.initialize_simulator (*this);

    geometry_model->create_coarse_mesh (triangulation);
    Assert (triangulation.all_reference_cells_are_hyper_cube(),
            ExcMessage ("ASPECT only supports meshes that are composed of quadrilateral "
                        "or hexahedral cells."));
    global_Omega_diameter = GridTools::diameter (triangulation);

    // After creating the coarse mesh, initialize mapping cache if one is used
    if (MappingQCache<dim> *map = dynamic_cast<MappingQCache<dim>*>(&(*mapping)))
      map->initialize(MappingQGeneric<dim>(4), triangulation);

    bool dg_limiter_enabled = parameters.use_limiter_for_discontinuous_temperature_solution;
    for (unsigned int c=0; c<introspection.n_compositional_fields; ++c)
      dg_limiter_enabled = dg_limiter_enabled || parameters.use_limiter_for_discontinuous_composition_solution[c];

    // Check that DG limiters are only used with cartesian mapping
    if (dg_limiter_enabled)
      AssertThrow(geometry_model->natural_coordinate_system() == Utilities::Coordinates::CoordinateSystem::cartesian,
                  ExcMessage("The limiter for the discontinuous temperature and composition solutions "
                             "has not been tested in non-Cartesian geometries and currently requires "
                             "the use of a Cartesian geometry model."));

    std::set<types::boundary_id> open_velocity_boundary_indicators
      = geometry_model->get_used_boundary_indicators();
    for (const auto p : boundary_velocity_manager.get_prescribed_boundary_velocity_indicators())
      open_velocity_boundary_indicators.erase (p);
    for (const auto p : boundary_velocity_manager.get_zero_boundary_velocity_indicators())
      open_velocity_boundary_indicators.erase (p);
    for (const auto p : boundary_velocity_manager.get_tangential_boundary_velocity_indicators())
      open_velocity_boundary_indicators.erase (p);

    // Make sure that we do the pressure right-hand side modification correctly for periodic boundaries
    using periodic_boundary_set
      = std::set<std::pair<std::pair<types::boundary_id, types::boundary_id>, unsigned int>>;
    periodic_boundary_set pbs = geometry_model->get_periodic_boundary_pairs();
    for (const auto &pb : pbs)
      {
        open_velocity_boundary_indicators.erase (pb.first.first);
        open_velocity_boundary_indicators.erase (pb.first.second);
      }

    // We need to do the RHS compatibility modification, if the model is
    // compressible or compatible (in the case of melt transport), and
    // there is no open boundary to balance the pressure.
    do_pressure_rhs_compatibility_modification = ((material_model->is_compressible() && !parameters.include_melt_transport)
                                                  ||
                                                  (parameters.include_melt_transport && !material_model->is_compressible())
                                                  || parameters.enable_prescribed_dilation)
                                                 &&
                                                 (open_velocity_boundary_indicators.size() == 0);

    // make sure that we don't have to fill every column of the statistics
    // object in each time step.
    statistics.set_auto_fill_mode(true);

    // finally produce a record of the run-time parameters by writing
    // the currently used values into a file
    // Only write the parameter files on the root node to avoid file system conflicts
    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
      {
        std::ofstream prm_out ((parameters.output_directory + "parameters.prm"));
        AssertThrow (prm_out,
                     ExcMessage (std::string("Could not open file <") +
                                 parameters.output_directory + "parameters.prm>."));
        prm.print_parameters(prm_out, ParameterHandler::PRM);

        std::ofstream json_out ((parameters.output_directory + "parameters.json"));
        AssertThrow (json_out,
                     ExcMessage (std::string("Could not open file <") +
                                 parameters.output_directory + "parameters.json>."));
        prm.print_parameters(json_out, ParameterHandler::JSON);
      }

    // check that the boundary condition selection is consistent
    check_consistency_of_boundary_conditions();

    // check that the setup of equations, material models, and heating terms is consistent
    check_consistency_of_formulation();

    if (parameters.use_discontinuous_temperature_discretization || parameters.have_discontinuous_composition_discretization)
      CitationInfo::add("dg");

    // now that all member variables have been set up, also
    // connect the functions that will actually do the assembly
    set_assemblers();
  }


  /**
   * Destructor.
   */
  template <int dim>
  Simulator<dim>::~Simulator ()
  {
    // The particle_manager object is declared before the triangulation, and so
    // is destroyed after the latter. But it stores a pointer to the
    // triangulation and uses it during destruction. This results in
    // trouble. So destroy it first.
    particle_managers.clear();

    // wait if there is a thread that's still writing the statistics
    // object (set from the output_statistics() function)
    if (output_statistics_thread.joinable())
      output_statistics_thread.join();

    // If an exception is being thrown (for example due to AssertThrow()), we
    // might end up here with currently active timing sections. The destructor
    // of TimerOutput does MPI communication, which can lead to deadlocks,
    // hangs, or confusing MPI error messages. To avoid this, we can call
    // reset() to remove all open sections. In a normal run, we won't have any
    // active sessions, so this won't hurt to do:
    computing_timer.reset();
  }



  template <int dim>
  void
  Simulator<dim>::
  start_timestep ()
  {
    // first produce some output for the screen to show where we are
    {
      const char *unit = (parameters.convert_to_years ? "years" : "seconds");
      const double multiplier = (parameters.convert_to_years ? 1./year_in_seconds : 1.0);

      pcout << "*** Timestep " << timestep_number
            << ":  t=" << (time * multiplier) << ' ' << unit
            << ", dt=" << (time_step * multiplier) << ' ' << unit
            << std::endl;
    }

    nonlinear_iteration = 0;

    signals.start_timestep(*this);

    // Copy particle handler to restore particle location and properties
    // before repeating a timestep
    for (auto &particle_manager : particle_managers)
      particle_manager.backup_particles();


    // then interpolate the current boundary velocities. copy constraints
    // into current_constraints and then add to current_constraints
    compute_current_constraints ();

    // If needed, construct sparsity patterns and matrices with the current
    // constraints. Of course we need to force assembly too.
    if (rebuild_sparsity_and_matrices)
      {
        TimerOutput::Scope timer (computing_timer, "Setup matrices");

        rebuild_sparsity_and_matrices = false;
        setup_system_matrix (introspection.index_sets.system_partitioning);
        setup_system_preconditioner (introspection.index_sets.system_partitioning);
        rebuild_stokes_matrix = rebuild_stokes_preconditioner = true;
      }

    // notify different system components that we started the next time step
    // TODO: implement this for all plugins that might need it at one place.
    // Temperature BC are currently updated in compute_current_constraints
    geometry_model->update();
    material_model->update();
    gravity_model->update();
    heating_model_manager.update();
    adiabatic_conditions->update();
    mesh_refinement_manager.update();
    if (parameters.mesh_deformation_enabled)
      mesh_deformation->update();

    if (prescribed_stokes_solution.get())
      prescribed_stokes_solution->update();

    for (auto &particle_manager : particle_managers)
      particle_manager.update();

    // do the same for the traction boundary conditions and other things
    // that end up in the bilinear form. we update those that end up in
    // the constraints object when calling compute_current_constraints()
    // above
    boundary_traction_manager.update();
  }



  template <int dim>
  void
  Simulator<dim>::
  compute_current_constraints ()
  {
    // We put the constraints we compute into a separate AffineConstraints<double> so we can check
    // if the set of constraints has changed. If it did, we need to update the sparsity patterns.
    AffineConstraints<double> new_current_constraints(
#if DEAL_II_VERSION_GTE(9,6,0)
      dof_handler.locally_owned_dofs(),
#endif
      introspection.index_sets.system_relevant_set
    );
    new_current_constraints.merge (constraints);
    compute_current_velocity_boundary_constraints(new_current_constraints);

    // If there is a fixed boundary temperature or heat flux,
    // update the temperature boundary condition.
    boundary_temperature_manager.update();
    boundary_heat_flux->update();

    // If we do not want to prescribe Dirichlet boundary conditions on outflow boundaries,
    // we update the boundary indicators of all faces that belong to outflow boundaries
    // so that they are not in the list of fixed temperature boundary indicators any more.
    // We will undo this change in a later step, after the constraints have been set.
    // As long as we allow deal.II 8.5, we can not have boundary ids of more than 256,
    // so we want to offset them by 128 and not allow more than 128 boundary ids.
    const unsigned int boundary_id_offset = 128;
    if (!boundary_temperature_manager.allows_fixed_temperature_on_outflow_boundaries())
      replace_outflow_boundary_ids(boundary_id_offset);

    // if using continuous temperature FE, do the same for the temperature variable:
    // evaluate the current boundary temperature and add these constraints as well
    if (!parameters.use_discontinuous_temperature_discretization)
      {
        // obtain the boundary indicators that belong to Dirichlet-type
        // temperature boundary conditions and interpolate the temperature
        // there
        for (const auto p : boundary_temperature_manager.get_fixed_temperature_boundary_indicators())
          {
            VectorFunctionFromScalarFunctionObject<dim> vector_function_object(
              [&] (const dealii::Point<dim> &x) -> double
            {
              return boundary_temperature_manager.boundary_temperature(p, x);
            },
            introspection.component_masks.temperature.first_selected_component(),
            introspection.n_components);

            VectorTools::interpolate_boundary_values (*mapping,
                                                      dof_handler,
                                                      p,
                                                      vector_function_object,
                                                      new_current_constraints,
                                                      introspection.component_masks.temperature);
          }
      }

    if (!boundary_temperature_manager.allows_fixed_temperature_on_outflow_boundaries())
      restore_outflow_boundary_ids(boundary_id_offset);

    // If there are fixed boundary compositions,
    // update the composition boundary condition.
    boundary_composition_manager.update();

    // If we do not want to prescribe Dirichlet boundary conditions on outflow boundaries,
    // use the same trick for marking up outflow boundary conditions for compositional fields
    // as we did above already for the temperature.
    if (!boundary_composition_manager.allows_fixed_composition_on_outflow_boundaries())
      replace_outflow_boundary_ids(boundary_id_offset);

    // now do the same for the composition variables:
    {
      // obtain the boundary indicators that belong to Dirichlet-type
      // composition boundary conditions and interpolate the composition
      // there
      for (unsigned int c=0; c<introspection.n_compositional_fields; ++c)
        if (parameters.use_discontinuous_composition_discretization[c] == false)
          for (const auto p : boundary_composition_manager.get_fixed_composition_boundary_indicators())
            {
              VectorFunctionFromScalarFunctionObject<dim> vector_function_object(
                [&] (const Point<dim> &x) -> double
              {
                return boundary_composition_manager.boundary_composition(p, x, c);
              },
              introspection.component_masks.compositional_fields[c].first_selected_component(),
              introspection.n_components);

              VectorTools::interpolate_boundary_values (*mapping,
                                                        dof_handler,
                                                        p,
                                                        vector_function_object,
                                                        new_current_constraints,
                                                        introspection.component_masks.compositional_fields[c]);
            }
    }

    if (!boundary_composition_manager.allows_fixed_composition_on_outflow_boundaries())
      restore_outflow_boundary_ids(boundary_id_offset);

    if (parameters.include_melt_transport)
      melt_handler->add_current_constraints (new_current_constraints);

    // let plugins add more constraints if they so choose, then close the
    // constraints object
    signals.post_constraints_creation(*this, new_current_constraints);

    new_current_constraints.close();

    // Now check if the current_constraints we just computed changed from before. Do this before the melt handler
    // adds constraints for the melt cells, because they are allowed to change (and often do) without us having
    // to reinit the sparsity pattern.
    // If the mesh got refined and the size of the linear system changed, the old and new constraint
    // matrices will have different entries, and we can not easily compare them.
    // Note that we can have a situation where mesh_has_changed==true on only some of the machines,
    // so we need to do the MPI::sum below in all cases to avoid a deadlock.
    const bool mesh_has_changed = (current_constraints.get_local_lines().size()
                                   != new_current_constraints.get_local_lines().size())
                                  ||
                                  (current_constraints.get_local_lines()
                                   != new_current_constraints.get_local_lines());
    // Figure out if any entry that was constrained before is now no longer constrained:
    bool constraints_changed = false;

    if (mesh_has_changed)
      {
        // The mesh changed, so we know we need to update the constraints:
        constraints_changed = true;
      }
    else
      {
        // the mesh has not changed on our machine, so compare the constraints:
        for (auto &row: current_constraints.get_lines())
          {
            if (!new_current_constraints.is_constrained(row.index))
              {
                constraints_changed = true;
                break;
              }
          }
      }

    // If at least one processor has different constraints, force rebuilding the matrices:
    const bool any_constrained_dofs_set_changed = (Utilities::MPI::sum(constraints_changed ? 1 : 0,
                                                                       mpi_communicator)
                                                   > 0);
    if (any_constrained_dofs_set_changed)
      rebuild_sparsity_and_matrices = true;

#if DEAL_II_VERSION_GTE(9,6,0)
    current_constraints = std::move(new_current_constraints);
#else
    current_constraints.reinit (introspection.index_sets.system_relevant_set);
    current_constraints.copy_from(new_current_constraints);
#endif
    current_constraints.close();

    // TODO: We should use current_constraints.is_consistent_in_parallel()
    // here to assert that our constraints are consistent between
    // processors. This got removed in
    // https://github.com/geodynamics/aspect/pull/3282 because it triggered in
    // normal computations due to small floating point differences. See
    // https://github.com/geodynamics/aspect/issues/3248 for the discussion.
  }



  namespace
  {
    template <int dim>
    bool solver_scheme_solves_advection_equations(const Parameters<dim> &parameters)
    {
      // Check if we use a solver scheme that solves the advection equations
      switch (parameters.nonlinear_solver)
        {
          case Parameters<dim>::NonlinearSolver::Kind::single_Advection_single_Stokes:
          case Parameters<dim>::NonlinearSolver::Kind::iterated_Advection_and_Stokes:
          case Parameters<dim>::NonlinearSolver::Kind::single_Advection_iterated_Stokes:
          case Parameters<dim>::NonlinearSolver::Kind::no_Advection_iterated_defect_correction_Stokes:
          case Parameters<dim>::NonlinearSolver::Kind::single_Advection_iterated_defect_correction_Stokes:
          case Parameters<dim>::NonlinearSolver::Kind::iterated_Advection_and_defect_correction_Stokes:
          case Parameters<dim>::NonlinearSolver::Kind::iterated_Advection_and_Newton_Stokes:
          case Parameters<dim>::NonlinearSolver::Kind::single_Advection_iterated_Newton_Stokes:
          case Parameters<dim>::NonlinearSolver::Kind::single_Advection_no_Stokes:
            return true;

          case Parameters<dim>::NonlinearSolver::Kind::no_Advection_iterated_Stokes:
          case Parameters<dim>::NonlinearSolver::Kind::no_Advection_single_Stokes:
          case Parameters<dim>::NonlinearSolver::Kind::first_timestep_only_single_Stokes:
          case Parameters<dim>::NonlinearSolver::Kind::no_Advection_no_Stokes:
            return false;
        }
      Assert(false, ExcNotImplemented());
      return false;
    }



    template <int dim>
    bool solver_scheme_solves_stokes_equations(const Parameters<dim> &parameters)
    {
      // Check if we use a solver scheme that solves the Stokes equations
      switch (parameters.nonlinear_solver)
        {
          case Parameters<dim>::NonlinearSolver::Kind::single_Advection_single_Stokes:
          case Parameters<dim>::NonlinearSolver::Kind::iterated_Advection_and_Stokes:
          case Parameters<dim>::NonlinearSolver::Kind::single_Advection_iterated_Stokes:
          case Parameters<dim>::NonlinearSolver::Kind::no_Advection_iterated_Stokes:
          case Parameters<dim>::NonlinearSolver::Kind::no_Advection_single_Stokes:
          case Parameters<dim>::NonlinearSolver::Kind::no_Advection_iterated_defect_correction_Stokes:
          case Parameters<dim>::NonlinearSolver::Kind::single_Advection_iterated_defect_correction_Stokes:
          case Parameters<dim>::NonlinearSolver::Kind::iterated_Advection_and_defect_correction_Stokes:
          case Parameters<dim>::NonlinearSolver::Kind::iterated_Advection_and_Newton_Stokes:
          case Parameters<dim>::NonlinearSolver::Kind::single_Advection_iterated_Newton_Stokes:
          case Parameters<dim>::NonlinearSolver::Kind::first_timestep_only_single_Stokes:
            return true;

          case Parameters<dim>::NonlinearSolver::Kind::single_Advection_no_Stokes:
          case Parameters<dim>::NonlinearSolver::Kind::no_Advection_no_Stokes:
            return false;
        }
      Assert(false, ExcNotImplemented());
      return false;
    }



    template <int dim>
    bool compositional_field_needs_matrix_block(const Introspection<dim> &introspection, const unsigned int composition_index)
    {
      const typename Simulator<dim>::AdvectionField adv_field (Simulator<dim>::AdvectionField::composition(composition_index));
      switch (adv_field.advection_method(introspection))
        {
          case Parameters<dim>::AdvectionFieldMethod::fem_field:
          case Parameters<dim>::AdvectionFieldMethod::fem_melt_field:
          case Parameters<dim>::AdvectionFieldMethod::fem_darcy_field:
          case Parameters<dim>::AdvectionFieldMethod::prescribed_field_with_diffusion:
            return true;
          case Parameters<dim>::AdvectionFieldMethod::particles:
          case Parameters<dim>::AdvectionFieldMethod::volume_of_fluid:
          case Parameters<dim>::AdvectionFieldMethod::static_field:
          case Parameters<dim>::AdvectionFieldMethod::prescribed_field:
            break;
          default:
            Assert (false, ExcNotImplemented());
        }
      return false;
    }



    template <int dim>
    bool compositional_fields_need_matrix_block(const Introspection<dim> &introspection)
    {
      // Check if any compositional field method actually requires a matrix block
      // (as opposed to all are advected by other means or prescribed fields)
      for (unsigned int c=0; c<introspection.n_compositional_fields; ++c)
        {
          if (compositional_field_needs_matrix_block(introspection, c))
            return true;
        }
      return false;
    }
  }



  template <int dim>
  Table<2,DoFTools::Coupling>
  Simulator<dim>::
  setup_system_matrix_coupling () const
  {
    Table<2,DoFTools::Coupling> coupling(introspection.n_components,
                                         introspection.n_components);

    // Start by assuming nothing couples
    coupling.fill (DoFTools::none);

    const typename Introspection<dim>::ComponentIndices &x
      = introspection.component_indices;

    // Determine which blocks in the Stokes blocks
    // of the matrix are in use. At the moment this distinguishes
    // 3 solvers: melt transport, matrix-free multigrid, and the default
    // matrix based algebraic multigrid.
    if (solver_scheme_solves_stokes_equations(parameters))
      {
        // The matrix-free solver does not work with melt transport
        Assert(!(parameters.include_melt_transport && stokes_matrix_free),
               ExcNotImplemented());

        if (stokes_matrix_free)
          {
            // nothing couples in the matrix free solver
          }
        else if (parameters.include_melt_transport)
          {
            // For the melt transport solver all velocities and pressures couple with themselves.
            // Additionally solid velocities couple with all pressures, and all pressures
            // couple with solid velocities.

            const unsigned int first_fluid_c_i = introspection.variable("fluid velocity").first_component_index;

            for (unsigned int d=0; d<dim; ++d)
              {
                for (unsigned int c=0; c<dim; ++c)
                  {
                    coupling[x.velocities[c]][x.velocities[d]] = DoFTools::always;
                    coupling[first_fluid_c_i+c][first_fluid_c_i+d] = DoFTools::always;
                  }

                coupling[x.velocities[d]][
                  introspection.variable("compaction pressure").first_component_index] = DoFTools::always;
                coupling[introspection.variable("compaction pressure").first_component_index]
                [x.velocities[d]]
                  = DoFTools::always;
                coupling[x.velocities[d]]
                [introspection.variable("fluid pressure").first_component_index]
                  = DoFTools::always;
                coupling[introspection.variable("fluid pressure").first_component_index]
                [x.velocities[d]]
                  = DoFTools::always;
              }

            coupling[introspection.variable("fluid pressure").first_component_index]
            [introspection.variable("fluid pressure").first_component_index]
              = DoFTools::always;
            coupling[introspection.variable("compaction pressure").first_component_index]
            [introspection.variable("compaction pressure").first_component_index]
              = DoFTools::always;
          }
        else
          {
            // The AMG matrix based solver: all velocities couple with all velocities,
            // pressure couples with all velocities and the other way around,
            // and pressures only couple with themselves for equal order elements
            for (unsigned int c=0; c<dim; ++c)
              for (unsigned int d=0; d<dim; ++d)
                coupling[x.velocities[c]][x.velocities[d]] = DoFTools::always;

            for (unsigned int d=0; d<dim; ++d)
              {
                coupling[x.velocities[d]][x.pressure] = DoFTools::always;
                coupling[x.pressure][x.velocities[d]] = DoFTools::always;
              }

            // For equal-order interpolation, we need a stabilization term
            // in the bottom right of Stokes matrix. Make sure we have the
            // necessary entries.
            if (parameters.use_equal_order_interpolation_for_stokes == true)
              coupling[x.pressure][x.pressure] = DoFTools::always;
          }
      }

    // Only enable temperature coupling if temperature block is needed
    if (solver_scheme_solves_advection_equations(parameters)
        &&
        parameters.temperature_method != Parameters<dim>::AdvectionFieldMethod::prescribed_field
        &&
        parameters.temperature_method != Parameters<dim>::AdvectionFieldMethod::static_field)
      coupling[x.temperature][x.temperature] = DoFTools::always;

    // Only enable composition coupling if a composition block is needed
    if (solver_scheme_solves_advection_equations(parameters)
        &&
        compositional_fields_need_matrix_block(introspection))
      {
        // We reuse matrix blocks for compositional fields, but we need different blocks for each base_element.
        // All other matrix blocks are left empty to save memory.
        for (const unsigned int base_element_index : introspection.get_composition_base_element_indices())
          {
            bool block_needed = false;
            for (const unsigned int c : introspection.get_compositional_field_indices_with_base_element(base_element_index))
              if (compositional_field_needs_matrix_block(introspection, c))
                block_needed = true;

            if (block_needed)
              {
                const unsigned int first_c = introspection.get_compositional_field_indices_with_base_element(base_element_index).front();
                coupling[x.compositional_fields[first_c]][x.compositional_fields[first_c]] = DoFTools::always;
              }
          }
      }

    // If we are using volume of fluid interface tracking, create a matrix block in the
    // field corresponding to the volume fraction.
    if (parameters.volume_of_fluid_tracking_enabled)
      {
        const unsigned int volume_of_fluid_block = volume_of_fluid_handler->field_struct_for_field_index(0)
                                                   .volume_fraction.first_component_index;
        coupling[volume_of_fluid_block][volume_of_fluid_block] = DoFTools::always;
      }

    return coupling;
  }


  template <int dim>
  void
  Simulator<dim>::
  setup_system_matrix (const std::vector<IndexSet> &system_partitioning)
  {
    system_matrix.clear ();

    const Table<2,DoFTools::Coupling> coupling = setup_system_matrix_coupling();
    LinearAlgebra::BlockDynamicSparsityPattern sp;

    sp.reinit (system_partitioning,
               system_partitioning,
               introspection.index_sets.system_relevant_partitioning,
               mpi_communicator);


    if ((parameters.use_discontinuous_temperature_discretization) ||
        (parameters.have_discontinuous_composition_discretization) ||
        (parameters.volume_of_fluid_tracking_enabled))
      {
        Table<2,DoFTools::Coupling> face_coupling (introspection.n_components,
                                                   introspection.n_components);
        face_coupling.fill (DoFTools::none);

        const typename Introspection<dim>::ComponentIndices &x
          = introspection.component_indices;
        if (parameters.use_discontinuous_temperature_discretization &&
            solver_scheme_solves_advection_equations(parameters) &&
            parameters.temperature_method != Parameters<dim>::AdvectionFieldMethod::prescribed_field &&
            parameters.temperature_method != Parameters<dim>::AdvectionFieldMethod::static_field)
          face_coupling[x.temperature][x.temperature] = DoFTools::always;

        for (const unsigned int base_element_index : introspection.get_composition_base_element_indices())
          {
            bool block_needed = false;
            for (const unsigned int c : introspection.get_compositional_field_indices_with_base_element(base_element_index))
              if (compositional_field_needs_matrix_block(introspection, c))
                block_needed = true;

            const unsigned int first_c = introspection.get_compositional_field_indices_with_base_element(base_element_index).front();

            if (parameters.use_discontinuous_composition_discretization[first_c]
                && solver_scheme_solves_advection_equations(parameters)
                && block_needed)
              face_coupling[x.compositional_fields[first_c]][x.compositional_fields[first_c]] = DoFTools::always;
          }

        if (parameters.volume_of_fluid_tracking_enabled)
          {
            const unsigned int volume_of_fluid_block = volume_of_fluid_handler->field_struct_for_field_index(0)
                                                       .volume_fraction.first_component_index;
            face_coupling[volume_of_fluid_block][volume_of_fluid_block] = DoFTools::always;
          }

        DoFTools::make_flux_sparsity_pattern (dof_handler,
                                              sp,
                                              current_constraints, false,
                                              coupling,
                                              face_coupling,
                                              Utilities::MPI::
                                              this_mpi_process(mpi_communicator));

        if (solver_scheme_solves_advection_equations(parameters))
          {
            // If we solve for more than one compositional field make sure we keep constrained entries
            // to allow different boundary conditions for different fields. In order to keep constrained
            // entries we need to repeat the call to DoFTools::make_flux_sparsity_pattern above,
            // but with 'true' as 4th argument, and only for the composition block.
            Table<2,DoFTools::Coupling> composition_coupling(introspection.n_components,
                                                             introspection.n_components);
            composition_coupling.fill (DoFTools::none);

            for (const unsigned int base_element_index : introspection.get_composition_base_element_indices())
              {
                const unsigned int first_c = introspection.get_compositional_field_indices_with_base_element(base_element_index).front();
                const unsigned int component = x.compositional_fields[first_c];
                composition_coupling[component][component] = coupling[component][component];

                const unsigned int block = introspection.block_indices.compositional_field_sparsity_pattern[first_c];
                sp.block(block,block).reinit(sp.block(block,block).locally_owned_range_indices(),
                                             sp.block(block,block).locally_owned_domain_indices());

                DoFTools::make_flux_sparsity_pattern (dof_handler,
                                                      sp,
                                                      current_constraints, true,
                                                      composition_coupling,
                                                      face_coupling,
                                                      Utilities::MPI::
                                                      this_mpi_process(mpi_communicator));
              }
          }
      }
    else
      {
        DoFTools::make_sparsity_pattern (dof_handler,
                                         coupling, sp,
                                         current_constraints, false,
                                         Utilities::MPI::
                                         this_mpi_process(mpi_communicator));

        // If we solve for more than one compositional field make sure we keep constrained entries
        // to allow different boundary conditions for different fields. In order to keep constrained
        // entries we need to repeat the call to DoFTools::make_sparsity_pattern above,
        // but with 'true' as 4th argument, and only for the composition block.
        if (solver_scheme_solves_advection_equations(parameters)
            &&
            compositional_fields_need_matrix_block(introspection))
          {
            Table<2,DoFTools::Coupling> composition_coupling(introspection.n_components,
                                                             introspection.n_components);
            composition_coupling.fill (DoFTools::none);

            const unsigned int component = introspection.component_indices.compositional_fields[0];
            composition_coupling[component][component] = coupling[component][component];

            const unsigned int block = introspection.get_components_to_blocks()[component];
            sp.block(block,block).reinit(sp.block(block,block).locally_owned_range_indices(),
                                         sp.block(block,block).locally_owned_domain_indices());

            DoFTools::make_sparsity_pattern (dof_handler,
                                             composition_coupling, sp,
                                             current_constraints, true,
                                             Utilities::MPI::
                                             this_mpi_process(mpi_communicator));
          }
      }

    sp.compress();

    // We may only allocate some of the matrix blocks, but the sparsity pattern
    // will still create entries for hanging nodes and boundary conditions.
    // These are unnecessary and are removed here.
    for (unsigned int i=0; i<introspection.n_components; ++i)
      {
        if (coupling[i][i] == DoFTools::none)
          {
            // The pressure is special, because it does not couple with itself,
            // but if the velocity block is assembled we also need to keep pressure
            // hanging node constraints. Thus skip clearing the
            // sparsity pattern in that case.
            if ((i == introspection.component_indices.pressure) &&
                coupling[introspection.component_indices.velocities[0]][introspection.component_indices.velocities[0]] != DoFTools::none)
              continue;

            const unsigned int block = introspection.get_components_to_blocks()[i];

            // TODO: using clear() would be nice here but clear() also resets the
            // size, so just reinit():
            sp.block(block,block).reinit(sp.block(block,block).locally_owned_range_indices(),
                                         sp.block(block,block).locally_owned_domain_indices());
            sp.block(block,block).compress();
          }
      }

    system_matrix.reinit (sp);
  }



  template <int dim>
  void Simulator<dim>::
  setup_system_preconditioner (const std::vector<IndexSet> &system_partitioning)
  {
    Amg_preconditioner.reset ();
    Mp_preconditioner.reset ();
    system_preconditioner_matrix.clear ();

    // The preconditioner matrix is only used for the Stokes block (velocity and Schur complement)
    // and only needed if we actually solve iteratively and matrix-based
    if (solver_scheme_solves_stokes_equations(parameters) == false)
      return;
    else if (parameters.stokes_solver_type == Parameters<dim>::StokesSolverType::block_gmg)
      return;
    else if (parameters.stokes_solver_type == Parameters<dim>::StokesSolverType::block_amg)
      {
        // continue below
      }
    else if (parameters.stokes_solver_type == Parameters<dim>::StokesSolverType::direct_solver)
      return;
    else
      AssertThrow(false, ExcNotImplemented());

    Table<2,DoFTools::Coupling> coupling (introspection.n_components,
                                          introspection.n_components);
    coupling.fill (DoFTools::none);

    const typename Introspection<dim>::ComponentIndices &x
      = introspection.component_indices;

    // velocity-velocity block (only block diagonal) is only
    // needed if we use the simplified A block preconditioner
    if (parameters.use_full_A_block_preconditioner == false)
      for (unsigned int d=0; d<dim; ++d)
        coupling[x.velocities[d]][x.velocities[d]] = DoFTools::always;

    // TODO: The newton handler always assembles the preconditioner matrix
    // even if it is not used. Fix that by modifying the newton assemblers.
    if (newton_handler.get() != nullptr)
      for (unsigned int d=0; d<dim; ++d)
        coupling[x.velocities[d]][x.velocities[d]] = DoFTools::always;

    // Schur complement block (pressure - pressure):
    if (parameters.include_melt_transport)
      {
        coupling[introspection.variable("fluid pressure").first_component_index]
        [introspection.variable("fluid pressure").first_component_index] = DoFTools::always;
        coupling[introspection.variable("compaction pressure").first_component_index]
        [introspection.variable("compaction pressure").first_component_index] = DoFTools::always;
        coupling[introspection.variable("compaction pressure").first_component_index]
        [introspection.variable("fluid pressure").first_component_index] = DoFTools::always;
        coupling[introspection.variable("fluid pressure").first_component_index]
        [introspection.variable("compaction pressure").first_component_index] = DoFTools::always;
      }
    else
      coupling[x.pressure][x.pressure] = DoFTools::always;
    // the system_preconditioner matrix is only used for the Stokes
    // system as there we use a separate matrix for preconditioning
    // than the system_matrix object). for temperature and
    // compositional fields, we precondition with the respective block
    // of system_matrix. consequently, there is no need to allocate
    // memory for these blocks in the system_preconditioner matrix and
    // its sparsity pattern here -- the corresponding entries of
    // 'coupling' simply remain at DoFTools::none

    LinearAlgebra::BlockDynamicSparsityPattern sp;

    sp.reinit (system_partitioning,
               system_partitioning,
               introspection.index_sets.system_relevant_partitioning,
               mpi_communicator);

    DoFTools::make_sparsity_pattern (dof_handler,
                                     coupling, sp,
                                     current_constraints, false,
                                     Utilities::MPI::
                                     this_mpi_process(mpi_communicator));


    sp.compress();

    // We are not interested in temperature and composition matrices for the
    // preconditioner matrix. But even though we specify a coupling of
    // DoFTools::none, entries for constrained entries for boundary conditions
    // and hanging nodes are being created by make_sparsity_pattern. These are
    // unnecessary, so we remove those entries here.
    {
      // temperature:
      const unsigned int block_idx = introspection.block_indices.temperature;
      // TODO: using clear() would be nice here but clear() also resets the
      // size, so just reinit():
      sp.block(block_idx, block_idx).reinit(sp.block(block_idx, block_idx).locally_owned_range_indices(),
                                            sp.block(block_idx, block_idx).locally_owned_domain_indices());
      sp.block(block_idx, block_idx).compress();
    }
    // compositions:
    for (unsigned int c=0; c<introspection.n_compositional_fields; ++c)
      {
        const unsigned int block_idx = introspection.block_indices.compositional_fields[c];
        // TODO: using clear() would be nice here but clear() also resets the
        // size, so just reinit():
        sp.block(block_idx, block_idx).reinit(sp.block(block_idx, block_idx).locally_owned_range_indices(),
                                              sp.block(block_idx, block_idx).locally_owned_domain_indices());
        sp.block(block_idx, block_idx).compress();
      }

    system_preconditioner_matrix.reinit (sp);
    if (parameters.use_bfbt)
      inverse_lumped_mass_matrix.reinit(introspection.index_sets.stokes_partitioning);
  }


  template <int dim>
  void Simulator<dim>::compute_initial_velocity_boundary_constraints (AffineConstraints<double> &constraints)
  {

    // This needs to happen after the periodic constraints are added:
    setup_nullspace_constraints(constraints);

    // then compute constraints for the velocity. the constraints we compute
    // here are the ones that are the same for all following time steps. in
    // addition, we may be computing constraints from boundary values for the
    // velocity that are different between time steps. these are then put
    // into current_constraints in start_timestep().
    signals.pre_compute_no_normal_flux_constraints(triangulation);
    {
      // do the interpolation for zero velocity
      for (const auto p : boundary_velocity_manager.get_zero_boundary_velocity_indicators())
        VectorTools::interpolate_boundary_values (*mapping,
                                                  dof_handler,
                                                  p,
                                                  Functions::ZeroFunction<dim>(introspection.n_components),
                                                  constraints,
                                                  introspection.component_masks.velocities);


      // do the same for no-normal-flux boundaries
      VectorTools::compute_no_normal_flux_constraints (dof_handler,
                                                       /* first_vector_component= */
                                                       introspection.component_indices.velocities[0],
                                                       boundary_velocity_manager.get_tangential_boundary_velocity_indicators(),
                                                       constraints,
                                                       *mapping);
    }


  }

  template <int dim>
  void Simulator<dim>::compute_current_velocity_boundary_constraints (AffineConstraints<double> &constraints)
  {
    // set the current time and do the interpolation
    // for the prescribed velocity fields
    boundary_velocity_manager.update();

    // put boundary conditions into constraints object for each boundary
    for (const auto boundary_id: boundary_velocity_manager.get_prescribed_boundary_velocity_indicators())
      {
        Utilities::VectorFunctionFromVelocityFunctionObject<dim> vel
        (introspection.n_components,
         [&] (const dealii::Point<dim> &x) -> Tensor<1,dim>
        {
          if (!assemble_newton_stokes_system || (assemble_newton_stokes_system && nonlinear_iteration == 0))
            return boundary_velocity_manager.boundary_velocity(boundary_id, x);
          else
            return Tensor<1,dim>();

          return Tensor<1,dim>();
        });

        VectorTools::interpolate_boundary_values (*mapping,
                                                  dof_handler,
                                                  boundary_id,
                                                  vel,
                                                  constraints,
                                                  boundary_velocity_manager.get_component_mask(boundary_id));
      }
  }


  template <int dim>
  void Simulator<dim>::setup_dofs ()
  {
    signals.edit_parameters_pre_setup_dofs(*this, parameters);

    TimerOutput::Scope timer (computing_timer, "Setup dof systems");

    dof_handler.distribute_dofs(finite_element);

    // Renumber the DoFs hierarchical so that we get the
    // same numbering if we resume the computation. This
    // is because the numbering depends on the order the
    // cells are created.
    DoFRenumbering::hierarchical (dof_handler);
    DoFRenumbering::component_wise (dof_handler,
                                    introspection.get_components_to_blocks());

    // set up the introspection object that stores all sorts of
    // information about components of the finite element, component
    // masks, etc
    setup_introspection();

    // print dof numbers. Do so with 1000s separator since they are frequently
    // large
    {
      std::locale s = pcout.get_stream().getloc();
      // Creating std::locale with an empty string previously caused problems
      // on some platforms, so the functionality to catch the exception and ignore
      // is kept here, even though explicitly setting a facet should always work.
      try
        {
          // Imbue the stream with a locale that does the right thing. The
          // locale is responsible for later deleting the object pointed
          // to by the last argument (the "facet"), see
          // https://en.cppreference.com/w/cpp/locale/locale/locale
          pcout.get_stream().imbue(std::locale(std::locale(),
                                               new aspect::Utilities::ThousandSep));
        }
      catch (const std::runtime_error &e)
        {
          // If the locale doesn't work, just give up
        }

      pcout << "Number of active cells: "
            << triangulation.n_global_active_cells()
            << " (on "
            << triangulation.n_global_levels()
            << " levels)"
            << std::endl
            << "Number of degrees of freedom: "
            << dof_handler.n_dofs()
            << " ("
            << introspection.system_dofs_per_block[0];
      for (unsigned int b=1; b<introspection.system_dofs_per_block.size(); ++b)
        pcout << '+' << introspection.system_dofs_per_block[b];
      pcout <<')'
            << std::endl
            << std::endl;

      pcout.get_stream().imbue(s);
    }

    // We need to set up the mesh deformation degrees of freedom first if mesh deformation
    // is active, since the mapping must be in place before applying boundary
    // conditions that rely on it (such as no flux BCs).
    if (parameters.mesh_deformation_enabled)
      mesh_deformation->setup_dofs();


    // Reconstruct the constraint-matrix:
#if DEAL_II_VERSION_GTE(9,6,0)
    constraints.reinit (dof_handler.locally_owned_dofs(), introspection.index_sets.system_relevant_set);
#else
    constraints.reinit(introspection.index_sets.system_relevant_set);
#endif

    // Set up the constraints for periodic boundary conditions:

    // Note: this has to happen _before_ we do hanging node constraints,
    // because inconsistent constraints could be generated in parallel otherwise.
    geometry_model->make_periodicity_constraints(dof_handler,
                                                 constraints);

    //  Make hanging node constraints:
    DoFTools::make_hanging_node_constraints (dof_handler,
                                             constraints);


    compute_initial_velocity_boundary_constraints(constraints);
    constraints.close();
    signals.post_compute_no_normal_flux_constraints(triangulation);

    // Finally initialize vectors. We delay construction of the sparsity
    // patterns and matrices until we have current_constraints.
    rebuild_sparsity_and_matrices = true;

    system_rhs.reinit(introspection.index_sets.system_partitioning, mpi_communicator);
    solution.reinit(introspection.index_sets.system_partitioning, introspection.index_sets.system_relevant_partitioning, mpi_communicator);
    old_solution.reinit(introspection.index_sets.system_partitioning, introspection.index_sets.system_relevant_partitioning, mpi_communicator);
    old_old_solution.reinit(introspection.index_sets.system_partitioning, introspection.index_sets.system_relevant_partitioning, mpi_communicator);
    current_linearization_point.reinit (introspection.index_sets.system_partitioning, introspection.index_sets.system_relevant_partitioning, mpi_communicator);

    if (parameters.use_operator_splitting)
      operator_split_reaction_vector.reinit (introspection.index_sets.system_partitioning, introspection.index_sets.system_relevant_partitioning, mpi_communicator);

    if (do_pressure_rhs_compatibility_modification)
      pressure_shape_function_integrals.reinit (introspection.index_sets.system_partitioning, mpi_communicator);

    rebuild_stokes_matrix         = true;
    rebuild_stokes_preconditioner = true;

    // Setup matrix-free dofs
    if (stokes_matrix_free)
      stokes_matrix_free->setup_dofs();
  }





  template <int dim>
  void Simulator<dim>::setup_introspection ()
  {
    // compute the various partitionings between processors and blocks
    // of vectors and matrices
    introspection.system_dofs_per_block = DoFTools::count_dofs_per_fe_block (dof_handler,
                                                                             introspection.get_components_to_blocks());

    {
      IndexSet system_index_set = dof_handler.locally_owned_dofs();
      introspection.index_sets.system_partitioning = system_index_set.split_by_block(introspection.system_dofs_per_block);

#if DEAL_II_VERSION_GTE(9,7,0)
      introspection.index_sets.system_relevant_set = DoFTools::extract_locally_relevant_dofs (dof_handler);
#else
      DoFTools::extract_locally_relevant_dofs (dof_handler,
                                               introspection.index_sets.system_relevant_set);
#endif
      introspection.index_sets.system_relevant_partitioning =
        introspection.index_sets.system_relevant_set.split_by_block(introspection.system_dofs_per_block);

      if (!parameters.include_melt_transport)
        {
          if (parameters.use_direct_stokes_solver)
            introspection.index_sets.locally_owned_pressure_dofs = system_index_set & Utilities::extract_locally_active_dofs_with_component(dof_handler, introspection.component_masks.pressure);
          else
            introspection.index_sets.locally_owned_pressure_dofs = introspection.index_sets.system_partitioning[introspection.block_indices.pressure];
        }

      if (parameters.include_melt_transport)
        {
          introspection.index_sets.locally_owned_melt_pressure_dofs = system_index_set & Utilities::extract_locally_active_dofs_with_component(dof_handler,
                                                                      introspection.variable("fluid pressure").component_mask|
                                                                      introspection.variable("compaction pressure").component_mask);
          introspection.index_sets.locally_owned_fluid_pressure_dofs = system_index_set & Utilities::extract_locally_active_dofs_with_component(dof_handler,
                                                                       introspection.variable("fluid pressure").component_mask);
        }



      introspection.index_sets.stokes_partitioning.clear ();
      introspection.index_sets.stokes_partitioning.push_back(introspection.index_sets.system_partitioning[introspection.block_indices.velocities]);
      if (!parameters.use_direct_stokes_solver)
        {
          if (parameters.include_melt_transport)
            // p_f and p_c are in the same block
            introspection.index_sets.stokes_partitioning.push_back(introspection.index_sets.system_partitioning[introspection.variable("fluid pressure").block_index]);
          else
            introspection.index_sets.stokes_partitioning.push_back(introspection.index_sets.system_partitioning[introspection.block_indices.pressure]);
        }

      Assert(!parameters.use_direct_stokes_solver ||
             parameters.include_melt_transport ||
             (introspection.block_indices.velocities == introspection.block_indices.pressure),
             ExcInternalError());
    }

  }



  template <int dim>
  void Simulator<dim>::postprocess ()
  {
    TimerOutput::Scope timer (computing_timer, "Postprocessing");
    pcout << "   Postprocessing:" << std::endl;

    // run all the postprocessing routines and then write
    // the current state of the statistics table to a file
    std::list<std::pair<std::string,std::string>>
    output_list = postprocess_manager.execute (statistics);

    // if we are on processor zero, print to screen
    // whatever the postprocessors have generated
    if (Utilities::MPI::this_mpi_process(mpi_communicator)==0)
      {
        // determine the width of the first column of text so that
        // everything gets nicely aligned; then output everything
        {
          unsigned int width = 0;
          for (const auto &p : output_list)
            width = std::max<unsigned int> (width, p.first.size());

          for (const auto &p : output_list)
            pcout << "     "
                  << std::left
                  << std::setw(width)
                  << p.first
                  << ' '
                  << p.second
                  << std::endl;
        }

        pcout << std::endl;
      }

    // finally, write the entire set of current results to disk
    output_statistics();
  }


  template <int dim>
  void Simulator<dim>::refine_mesh (const unsigned int max_grid_level)
  {
    parallel::distributed::SolutionTransfer<dim,LinearAlgebra::BlockVector>
    system_trans(dof_handler);

    std::unique_ptr<parallel::distributed::SolutionTransfer<dim,LinearAlgebra::Vector>>
    mesh_deformation_trans;

    {
      TimerOutput::Scope timer (computing_timer, "Refine mesh structure, part 1");

      Vector<float> estimated_error_per_cell (triangulation.n_active_cells());
      mesh_refinement_manager.execute (estimated_error_per_cell);

      if (parameters.adapt_by_fraction_of_cells)
        parallel::distributed::GridRefinement::
        refine_and_coarsen_fixed_number   (triangulation,
                                           estimated_error_per_cell,
                                           parameters.refinement_fraction,
                                           parameters.coarsening_fraction);
      else
        parallel::distributed::GridRefinement::
        refine_and_coarsen_fixed_fraction (triangulation,
                                           estimated_error_per_cell,
                                           parameters.refinement_fraction,
                                           parameters.coarsening_fraction);

      mesh_refinement_manager.tag_additional_cells ();


      // clear refinement flags if parameter.refinement_fraction=0.0
      if (parameters.refinement_fraction==0.0)
        {
          for (const auto &cell : triangulation.active_cell_iterators())
            cell->clear_refine_flag ();
        }
      else
        {
          // limit maximum refinement level
          if (triangulation.n_levels() > max_grid_level)
            for (typename Triangulation<dim>::active_cell_iterator
                 cell = triangulation.begin_active(max_grid_level);
                 cell != triangulation.end(); ++cell)
              cell->clear_refine_flag ();
        }

      // clear coarsening flags if parameter.coarsening_fraction=0.0
      if (parameters.coarsening_fraction==0.0)
        {
          for (const auto &cell : triangulation.active_cell_iterators())
            cell->clear_coarsen_flag ();
        }
      else
        {
          // limit minimum refinement level
          for (typename Triangulation<dim>::active_cell_iterator
               cell = triangulation.begin_active(0);
               cell != triangulation.end_active(parameters.min_grid_level); ++cell)
            cell->clear_coarsen_flag ();
        }


      // Next set up whatever is necessary to transfer the solution from old
      // to new mesh:
      std::vector<const LinearAlgebra::BlockVector *> x_system
        = { &solution, &old_solution };

      if (parameters.mesh_deformation_enabled)
        x_system.push_back(&mesh_deformation->mesh_velocity);

      std::vector<const LinearAlgebra::Vector *> x_fs_system;
      if (parameters.mesh_deformation_enabled)
        {
          x_fs_system.push_back (&mesh_deformation->mesh_displacements);
          x_fs_system.push_back (&mesh_deformation->old_mesh_displacements);
          x_fs_system.push_back (&mesh_deformation->initial_topography);
          mesh_deformation_trans
            = std::make_unique<parallel::distributed::SolutionTransfer<dim,LinearAlgebra::Vector>>
              (mesh_deformation->mesh_deformation_dof_handler);
        }


      // Possibly store data of plugins associated with cells
      signals.pre_refinement_store_user_data(triangulation);

      exchange_refinement_flags();

      triangulation.prepare_coarsening_and_refinement();
      bool any_flags_set = false;
      {
        for (const auto &cell:dof_handler.active_cell_iterators())
          {
            if (cell->refine_flag_set() || cell->coarsen_flag_set())
              {
                any_flags_set = true;
                break;
              }
          }
      }
      const bool mesh_changed = Utilities::MPI::max(any_flags_set?1:0,mpi_communicator) == 1 ? true : false;
      if (!mesh_changed)
        {
          pcout << "Skipping mesh refinement, because the mesh did not change.\n" << std::endl;
          return;
        }

      system_trans.prepare_for_coarsening_and_refinement(x_system);

      if (parameters.mesh_deformation_enabled)
        mesh_deformation_trans->prepare_for_coarsening_and_refinement(x_fs_system);

      triangulation.execute_coarsening_and_refinement ();
      if (MappingQCache<dim> *map = dynamic_cast<MappingQCache<dim>*>(&(*mapping)))
        map->initialize(MappingQGeneric<dim>(4), triangulation);
    } // leave the timed section

    setup_dofs ();

    {
      TimerOutput::Scope timer (computing_timer, "Refine mesh structure, part 2");

      LinearAlgebra::BlockVector distributed_system;
      LinearAlgebra::BlockVector old_distributed_system;
      LinearAlgebra::BlockVector distributed_mesh_velocity;

      distributed_system.reinit(introspection.index_sets.system_partitioning, mpi_communicator);
      old_distributed_system.reinit(introspection.index_sets.system_partitioning, mpi_communicator);
      if (parameters.mesh_deformation_enabled)
        distributed_mesh_velocity.reinit(introspection.index_sets.system_partitioning, mpi_communicator);

      std::vector<LinearAlgebra::BlockVector *> system_tmp
        = { &distributed_system, &old_distributed_system};

      if (parameters.mesh_deformation_enabled)
        system_tmp.push_back(&distributed_mesh_velocity);

      // transfer the data previously stored into the vectors indexed by
      // system_tmp. then ensure that the interpolated solution satisfies
      // hanging node constraints
      //
      // note that the 'constraints' variable contains hanging node constraints
      // and constraints from periodic boundary conditions (as well as from
      // zero and tangential velocity boundary conditions), but not from
      // non-homogeneous boundary conditions. the latter are added not in setup_dofs(),
      // which we call above, but added to 'current_constraints' in start_timestep(),
      // which we do not want to call here.
      //
      // however, what we have should be sufficient: we have everything that
      // is necessary to make the solution vectors *conforming* on the current
      // mesh.
      system_trans.interpolate (system_tmp);

      constraints.distribute (distributed_system);
      solution     = distributed_system;

      constraints.distribute (old_distributed_system);
      old_solution = old_distributed_system;

      // We need the current linearization point at the start of the new time step
      // when we set the boundary conditions for advected fields (to determine parts
      // of the boundary with outflow). Therefore, we here set it to the solution
      // vector, but it will be reinitialized the next time the equations are solved.
      current_linearization_point = distributed_system;

      // do the same as above, but for the mesh deformation solution
      if (parameters.mesh_deformation_enabled)
        {
          constraints.distribute (distributed_mesh_velocity);
          mesh_deformation->mesh_velocity = distributed_mesh_velocity;

          LinearAlgebra::Vector distributed_mesh_displacements,
                        distributed_old_mesh_displacements,
                        distributed_initial_topography;

          distributed_mesh_displacements.reinit(mesh_deformation->mesh_locally_owned,
                                                mpi_communicator);
          distributed_old_mesh_displacements.reinit(mesh_deformation->mesh_locally_owned,
                                                    mpi_communicator);
          distributed_initial_topography.reinit(mesh_deformation->mesh_locally_owned,
                                                mpi_communicator);

          std::vector<LinearAlgebra::Vector *> system_tmp
          = { &distributed_mesh_displacements,
              &distributed_old_mesh_displacements,
              &distributed_initial_topography
            };

          mesh_deformation_trans->interpolate (system_tmp);

          mesh_deformation->mesh_vertex_constraints.distribute (distributed_mesh_displacements);
          mesh_deformation->mesh_displacements = distributed_mesh_displacements;

          mesh_deformation->mesh_vertex_constraints.distribute (distributed_old_mesh_displacements);
          mesh_deformation->old_mesh_displacements = distributed_old_mesh_displacements;

          mesh_deformation->mesh_vertex_constraints.distribute (distributed_initial_topography);
          mesh_deformation->initial_topography = distributed_initial_topography;
        }

      // Possibly load data of plugins associated with cells
      signals.post_refinement_load_user_data(triangulation);

      // calculate global volume after refining mesh
      global_volume = GridTools::volume (triangulation, *mapping);
    }
  }



  template <int dim>
  void
  Simulator<dim>::
  solve_timestep ()
  {
    // start any scheme with an extrapolated value from the previous
    // two time steps if those are available
    initialize_current_linearization_point();

    // The mesh deformation scheme is currently not built to work inside a nonlinear solver.
    // We do the mesh deformation execution at the beginning of the timestep for a specific reason.
    // The time step size is calculated AFTER the whole solve_timestep() function.  If we call
    // mesh_deformation_execute() after the Stokes solve, it will be before we know what the appropriate
    // time step to take is, and we will timestep the boundary incorrectly.
    if (parameters.mesh_deformation_enabled)
      {
        mesh_deformation->execute ();

        // calculate global volume after deforming mesh
        global_volume = GridTools::volume (triangulation, *mapping);
        signals.post_mesh_deformation(*this);
      }

    // Compute the reactions of compositional fields and temperature in case of operator splitting.
    if (parameters.use_operator_splitting)
      compute_reactions ();

    try
      {
        switch (parameters.nonlinear_solver)
          {
            case NonlinearSolver::single_Advection_single_Stokes:
            {
              solve_single_advection_single_stokes();
              break;
            }

            case NonlinearSolver::no_Advection_iterated_Stokes:
            {
              solve_no_advection_iterated_stokes();
              break;
            }

            case NonlinearSolver::no_Advection_single_Stokes:
            {
              solve_no_advection_single_stokes();
              break;
            }

            case NonlinearSolver::iterated_Advection_and_Stokes:
            {
              solve_iterated_advection_and_stokes();
              break;
            }

            case NonlinearSolver::single_Advection_iterated_Stokes:
            {
              solve_single_advection_iterated_stokes();
              break;
            }

            case NonlinearSolver::no_Advection_iterated_defect_correction_Stokes:
            {
              solve_no_advection_iterated_defect_correction_stokes();
              break;
            }

            case NonlinearSolver::single_Advection_iterated_defect_correction_Stokes:
            {
              solve_single_advection_iterated_defect_correction_stokes();
              break;
            }

            case NonlinearSolver::iterated_Advection_and_defect_correction_Stokes:
            {
              solve_iterated_advection_and_defect_correction_stokes();
              break;
            }

            case NonlinearSolver::iterated_Advection_and_Newton_Stokes:
            {
              solve_iterated_advection_and_newton_stokes(/*use_newton_iterations =*/ true);
              break;
            }

            case NonlinearSolver::single_Advection_iterated_Newton_Stokes:
            {
              solve_single_advection_and_iterated_newton_stokes(/*use_newton_iterations =*/ true);
              break;
            }

            case NonlinearSolver::single_Advection_no_Stokes:
            {
              solve_single_advection_no_stokes();
              break;
            }

            case NonlinearSolver::first_timestep_only_single_Stokes:
            {
              solve_first_timestep_only_single_stokes();
              break;
            }

            case NonlinearSolver::no_Advection_no_Stokes:
            {
              solve_no_advection_no_stokes();
              break;
            }

            default:
              Assert (false, ExcNotImplemented());
          }
        pcout << std::endl;
      }
    catch (ExcNonlinearSolverNoConvergence &)
      {
        pcout << "\n"
              "   WARNING: The nonlinear solver in the current timestep failed to converge.\n"
              "   Acting according to the parameter 'Nonlinear solver failure strategy':"
              << std::endl;
        ++nonlinear_solver_failures;

        switch (parameters.nonlinear_solver_failure_strategy)
          {
            case Parameters<dim>::NonlinearSolverFailureStrategy::continue_with_next_timestep:
            {
              pcout << "   Continuing to the next timestep even though solution is not fully converged."
                    << std::endl
                    << std::endl;
              // swallow the exception and continue
              break;
            }
            case Parameters<dim>::NonlinearSolverFailureStrategy::cut_timestep_size:
            {
              if (timestep_number == 0)
                {
                  pcout << "   Error: Can not cut the timestep in step 0. Aborting."
                        << std::endl;
                  // Rethrow the current exception
                  throw;
                }
              else
                pcout << "   Trying to cut the current time step.\n"
                      << std::endl;

              time_stepping_manager.template get_matching_active_plugin<TimeStepping::RepeatOnNonlinearFail<dim>>()
              .nonlinear_solver_has_failed();
              break;
            }

            case Parameters<dim>::NonlinearSolverFailureStrategy::abort_program:
            {
              pcout << "   Aborting simulation as requested." << std::endl;
              // Rethrow the current exception
              throw;
            }

            default:
              AssertThrow(false, ExcNotImplemented());
          }
      }
  }


  /**
   * This is the main function of the program, containing the overall
   * logic which function is called when.
   */
  template <int dim>
  void Simulator<dim>::run ()
  {
    CitationInfo::print_info_block(pcout);

    unsigned int max_refinement_level = parameters.initial_global_refinement +
                                        parameters.initial_adaptive_refinement;
    pre_refinement_step = 0;

    // if we want to resume a computation from an earlier point
    // then reload it from a snapshot. otherwise do the basic
    // start-up
    if (parameters.resume_computation == true)
      {
        resume_from_snapshot();
        // we need to remove additional_refinement_times that are in the past
        // and adjust max_refinement_level which is not written to file
        while ((parameters.additional_refinement_times.size() > 0)
               &&
               (parameters.additional_refinement_times.front () < time+time_step))
          {
            ++max_refinement_level;
            parameters.additional_refinement_times
            .erase (parameters.additional_refinement_times.begin());
          }
      }
    else
      {
        time = parameters.start_time;

        // Instead of calling global_refine(n) we flag all cells for
        // refinement and then allow the mesh refinement plugins to unflag
        // the cells if desired. This procedure is repeated n times. If there
        // is no plugin that modifies the flags, it is equivalent to
        // refine_global(n).
        for (unsigned int n=0; n<parameters.initial_global_refinement; ++n)
          {
            for (const auto &cell : triangulation.active_cell_iterators())
              if (cell->is_locally_owned())
                cell->set_refine_flag ();

            mesh_refinement_manager.tag_additional_cells ();

            exchange_refinement_flags();

            triangulation.execute_coarsening_and_refinement();
            if (MappingQCache<dim> *map = dynamic_cast<MappingQCache<dim>*>(&(*mapping)))
              map->initialize(MappingQGeneric<dim>(4), triangulation);
          }

        setup_dofs();

        // calculate global volume after refining mesh
        global_volume = GridTools::volume (triangulation, *mapping);
      }

    // start the timer for periodic checkpoints after the setup above
    time_t last_checkpoint_time = std::time(nullptr);

  start_time_iteration:

    if (parameters.resume_computation == false)
      {
        TimerOutput::Scope timer (computing_timer, "Setup initial conditions");

        timestep_number           = 0;
        time_step = old_time_step = 0;

        if (! parameters.skip_setup_initial_conditions_on_initial_refinement
            ||
            ! (pre_refinement_step < parameters.initial_adaptive_refinement))
          {
            // Add topography to box models after all initial refinement
            // is completed.
            if (pre_refinement_step == parameters.initial_adaptive_refinement)
              signals.pre_set_initial_state (triangulation);

            set_initial_temperature_and_compositional_fields ();
            compute_initial_pressure_field ();

            signals.post_set_initial_state (*this);
          }
      }

    // Start the principal loop over time steps. At this point, everything
    // is completely initialized, so set that status as well
    simulator_is_past_initialization = true;
    do
      {
        // During pre-refinement, do not solve if we are asked to skip it:
        if (! (parameters.skip_solvers_on_initial_refinement
               && pre_refinement_step < parameters.initial_adaptive_refinement))
          {
            start_timestep ();

            // then do the core work: assemble systems and solve
            solve_timestep ();
          }

        // See if we have to start over with a new adaptive refinement cycle
        // at the beginning of the simulation. If so, set the
        // simulator_is_past_initialization variable back to false because we will
        // have to re-initialize some variables such as the size of vectors,
        // the initial state, etc.
        if (timestep_number == 0)
          {
            const bool initial_refinement_done = maybe_do_initial_refinement(max_refinement_level);
            if (initial_refinement_done)
              {
                simulator_is_past_initialization = false;
                goto start_time_iteration;
              }
          }

        // We're now definitely past the point where we need the initial
        // conditions objects, and we can release the pointers to these objects
        // that we have created in the constructor of this class. If some of the
        // other plugins created there still need access to these initial
        // conditions, they will have created their own shared pointers.
        initial_temperature_manager.reset();
        initial_composition_manager.reset();
#ifdef ASPECT_WITH_WORLD_BUILDER
        // The same applies to the world builder object:
        world_builder.reset();
#endif
        // Prepare the next time step:
        time_stepping_manager.update();

        const double new_time_step_size = time_stepping_manager.get_next_time_step_size();

        // if we postprocess nonlinear iterations, this function is called within
        // solve_timestep () in the individual solver schemes
        if (!time_stepping_manager.should_repeat_time_step()
            && !parameters.run_postprocessors_on_nonlinear_iterations)
          postprocess ();

        // if the time stepping manager tells us to refine the mesh,
        // we need to do this before going to the next time step
        if (time_stepping_manager.should_refine_mesh())
          {
            pcout << "Refining the mesh based on the time stepping manager ...\n" << std::endl;
            refine_mesh(max_refinement_level);
          }

        if (time_stepping_manager.should_repeat_time_step())
          {
            pcout << "Repeating the current time step based on the time stepping manager ...\n" << std::endl;

            if (mesh_deformation)
              mesh_deformation->mesh_displacements = mesh_deformation->old_mesh_displacements;

            // adjust time and time_step size:
            time = time - time_step + new_time_step_size;
            time_step = new_time_step_size;

            // Restore particles through stored copy of particle handler,
            // created in start_timestep(),
            // but only if this timestep is to be repeated.
            for (auto &particle_manager : particle_managers)
              particle_manager.restore_particles();

            continue; // repeat time step loop
          }

        if (!time_stepping_manager.should_refine_mesh())
          maybe_refine_mesh(new_time_step_size, max_refinement_level);

        // see if we want to write a timing summary
        maybe_write_timing_output();

        // update values for timestep, increment time step by one.
        advance_time(new_time_step_size);

        // Check whether to terminate the simulation:
        const bool should_terminate = time_stepping_manager.should_simulation_terminate_now();

        const bool write_checkpoint_due_to_termination
          = should_terminate && time_stepping_manager.need_checkpoint_on_terminate();

        const bool checkpoint_written = maybe_write_checkpoint(last_checkpoint_time,
                                                               write_checkpoint_due_to_termination);

        if (checkpoint_written)
          last_checkpoint_time = std::time(nullptr);

        if (should_terminate)
          break;
      }
    while (true);

    // we disable automatic summary printing so that it won't happen when
    // throwing an exception. Therefore, we have to do this manually here:
    computing_timer.print_summary ();

    if (nonlinear_solver_failures > 0)
      pcout << "\nWARNING: During this computation " << nonlinear_solver_failures << " nonlinear solver failures occurred!" << std::endl;

    pcout << "-- Total wallclock time elapsed including restarts: "
          << std::round(wall_timer.wall_time()+total_walltime_until_last_snapshot)
          << 's' << std::endl;

    CitationInfo::print_info_block (pcout);

    stokes_matrix_free.reset();
  }
}



// explicit instantiation of the functions we implement in this file
namespace aspect
{
#define INSTANTIATE(dim) \
  template class Simulator<dim>;

  ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
}
