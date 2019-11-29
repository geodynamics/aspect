/*
  Copyright (C) 2011 - 2019 by the authors of the ASPECT code.

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
#include <aspect/citation_info.h>
#include <aspect/postprocess/particles.h>

#ifdef ASPECT_WITH_WORLD_BUILDER
#include <world_builder/world.h>
#endif

#include <aspect/simulator/assemblers/interface.h>
#include <aspect/geometry_model/initial_topography_model/zero_topography.h>

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

#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/derivative_approximation.h>
#include <deal.II/numerics/vector_tools.h>

#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/distributed/grid_refinement.h>

#include <fstream>
#include <iostream>
#include <iomanip>
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
    std::vector<VariableDeclaration<dim> > construct_variables(const Parameters<dim> &parameters,
                                                               SimulatorSignals<dim> &signals,
                                                               std::unique_ptr<MeltHandler<dim> > &melt_handler)
    {
      std::vector<VariableDeclaration<dim> > variables
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
        return std_cxx14::make_unique<MappingQ<dim>>(4, true);
      if (dynamic_cast<const InitialTopographyModel::ZeroTopography<dim>*>(&initial_topography_model) != nullptr)
        return std_cxx14::make_unique<MappingCartesian<dim>>();

      return std_cxx14::make_unique<MappingQ1<dim>>();
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
    assemblers (std_cxx14::make_unique<Assemblers::Manager<dim>>()),
    parameters (prm, mpi_communicator_),
    melt_handler (parameters.include_melt_transport ?
                  std_cxx14::make_unique<MeltHandler<dim>>(prm) :
                  nullptr),
    newton_handler ((parameters.nonlinear_solver == NonlinearSolver::iterated_Advection_and_Newton_Stokes ||
                     parameters.nonlinear_solver == NonlinearSolver::single_Advection_iterated_Newton_Stokes) ?
                    std_cxx14::make_unique<NewtonHandler<dim>>() :
                    nullptr),
    post_signal_creation(
      std::bind (&internals::SimulatorSignals::call_connector_functions<dim>,
                 std::ref(signals))),
    volume_of_fluid_handler (parameters.volume_of_fluid_tracking_enabled ?
                             std_cxx14::make_unique<VolumeOfFluidHandler<dim>> (*this, prm) :
                             nullptr),
    introspection (construct_variables<dim>(parameters, signals, melt_handler), parameters),
    mpi_communicator (Utilities::MPI::duplicate_communicator (mpi_communicator_)),
    iostream_tee_device(std::cout, log_file_stream),
    iostream_tee_stream(iostream_tee_device),
    pcout (iostream_tee_stream,
           (Utilities::MPI::
            this_mpi_process(mpi_communicator)
            == 0)),

    computing_timer (mpi_communicator,
                     pcout,
                     TimerOutput::never,
                     TimerOutput::wall_times),
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
                   std_cxx14::make_unique<WorldBuilder::World>(parameters.world_builder_file) :
                   nullptr),
#endif
    boundary_heat_flux (BoundaryHeatFlux::create_boundary_heat_flux<dim>(prm)),
    particle_world(nullptr),
    time (numbers::signaling_nan<double>()),
    time_step (numbers::signaling_nan<double>()),
    old_time_step (numbers::signaling_nan<double>()),
    timestep_number (numbers::invalid_unsigned_int),
    nonlinear_iteration (numbers::invalid_unsigned_int),

    triangulation (mpi_communicator,
                   (parameters.stokes_solver_type == Parameters<dim>::StokesSolverType::block_gmg
                    ?
                    typename Triangulation<dim>::MeshSmoothing
                    (Triangulation<dim>::smoothing_on_refinement |
                     Triangulation<dim>::smoothing_on_coarsening |
                     Triangulation<dim>::limit_level_difference_at_vertices)
                    :
                    typename Triangulation<dim>::MeshSmoothing
                    (Triangulation<dim>::smoothing_on_refinement |
                     Triangulation<dim>::smoothing_on_coarsening))
                   ,
                   (parameters.stokes_solver_type == Parameters<dim>::StokesSolverType::block_gmg
                    ?
                    typename parallel::distributed::Triangulation<dim>::Settings
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
    assemble_newton_stokes_system ((parameters.nonlinear_solver == NonlinearSolver::iterated_Advection_and_Newton_Stokes ||
                                    parameters.nonlinear_solver == NonlinearSolver::single_Advection_iterated_Newton_Stokes) ? true : false),
    rebuild_stokes_preconditioner (true)
  {
    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
      {
        // only open the log file on processor 0, the other processors won't be
        // writing into the stream anyway
        log_file_stream.open((parameters.output_directory + "log.txt").c_str(),
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

    heating_model_manager.initialize_simulator (*this);
    heating_model_manager.parse_parameters (prm);

    if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(gravity_model.get()))
      sim->initialize_simulator (*this);
    gravity_model->parse_parameters (prm);
    gravity_model->initialize ();

    // Create the initial condition plugins
    initial_temperature_manager.initialize_simulator(*this);
    initial_temperature_manager.parse_parameters (prm);

    // Create the initial composition plugins
    initial_composition_manager.initialize_simulator(*this);
    initial_composition_manager.parse_parameters (prm);

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

    // Initialize the mesh deformation handler
    if (parameters.mesh_deformation_enabled)
      {
        // Allocate the MeshDeformationHandler object
        mesh_deformation = std_cxx14::make_unique<MeshDeformation::MeshDeformationHandler<dim>>(*this);
        mesh_deformation->initialize_simulator(*this);
        mesh_deformation->parse_parameters(prm);
      }

    // Initialize the melt handler
    if (parameters.include_melt_transport)
      {
        melt_handler->initialize_simulator (*this);
        melt_handler->initialize();
      }

    // If the solver type is a Newton type of solver, we need to set make sure
    // assemble_newton_stokes_system set to true.
    if (parameters.nonlinear_solver == NonlinearSolver::iterated_Advection_and_Newton_Stokes ||
        parameters.nonlinear_solver == NonlinearSolver::single_Advection_iterated_Newton_Stokes)
      {
        assemble_newton_stokes_system = true;
        newton_handler->initialize_simulator(*this);
        newton_handler->parameters.parse_parameters(prm);
      }

    if (parameters.stokes_solver_type == Parameters<dim>::StokesSolverType::block_gmg)
      {
        switch (parameters.stokes_velocity_degree)
          {
            case 2:
              stokes_matrix_free = std_cxx14::make_unique<StokesMatrixFreeHandlerImplementation<dim,2>>(*this, prm);
              break;
            case 3:
              stokes_matrix_free = std_cxx14::make_unique<StokesMatrixFreeHandlerImplementation<dim,3>>(*this, prm);
              break;
            default:
              AssertThrow(false, ExcMessage("The finite element degree for the Stokes system you selected is not supported yet."));
          }

      }

    postprocess_manager.initialize_simulator (*this);
    postprocess_manager.parse_parameters (prm);

    if (postprocess_manager.template has_matching_postprocessor<Postprocess::Particles<dim> >())
      {
        particle_world.reset(new Particle::World<dim>());
        if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(particle_world.get()))
          sim->initialize_simulator (*this);

        particle_world->parse_parameters(prm);
        particle_world->initialize();
      }

    mesh_refinement_manager.initialize_simulator (*this);
    mesh_refinement_manager.parse_parameters (prm);

    // VoF Must be initialized after mesh_refinement_manager, due to needing to check
    // for a mesh refinement strategy
    if (parameters.volume_of_fluid_tracking_enabled)
      {
        volume_of_fluid_handler->initialize (prm);
      }

    termination_manager.initialize_simulator (*this);
    termination_manager.parse_parameters (prm);

    lateral_averaging.initialize_simulator (*this);

    geometry_model->create_coarse_mesh (triangulation);
    global_Omega_diameter = GridTools::diameter (triangulation);

    for (std::map<types::boundary_id,std::pair<std::string,std::string> >::const_iterator
         p = parameters.prescribed_traction_boundary_indicators.begin();
         p != parameters.prescribed_traction_boundary_indicators.end();
         ++p)
      {
        BoundaryTraction::Interface<dim> *bv
          = BoundaryTraction::create_boundary_traction<dim>
            (p->second.second);
        boundary_traction[p->first].reset (bv);
        if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(bv))
          sim->initialize_simulator(*this);
        bv->parse_parameters (prm);
        bv->initialize ();
      }

    std::set<types::boundary_id> open_velocity_boundary_indicators
      = geometry_model->get_used_boundary_indicators();
    for (std::map<types::boundary_id,std::pair<std::string,std::vector<std::string> > >::const_iterator
         p = boundary_velocity_manager.get_active_boundary_velocity_names().begin();
         p != boundary_velocity_manager.get_active_boundary_velocity_names().end();
         ++p)
      open_velocity_boundary_indicators.erase (p->first);
    for (std::set<types::boundary_id>::const_iterator
         p = boundary_velocity_manager.get_zero_boundary_velocity_indicators().begin();
         p != boundary_velocity_manager.get_zero_boundary_velocity_indicators().end();
         ++p)
      open_velocity_boundary_indicators.erase (*p);
    for (std::set<types::boundary_id>::const_iterator
         p = boundary_velocity_manager.get_tangential_boundary_velocity_indicators().begin();
         p != boundary_velocity_manager.get_tangential_boundary_velocity_indicators().end();
         ++p)
      open_velocity_boundary_indicators.erase (*p);

    // We need to do the RHS compatibility modification, if the model is
    // compressible or compatible (in the case of melt transport), and
    // there is no open boundary to balance the pressure.
    do_pressure_rhs_compatibility_modification = ((material_model->is_compressible() && !parameters.include_melt_transport)
                                                  ||
                                                  (parameters.include_melt_transport && !material_model->is_compressible()))
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
        std::ofstream prm_out ((parameters.output_directory + "parameters.prm").c_str());
        AssertThrow (prm_out,
                     ExcMessage (std::string("Could not open file <") +
                                 parameters.output_directory + "parameters.prm>."));
        prm.print_parameters(prm_out, ParameterHandler::Text);

        std::ofstream tex_out ((parameters.output_directory + "parameters.tex").c_str());
        AssertThrow (tex_out,
                     ExcMessage (std::string("Could not open file <") +
                                 parameters.output_directory + "parameters.tex>."));
        prm.print_parameters(tex_out, ParameterHandler::LaTeX);

        std::ofstream json_out ((parameters.output_directory + "parameters.json").c_str());
        AssertThrow (json_out,
                     ExcMessage (std::string("Could not open file <") +
                                 parameters.output_directory + "parameters.json>."));
        prm.print_parameters(json_out, ParameterHandler::JSON);
      }

    // check that the boundary condition selection is consistent
    check_consistency_of_boundary_conditions();

    // check that the setup of equations, material models, and heating terms is consistent
    check_consistency_of_formulation();

    if (parameters.use_discontinuous_temperature_discretization || parameters.use_discontinuous_composition_discretization)
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
    particle_world.reset(nullptr);
    // wait if there is a thread that's still writing the statistics
    // object (set from the output_statistics() function)
    output_statistics_thread.join();

    // If an exception is being thrown (for example due to AssertThrow()), we
    // might end up here with currently active timing sections. The destructor
    // of TimerOutput does MPI communication, which can lead to deadlocks,
    // hangs, or confusing MPI error messages. To avoid this, we can call
    // reset() to remove all open sections. In a normal run, we won't have any
    // active sessions, so this won't hurt to do:
    computing_timer.reset();
  }


  namespace
  {
    /**
     * Conversion object where one can provide a function that returns
     * a tensor for the velocity at a given point and it returns something
     * that matches the dealii::Function interface with a number of output
     * components equal to the number of components of the finite element
     * in use.
     */
    template <int dim>
    class VectorFunctionFromVelocityFunctionObject : public Function<dim>
    {
      public:
        /**
         * Given a function object that takes a Point and returns a Tensor<1,dim>,
         * convert this into an object that matches the Function@<dim@>
         * interface.
         *
         * @param n_components total number of components of the finite element system.
         * @param function_object The function that will form one component
         *     of the resulting Function object.
         */
        VectorFunctionFromVelocityFunctionObject (const unsigned int n_components,
                                                  const std::function<Tensor<1,dim> (const Point<dim> &)> &function_object);

        /**
         * Return the value of the
         * function at the given
         * point. Returns the value the
         * function given to the constructor
         * produces for this point.
         */
        double value (const Point<dim>   &p,
                      const unsigned int  component = 0) const override;

        /**
         * Return all components of a
         * vector-valued function at a
         * given point.
         *
         * <tt>values</tt> shall have the right
         * size beforehand,
         * i.e. #n_components.
         */
        void vector_value (const Point<dim>   &p,
                           Vector<double>     &values) const override;

      private:
        /**
         * The function object which we call when this class's value() or
         * value_list() functions are called.
         **/
        const std::function<Tensor<1,dim> (const Point<dim> &)> function_object;
    };


    template <int dim>
    VectorFunctionFromVelocityFunctionObject<dim>::
    VectorFunctionFromVelocityFunctionObject
    (const unsigned int n_components,
     const std::function<Tensor<1,dim> (const Point<dim> &)> &function_object)
      :
      Function<dim>(n_components),
      function_object (function_object)
    {
    }



    template <int dim>
    double
    VectorFunctionFromVelocityFunctionObject<dim>::value (const Point<dim> &p,
                                                          const unsigned int component) const
    {
      Assert (component < this->n_components,
              ExcIndexRange (component, 0, this->n_components));

      if (component < dim)
        {
          const Tensor<1,dim> v = function_object(p);
          return v[component];
        }
      else
        return 0;
    }



    template <int dim>
    void
    VectorFunctionFromVelocityFunctionObject<dim>::
    vector_value (const Point<dim>   &p,
                  Vector<double>     &values) const
    {
      AssertDimension(values.size(), this->n_components);

      // set everything to zero, and then the right components to their correct values
      values = 0;

      const Tensor<1,dim> v = function_object(p);
      for (unsigned int d=0; d<dim; ++d)
        values(d) = v[d];
    }
  }


  template <int dim>
  void
  Simulator<dim>::
  start_timestep ()
  {
    // first produce some output for the screen to show where we are
    if (parameters.convert_to_years == true)
      pcout << "*** Timestep " << timestep_number
            << ":  t=" << time/year_in_seconds
            << " years"
            << std::endl;
    else
      pcout << "*** Timestep " << timestep_number
            << ":  t=" << time
            << " seconds"
            << std::endl;

    nonlinear_iteration = 0;

    // Re-compute the pressure scaling factor. In some sense, it would be nice
    // if we did this not just once per time step, but once for each solve --
    // i.e., multiple times per time step if we iterate out the nonlinearity
    // during a Newton or Picard iteration. But that's more work,
    // and the function would need to be called in more different places.
    // Unless we have evidence that that's necessary, let's assume that
    // the reference viscosity does not change too much between nonlinear
    // iterations and that it's ok to update it only once per time step.
    //
    // The call to this function must precede the one to the computation
    // of constraints because some kinds of constraints require scaling
    // pressure degrees of freedom to a size adjusted by the pressure_scaling
    // factor.
    compute_pressure_scaling_factor();

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
    material_model->update();
    gravity_model->update();
    heating_model_manager.update();
    adiabatic_conditions->update();
    mesh_refinement_manager.update();
    if (parameters.mesh_deformation_enabled)
      mesh_deformation->update();

    if (prescribed_stokes_solution.get())
      prescribed_stokes_solution->update();

    // do the same for the traction boundary conditions and other things
    // that end up in the bilinear form. we update those that end up in
    // the constraints object when calling compute_current_constraints()
    // above
    for (auto &p : boundary_traction)
      p.second->update ();
  }



  template <int dim>
  void
  Simulator<dim>::
  compute_current_constraints ()
  {
    // We put the constraints we compute into a separate ConstraintMatrix so we can check
    // if the set of constraints has changed. If it did, we need to update the sparsity patterns.
    ConstraintMatrix new_current_constraints;
    new_current_constraints.clear ();
    new_current_constraints.reinit (introspection.index_sets.system_relevant_set);
    new_current_constraints.merge (constraints);
    compute_current_velocity_boundary_constraints(new_current_constraints);

    // If there is a fixed boundary temperature or heat flux,
    // update the temperature boundary condition.
    boundary_temperature_manager.update();
    boundary_heat_flux->update();

    // If we do not want to prescribe Dirichlet boundary conditions on outflow boundaries,
    // we update the boundary indicators of all faces that belong to ouflow boundaries
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
        for (std::set<types::boundary_id>::const_iterator
             p = boundary_temperature_manager.get_fixed_temperature_boundary_indicators().begin();
             p != boundary_temperature_manager.get_fixed_temperature_boundary_indicators().end(); ++p)
          {
            auto lambda = [&] (const Point<dim> &x) -> double
            {
              return boundary_temperature_manager.boundary_temperature(*p, x);
            };

            VectorFunctionFromScalarFunctionObject<dim> vector_function_object(
              lambda,
              introspection.component_masks.temperature.first_selected_component(),
              introspection.n_components);

            VectorTools::interpolate_boundary_values (*mapping,
                                                      dof_handler,
                                                      *p,
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

    // now do the same for the composition variable:
    if (!parameters.use_discontinuous_composition_discretization)
      {
        // obtain the boundary indicators that belong to Dirichlet-type
        // composition boundary conditions and interpolate the composition
        // there
        for (unsigned int c=0; c<introspection.n_compositional_fields; ++c)
          for (std::set<types::boundary_id>::const_iterator
               p = boundary_composition_manager.get_fixed_composition_boundary_indicators().begin();
               p != boundary_composition_manager.get_fixed_composition_boundary_indicators().end(); ++p)
            {
              auto lambda = [&] (const Point<dim> &x) -> double
              {
                return boundary_composition_manager.boundary_composition(*p, x, c);
              };

              VectorFunctionFromScalarFunctionObject<dim> vector_function_object(
                lambda,
                introspection.component_masks.compositional_fields[c].first_selected_component(),
                introspection.n_components);

              VectorTools::interpolate_boundary_values (*mapping,
                                                        dof_handler,
                                                        *p,
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
    const bool any_constrained_dofs_set_changed = Utilities::MPI::sum(constraints_changed ? 1 : 0,
                                                                      mpi_communicator) > 0;
    if (any_constrained_dofs_set_changed)
      rebuild_sparsity_and_matrices = true;

    current_constraints.copy_from(new_current_constraints);

    // TODO: We should use current_constraints.is_consistent_in_parallel()
    // here to assert that our constraints are consistent between
    // processors. This got removed in
    // https://github.com/geodynamics/aspect/pull/3282 because it triggered in
    // normal computations due to small floating point differences. See
    // https://github.com/geodynamics/aspect/issues/3248 for the discussion.
  }



  template <int dim>
  void
  Simulator<dim>::
  setup_system_matrix (const std::vector<IndexSet> &system_partitioning)
  {
    system_matrix.clear ();

    bool have_fem_compositional_field = false;
    for (unsigned int c=0; c<introspection.n_compositional_fields; ++c)
      {
        const AdvectionField adv_field (AdvectionField::composition(c));
        if (adv_field.advection_method(introspection)==Parameters<dim>::AdvectionFieldMethod::fem_field
            || adv_field.advection_method(introspection)==Parameters<dim>::AdvectionFieldMethod::fem_melt_field)
          {
            have_fem_compositional_field = true;
            break;
          }
      }

    Table<2,DoFTools::Coupling> coupling (introspection.n_components,
                                          introspection.n_components);
    coupling.fill (DoFTools::none);

    // determine which blocks should be fillable in the matrix.
    // note:
    // - all velocities couple with all velocities
    // - pressure couples with all velocities and the other way
    //   around
    // - temperature only couples with itself
    // - compositional fields only couple with themselves
    // - additionally, in models with melt transport fluid pressure
    //   and compaction pressures couple with themselves
    {
      const typename Introspection<dim>::ComponentIndices &x
        = introspection.component_indices;

      for (unsigned int c=0; c<dim; ++c)
        for (unsigned int d=0; d<dim; ++d)
          coupling[x.velocities[c]][x.velocities[d]] = DoFTools::always;

      if (parameters.include_melt_transport)
        {
          for (unsigned int d=0; d<dim; ++d)
            {
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
          for (unsigned int d=0; d<dim; ++d)
            {
              coupling[x.velocities[d]][x.pressure] = DoFTools::always;
              coupling[x.pressure][x.velocities[d]] = DoFTools::always;
            }
        }
      // Do not allocate a temperature matrix if no temperature
      // solves are going to be performed.
      if (!(parameters.nonlinear_solver == NonlinearSolver::Kind::no_Advection_iterated_Stokes
            ||
            parameters.nonlinear_solver == NonlinearSolver::Kind::no_Advection_no_Stokes
            ||
            parameters.nonlinear_solver == NonlinearSolver::Kind::first_timestep_only_single_Stokes))
        coupling[x.temperature][x.temperature] = DoFTools::always;

      // For equal-order interpolation, we need a stabilization term
      // in the bottom right of Stokes matrix. Make sure we have the
      // necessary entries.
      if (parameters.use_equal_order_interpolation_for_stokes == true)
        coupling[x.pressure][x.pressure] = DoFTools::always;

      // If we have at least one compositional field that is a FEM field, we
      // create a matrix block in the first compositional block. Its sparsity
      // pattern will later be used to allocate composition matrices as
      // needed.  All other matrix blocks are left empty here.
      if (have_fem_compositional_field)
        coupling[x.compositional_fields[0]][x.compositional_fields[0]] = DoFTools::always;

      // If we are using VolumeOfFluid interface tracking, create a matrix block in the
      // field corresponding to the volume fraction.
      if (parameters.volume_of_fluid_tracking_enabled)
        {
          const unsigned int volume_of_fluid_block = volume_of_fluid_handler->field_struct_for_field_index(0)
                                                     .volume_fraction.first_component_index;
          coupling[volume_of_fluid_block][volume_of_fluid_block] = DoFTools::always;
        }
      if (stokes_matrix_free)
        {
          // do not allocate memory for the Stokes matrix:
          Assert(!parameters.include_melt_transport, ExcNotImplemented());
          for (unsigned int c=0; c<dim; ++c)
            for (unsigned int d=0; d<dim; ++d)
              coupling[x.velocities[c]][x.velocities[d]] = DoFTools::none;
          for (unsigned int d=0; d<dim; ++d)
            {
              coupling[x.velocities[d]][x.pressure] = DoFTools::none;
              coupling[x.pressure][x.velocities[d]] = DoFTools::none;
            }
        }
    }

    LinearAlgebra::BlockDynamicSparsityPattern sp;
#ifdef ASPECT_USE_PETSC
    sp.reinit (introspection.index_sets.system_relevant_partitioning);
#else
    sp.reinit (system_partitioning,
               system_partitioning,
               introspection.index_sets.system_relevant_partitioning,
               mpi_communicator);
#endif

    if ((parameters.use_discontinuous_temperature_discretization) ||
        (parameters.use_discontinuous_composition_discretization) ||
        (parameters.volume_of_fluid_tracking_enabled))
      {
        Table<2,DoFTools::Coupling> face_coupling (introspection.n_components,
                                                   introspection.n_components);
        face_coupling.fill (DoFTools::none);

        const typename Introspection<dim>::ComponentIndices &x
          = introspection.component_indices;
        if (parameters.use_discontinuous_temperature_discretization)
          face_coupling[x.temperature][x.temperature] = DoFTools::always;

        // Only allocate composition 0 matrix if needed. Same as the non-DG case (see above)
        if (parameters.use_discontinuous_composition_discretization && have_fem_compositional_field)
          face_coupling[x.compositional_fields[0]][x.compositional_fields[0]] = DoFTools::always;

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
      }
    else
      DoFTools::make_sparsity_pattern (dof_handler,
                                       coupling, sp,
                                       current_constraints, false,
                                       Utilities::MPI::
                                       this_mpi_process(mpi_communicator));

#ifdef ASPECT_USE_PETSC
    SparsityTools::distribute_sparsity_pattern(sp,
                                               dof_handler.locally_owned_dofs_per_processor(),
                                               mpi_communicator, introspection.index_sets.system_relevant_set);

    sp.compress();

    system_matrix.reinit (system_partitioning, system_partitioning, sp, mpi_communicator);
#else
    sp.compress();

    // We only allocate a composition matrix block for composition 0 (see
    // above). But even though we specify a coupling of DoFTools::none for the
    // other composition blocks, entries for constrained entries for boundary
    // conditions and hanging nodes are being created by
    // make_sparsity_pattern. These are unnecessary, so we remove those
    // entries here.
    for (unsigned int c=1; c<introspection.n_compositional_fields; ++c)
      {
        const unsigned int block_idx = introspection.block_indices.compositional_fields[c];
        // TODO: using clear() would be nice here but clear() also resets the
        // size, so just reinit():
        sp.block(block_idx, block_idx).reinit(sp.block(block_idx, block_idx).locally_owned_range_indices(),sp.block(block_idx, block_idx).locally_owned_domain_indices());
        sp.block(block_idx, block_idx).compress();
      }

    system_matrix.reinit (sp);
#endif
  }



  template <int dim>
  void Simulator<dim>::
  setup_system_preconditioner (const std::vector<IndexSet> &system_partitioning)
  {
    Amg_preconditioner.reset ();
    Mp_preconditioner.reset ();
    system_preconditioner_matrix.clear ();

    // The preconditioner matrix is only used for the Stokes block (velocity and Schur complement) and only needed if we actually solve iteratively and matrix-based
    if (parameters.stokes_solver_type == Parameters<dim>::StokesSolverType::block_gmg)
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

    // velocity-velocity block (only block diagonal):
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

#ifdef ASPECT_USE_PETSC
    sp.reinit (introspection.index_sets.system_relevant_partitioning);
#else
    sp.reinit (system_partitioning,
               system_partitioning,
               introspection.index_sets.system_relevant_partitioning,
               mpi_communicator);
#endif

    DoFTools::make_sparsity_pattern (dof_handler,
                                     coupling, sp,
                                     current_constraints, false,
                                     Utilities::MPI::
                                     this_mpi_process(mpi_communicator));

#ifdef ASPECT_USE_PETSC
    SparsityTools::distribute_sparsity_pattern(sp,
                                               dof_handler.locally_owned_dofs_per_processor(),
                                               mpi_communicator, introspection.index_sets.system_relevant_set);

    sp.compress();

    system_preconditioner_matrix.reinit (system_partitioning, system_partitioning, sp, mpi_communicator);
#else
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
      sp.block(block_idx, block_idx).reinit(sp.block(block_idx, block_idx).locally_owned_range_indices(),sp.block(block_idx, block_idx).locally_owned_domain_indices());
      sp.block(block_idx, block_idx).compress();
    }
    // compositions:
    for (unsigned int c=0; c<introspection.n_compositional_fields; ++c)
      {
        const unsigned int block_idx = introspection.block_indices.compositional_fields[c];
        // TODO: using clear() would be nice here but clear() also resets the
        // size, so just reinit():
        sp.block(block_idx, block_idx).reinit(sp.block(block_idx, block_idx).locally_owned_range_indices(),sp.block(block_idx, block_idx).locally_owned_domain_indices());
        sp.block(block_idx, block_idx).compress();
      }

    system_preconditioner_matrix.reinit (sp);
#endif
  }


  template <int dim>
  void Simulator<dim>::compute_initial_velocity_boundary_constraints (ConstraintMatrix &constraints)
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
      for (std::set<types::boundary_id>::const_iterator
           p = boundary_velocity_manager.get_zero_boundary_velocity_indicators().begin();
           p != boundary_velocity_manager.get_zero_boundary_velocity_indicators().end(); ++p)
        VectorTools::interpolate_boundary_values (*mapping,
                                                  dof_handler,
                                                  *p,
                                                  ZeroFunction<dim>(introspection.n_components),
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
  void Simulator<dim>::compute_current_velocity_boundary_constraints (ConstraintMatrix &constraints)
  {
    // set the current time and do the interpolation
    // for the prescribed velocity fields
    boundary_velocity_manager.update();
    for (typename std::map<types::boundary_id,std::pair<std::string, std::vector<std::string> > >::const_iterator
         p = boundary_velocity_manager.get_active_boundary_velocity_names().begin();
         p != boundary_velocity_manager.get_active_boundary_velocity_names().end(); ++p)
      {
        VectorFunctionFromVelocityFunctionObject<dim> vel
        (introspection.n_components,
         [&] (const Point<dim> &x) -> Tensor<1,dim>
        {
          return boundary_velocity_manager.boundary_velocity(p->first, x);
        });

        // here we create a mask for interpolate_boundary_values out of the 'selector'
        std::vector<bool> mask(introspection.component_masks.velocities.size(), false);
        const std::string &comp = p->second.first;

        if (comp.length()>0)
          {
            for (std::string::const_iterator direction=comp.begin(); direction!=comp.end(); ++direction)
              {
                switch (*direction)
                  {
                    case 'x':
                      mask[introspection.component_indices.velocities[0]] = true;
                      break;
                    case 'y':
                      mask[introspection.component_indices.velocities[1]] = true;
                      break;
                    case 'z':
                      // we must be in 3d, or 'z' should never have gotten through
                      Assert (dim==3, ExcInternalError());
                      if (dim==3)
                        mask[introspection.component_indices.velocities[dim-1]] = true;
                      break;
                    default:
                      Assert (false, ExcInternalError());
                  }
              }
          }
        else
          {
            // no mask given -- take all velocities
            for (unsigned int i=0; i<introspection.component_masks.velocities.size(); ++i)
              mask[i]=introspection.component_masks.velocities[i];
          }

        if (!assemble_newton_stokes_system || (assemble_newton_stokes_system && nonlinear_iteration == 0))
          {
            VectorTools::interpolate_boundary_values (*mapping,
                                                      dof_handler,
                                                      p->first,
                                                      vel,
                                                      constraints,
                                                      mask);
          }
        else
          {
            VectorTools::interpolate_boundary_values (*mapping,
                                                      dof_handler,
                                                      p->first,
                                                      ZeroFunction<dim>(introspection.n_components),
                                                      constraints,
                                                      mask);
          }
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
    constraints.clear();
    constraints.reinit(introspection.index_sets.system_relevant_set);

    // Set up the constraints for periodic boundary conditions:

    // Note: this has to happen _before_ we do hanging node constraints,
    // because inconsistent contraints could be generated in parallel otherwise.
    {
      typedef std::set< std::pair< std::pair< types::boundary_id, types::boundary_id>, unsigned int> >
      periodic_boundary_set;
      periodic_boundary_set pbs = geometry_model->get_periodic_boundary_pairs();

      for (periodic_boundary_set::iterator p = pbs.begin(); p != pbs.end(); ++p)
        {
          DoFTools::make_periodicity_constraints(dof_handler,
                                                 (*p).first.first,  // first boundary id
                                                 (*p).first.second, // second boundary id
                                                 (*p).second,       // cartesian direction for translational symmetry
                                                 constraints);
        }
    }

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
    DoFTools::count_dofs_per_block (dof_handler,
                                    introspection.system_dofs_per_block,
                                    introspection.get_components_to_blocks());
    {
      IndexSet system_index_set = dof_handler.locally_owned_dofs();
      aspect::Utilities::split_by_block (introspection.system_dofs_per_block,
                                         system_index_set,
                                         introspection.index_sets.system_partitioning);

      DoFTools::extract_locally_relevant_dofs (dof_handler,
                                               introspection.index_sets.system_relevant_set);
      aspect::Utilities::split_by_block (introspection.system_dofs_per_block,
                                         introspection.index_sets.system_relevant_set,
                                         introspection.index_sets.system_relevant_partitioning);

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
    std::list<std::pair<std::string,std::string> >
    output_list = postprocess_manager.execute (statistics);

    // if we are on processor zero, print to screen
    // whatever the postprocessors have generated
    if (Utilities::MPI::this_mpi_process(mpi_communicator)==0)
      {
        // determine the width of the first column of text so that
        // everything gets nicely aligned; then output everything
        {
          unsigned int width = 0;
          for (std::list<std::pair<std::string,std::string> >::const_iterator
               p = output_list.begin();
               p != output_list.end(); ++p)
            width = std::max<unsigned int> (width, p->first.size());

          for (std::list<std::pair<std::string,std::string> >::const_iterator
               p = output_list.begin();
               p != output_list.end(); ++p)
            pcout << "     "
                  << std::left
                  << std::setw(width)
                  << p->first
                  << " "
                  << p->second
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

    std::unique_ptr<parallel::distributed::SolutionTransfer<dim,LinearAlgebra::Vector> >
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

      std::vector<const LinearAlgebra::BlockVector *> x_system (2);
      x_system[0] = &solution;
      x_system[1] = &old_solution;

      if (parameters.mesh_deformation_enabled)
        x_system.push_back( &mesh_deformation->mesh_velocity );

      std::vector<const LinearAlgebra::Vector *> x_fs_system (1);

      if (parameters.mesh_deformation_enabled)
        {
          x_fs_system[0] = &mesh_deformation->mesh_displacements;
          mesh_deformation_trans
            = std_cxx14::make_unique<parallel::distributed::SolutionTransfer<dim,LinearAlgebra::Vector>>
              (mesh_deformation->mesh_deformation_dof_handler);
        }


      // Possibly store data of plugins associated with cells
      signals.pre_refinement_store_user_data(triangulation);

      triangulation.prepare_coarsening_and_refinement();
      system_trans.prepare_for_coarsening_and_refinement(x_system);

      if (parameters.mesh_deformation_enabled)
        mesh_deformation_trans->prepare_for_coarsening_and_refinement(x_fs_system);

      triangulation.execute_coarsening_and_refinement ();
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

      std::vector<LinearAlgebra::BlockVector *> system_tmp (2);
      system_tmp[0] = &distributed_system;
      system_tmp[1] = &old_distributed_system;

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

      // do the same as above, but for the mesh deformation solution
      if (parameters.mesh_deformation_enabled)
        {
          constraints.distribute (distributed_mesh_velocity);
          mesh_deformation->mesh_velocity = distributed_mesh_velocity;

          LinearAlgebra::Vector distributed_mesh_displacements;

          distributed_mesh_displacements.reinit(mesh_deformation->mesh_locally_owned,
                                                mpi_communicator);

          std::vector<LinearAlgebra::Vector *> system_tmp (1);
          system_tmp[0] = &distributed_mesh_displacements;

          mesh_deformation_trans->interpolate (system_tmp);
          mesh_deformation->mesh_vertex_constraints.distribute (distributed_mesh_displacements);
          mesh_deformation->mesh_displacements = distributed_mesh_displacements;
        }

      // Possibly load data of plugins associated with cells
      signals.post_refinement_load_user_data(triangulation);

      // calculate global volume after displacing mesh (if we have, in fact, displaced it)
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
      mesh_deformation->execute ();

    // Compute the reactions of compositional fields and temperature in case of operator splitting.
    if (parameters.use_operator_splitting)
      compute_reactions ();

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

        case NonlinearSolver::iterated_Advection_and_Newton_Stokes:
        {
          solve_iterated_advection_and_newton_stokes();
          break;
        }

        case NonlinearSolver::single_Advection_iterated_Newton_Stokes:
        {
          solve_single_advection_iterated_newton_stokes();
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

    if (particle_world.get() != nullptr)
      {
        // Do not advect the particles in the initial refinement stage
        const bool in_initial_refinement = (timestep_number == 0)
                                           && (pre_refinement_step < parameters.initial_adaptive_refinement);
        if (!in_initial_refinement)
          // Advance the particles in the world to the current time
          particle_world->advance_timestep();

        if (particle_world->get_property_manager().need_update() == Particle::Property::update_output_step)
          particle_world->update_particles();
      }
    pcout << std::endl;
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
              cell->set_refine_flag ();

            mesh_refinement_manager.tag_additional_cells ();
            triangulation.execute_coarsening_and_refinement();
          }

        setup_dofs();

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
        // Only solve if we are not in pre-refinement, or we do not want to skip
        // solving in pre-refinement.
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

        // if we postprocess nonlinear iterations, this function is called within
        // solve_timestep () in the individual solver schemes
        if (!parameters.run_postprocessors_on_nonlinear_iterations)
          postprocess ();

        // get new time step size
        const double new_time_step = compute_time_step();

        // see if we want to refine the mesh
        maybe_refine_mesh(new_time_step,max_refinement_level);

        // see if we want to write a timing summary
        maybe_write_timing_output();

        // update values for timestep, increment time step by one.
        old_time_step = time_step;
        time_step = new_time_step;
        time += time_step;
        ++timestep_number;

        // prepare for the next time step by shifting solution vectors
        // by one time step. In timestep 0 (just increased in the
        // line above) initialize both old_solution
        // and old_old_solution with the currently computed solution.
        if (timestep_number == 1)
          {
            old_old_solution      = solution;
            old_solution          = solution;
          }
        else
          {
            old_old_solution      = old_solution;
            old_solution          = solution;
          }

        // check whether to terminate the simulation. the
        // first part of the pair indicates whether to terminate
        // the execution; the second indicates whether to do one
        // more checkpoint
        const std::pair<bool,bool> termination = termination_manager.execute();

        const bool checkpoint_written = maybe_write_checkpoint(last_checkpoint_time,termination);
        if (checkpoint_written)
          last_checkpoint_time = std::time(nullptr);

        // see if we want to terminate
        if (termination.first)
          break;
      }
    while (true);

    // we disable automatic summary printing so that it won't happen when
    // throwing an exception. Therefore, we have to do this manually here:
    computing_timer.print_summary ();

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
}
