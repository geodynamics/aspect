/* $Id$ */
/* Author: Martin Kronbichler, Uppsala University,
           Wolfgang Bangerth, Texas A&M University,
     Timo Heister, University of Goettingen, 2008-2011 */
/*                                                                */
/*    Copyright (C) 2008, 2009, 2010, 2011 by the deal.II authors */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */

#include <aspect/simulator.h>
#include <aspect/equation_data.h>
#include <aspect/global.h>

#include <deal.II/base/index_set.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vectors.h>

#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/distributed/grid_refinement.h>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <locale>
#include <string>


using namespace dealii;

// In the following namespace, we define the
// various pieces of equation data. All of
// these are exhaustively discussed in the
// description of the testcase in the
// introduction:
namespace EquationData
{
  double kappa                 = 1e-6;
}



namespace aspect
{
  /**
   * Constructor. Initialize all member variables.
   **/
  template <int dim>
  Simulator<dim>::Simulator (ParameterHandler &prm)
    :
    parameters (prm),
    pcout (std::cout,
           (Utilities::MPI::
            this_mpi_process(MPI_COMM_WORLD)
            == 0)),

    geometry_model (GeometryModel::create_geometry_model<dim>(prm)),
    material_model (MaterialModel::create_material_model<dim>(prm)),
    gravity_model (GravityModel::create_gravity_model<dim>(prm)),
    boundary_temperature (BoundaryTemperature::create_boundary_temperature<dim>(prm)),
    initial_conditions (InitialConditions::create_initial_conditions (prm,
                                                                      *geometry_model,
                                                                      *boundary_temperature,
                                                                      *adiabatic_conditions)),

    triangulation (MPI_COMM_WORLD,
                   typename Triangulation<dim>::MeshSmoothing
                   (Triangulation<dim>::smoothing_on_refinement |
                    Triangulation<dim>::smoothing_on_coarsening),
                   parallel::distributed::Triangulation<dim>::mesh_reconstruction_after_repartitioning),

    mapping (4),

    stokes_fe (FE_Q<dim>(parameters.stokes_velocity_degree),
               dim,
               (parameters.use_locally_conservative_discretization
                ?
                static_cast<const FiniteElement<dim> &>
                (FE_DGP<dim>(parameters.stokes_velocity_degree-1))
                :
                static_cast<const FiniteElement<dim> &>
                (FE_Q<dim>(parameters.stokes_velocity_degree-1))),
               1),

    stokes_dof_handler (triangulation),

    temperature_fe (parameters.temperature_degree),
    temperature_dof_handler (triangulation),

    time (0),
    time_step (0),
    old_time_step (0),
    timestep_number (0),
    rebuild_stokes_matrix (true),
    rebuild_stokes_preconditioner (true),

    computing_timer (pcout, TimerOutput::summary,
                     TimerOutput::wall_times)
  {
    // continue with initializing members that can't be initialized for one reason
    // or another in the member initializer list above
    postprocess_manager.parse_parameters (prm);
    postprocess_manager.initialize (*this);

    geometry_model->create_coarse_mesh (triangulation);
    global_Omega_diameter = GridTools::diameter (triangulation);

    adiabatic_conditions.reset (new AdiabaticConditions<dim>(*geometry_model,
                                                             *gravity_model,
                                                             *material_model));

    pressure_scaling = material_model->reference_viscosity() / geometry_model->length_scale();

    // make sure that we don't have to fill every column of the statistics
    // object in each time step.
    statistics.set_auto_fill_mode(true);

    // finally produce a record of the run-time parameters by writing
    // the currently used values into a file
    std::ofstream prm_out ((parameters.output_directory + "parameters.prm").c_str());
    AssertThrow (prm_out, ExcIO());
    prm.print_parameters(prm_out, ParameterHandler::Text);
  }



  template <int dim>
  void Simulator<dim>::
  setup_stokes_matrix (const std::vector<IndexSet> &stokes_partitioning)
  {
    stokes_matrix.clear ();

    TrilinosWrappers::BlockSparsityPattern sp (stokes_partitioning,
                                               MPI_COMM_WORLD);

    Table<2,DoFTools::Coupling> coupling (dim+1, dim+1);

    for (unsigned int c=0; c<dim+1; ++c)
      for (unsigned int d=0; d<dim+1; ++d)
        if (! ((c==dim) && (d==dim)))
          coupling[c][d] = DoFTools::always;
        else
          coupling[c][d] = DoFTools::none;

    DoFTools::make_sparsity_pattern (stokes_dof_handler,
                                     coupling, sp,
                                     stokes_constraints, false,
                                     Utilities::MPI::
                                     this_mpi_process(MPI_COMM_WORLD));
    sp.compress();

    stokes_matrix.reinit (sp);
  }



  template <int dim>
  void Simulator<dim>::
  setup_stokes_preconditioner (const std::vector<IndexSet> &stokes_partitioning)
  {
    Amg_preconditioner.reset ();
    Mp_preconditioner.reset ();

    stokes_preconditioner_matrix.clear ();

    TrilinosWrappers::BlockSparsityPattern sp (stokes_partitioning,
                                               MPI_COMM_WORLD);

    Table<2,DoFTools::Coupling> coupling (dim+1, dim+1);
    for (unsigned int c=0; c<dim+1; ++c)
      for (unsigned int d=0; d<dim+1; ++d)
        if (c == d)
          coupling[c][d] = DoFTools::always;
        else
          coupling[c][d] = DoFTools::none;

    DoFTools::make_sparsity_pattern (stokes_dof_handler,
                                     coupling, sp,
                                     stokes_constraints, false,
                                     Utilities::MPI::
                                     this_mpi_process(MPI_COMM_WORLD));
    sp.compress();

    stokes_preconditioner_matrix.reinit (sp);
  }


  template <int dim>
  void Simulator<dim>::
  setup_temperature_matrix (const IndexSet &temperature_partitioner)
  {
    T_preconditioner.reset ();
    temperature_matrix.clear ();

    TrilinosWrappers::SparsityPattern sp (temperature_partitioner,
                                          MPI_COMM_WORLD);
    DoFTools::make_sparsity_pattern (temperature_dof_handler, sp,
                                     temperature_constraints, false,
                                     Utilities::MPI::
                                     this_mpi_process(MPI_COMM_WORLD));
    sp.compress();

    temperature_matrix.reinit (sp);
  }



  template <int dim>
  void Simulator<dim>::setup_dofs ()
  {
    computing_timer.enter_section("Setup dof systems");

    stokes_dof_handler.distribute_dofs (stokes_fe);

    // Renumber the DoFs hierarchical so that we get the
    // same numbering if we resume the computation. This
    // is because the numbering depends on the order the
    // cells are created.
    DoFRenumbering::hierarchical (stokes_dof_handler);
    std::vector<unsigned int> stokes_sub_blocks (dim+1,0);
    stokes_sub_blocks[dim] = 1;
    DoFRenumbering::component_wise (stokes_dof_handler, stokes_sub_blocks);

    temperature_dof_handler.distribute_dofs (temperature_fe);

    DoFRenumbering::hierarchical (temperature_dof_handler);
    std::vector<unsigned int> stokes_dofs_per_block (2);
    DoFTools::count_dofs_per_block (stokes_dof_handler, stokes_dofs_per_block,
                                    stokes_sub_blocks);

    const unsigned int n_u = stokes_dofs_per_block[0],
                       n_p = stokes_dofs_per_block[1],
                       n_T = temperature_dof_handler.n_dofs();

    // print dof numbers with 1000s
    // separator since they are frequently
    // large
    std::locale s = pcout.get_stream().getloc();
    // Creating std::locale with an empty string causes problems
    // on some platforms, so catch the exception and ignore
    try
      {
        pcout.get_stream().imbue(std::locale(""));
      }
    catch (std::runtime_error e)
      {
        // If the locale doesn't work, just give up
      }
    pcout << "Number of active cells: "
          << triangulation.n_global_active_cells()
          << " (on "
          << triangulation.n_levels()
          << " levels)"
          << std::endl
          << "Number of degrees of freedom: "
          << n_u + n_p + n_T
          << " (" << n_u << '+' << n_p << '+'<< n_T <<')'
          << std::endl
          << std::endl;
    pcout.get_stream().imbue(s);



    std::vector<IndexSet> stokes_partitioning, stokes_relevant_partitioning;
    IndexSet temperature_partitioning (n_T), temperature_relevant_partitioning (n_T);
    IndexSet stokes_relevant_set;
    {
      IndexSet stokes_index_set = stokes_dof_handler.locally_owned_dofs();
      stokes_partitioning.push_back(stokes_index_set.get_view(0,n_u));
      stokes_partitioning.push_back(stokes_index_set.get_view(n_u,n_u+n_p));

      DoFTools::extract_locally_relevant_dofs (stokes_dof_handler,
                                               stokes_relevant_set);
      stokes_relevant_partitioning.push_back(stokes_relevant_set.get_view(0,n_u));
      stokes_relevant_partitioning.push_back(stokes_relevant_set.get_view(n_u,n_u+n_p));

      temperature_partitioning = temperature_dof_handler.locally_owned_dofs();
      DoFTools::extract_locally_relevant_dofs (temperature_dof_handler,
                                               temperature_relevant_partitioning);
    }

    {

      stokes_constraints.clear ();
      stokes_constraints.reinit (stokes_relevant_set);

      DoFTools::make_hanging_node_constraints (stokes_dof_handler,
                                               stokes_constraints);

      // obtain the boundary indicators that belong to zero velocity
      // and no-normal-flux type
      typedef unsigned char boundary_indicator_t;
      const std::set<boundary_indicator_t>
      zero_boundary_indicators
        = geometry_model->get_zero_velocity_boundary_indicators ();

      const std::set<boundary_indicator_t>
      no_normal_flux_boundary_indicators
        = geometry_model->get_tangential_velocity_boundary_indicators ();

      // also make sure that the same indicator isn't used for both
      {
        std::vector<boundary_indicator_t> intersection;
        std::set_intersection (zero_boundary_indicators.begin(),
                               zero_boundary_indicators.end(),
                               no_normal_flux_boundary_indicators.begin(),
                               no_normal_flux_boundary_indicators.end(),
                               std::back_inserter(intersection));
        Assert (intersection.size() == 0,
                ExcMessage ("Some boundary indicators are declared to be both "
                            "zero velocity and tangential velocity. That can't work!"));
      }

      // do the interpolation for zero velocity
      std::vector<bool> velocity_mask (dim+1, true);
      velocity_mask[dim] = false;
      for (std::set<boundary_indicator_t>::const_iterator
           p = zero_boundary_indicators.begin();
           p != zero_boundary_indicators.end(); ++p)
        VectorTools::interpolate_boundary_values (stokes_dof_handler,
                                                  *p,
                                                  ZeroFunction<dim>(dim+1),
                                                  stokes_constraints,
                                                  velocity_mask);

      // do the same for no-normal-flux boundaries
      VectorTools::compute_no_normal_flux_constraints (stokes_dof_handler,
                                                       /* first_vector_component= */ 0,
                                                       no_normal_flux_boundary_indicators,
                                                       stokes_constraints,
                                                       mapping);
      stokes_constraints.close ();
    }
    {
      temperature_constraints.clear ();
      temperature_constraints.reinit (temperature_relevant_partitioning);

      DoFTools::make_hanging_node_constraints (temperature_dof_handler,
                                               temperature_constraints);

      // obtain the boundary indicators that belong to Dirichlet-type
      // temperature boundary conditions and interpolate the temperature
      // there
      const std::set<unsigned char>
      temperature_dirichlet_boundary_indicators
        = geometry_model->get_temperature_dirichlet_boundary_indicators ();
      for (std::set<unsigned char>::const_iterator
           p = temperature_dirichlet_boundary_indicators.begin();
           p != temperature_dirichlet_boundary_indicators.end(); ++p)
        VectorTools::interpolate_boundary_values (temperature_dof_handler,
                                                  *p,
                                                  ScalarFunctionFromFunctionObject<dim>(std_cxx1x::bind (&BoundaryTemperature::Interface<dim>::temperature,
                                                                                        std_cxx1x::cref(*boundary_temperature),
                                                                                        std_cxx1x::cref(*geometry_model),
                                                                                        *p,
                                                                                        std_cxx1x::_1)),
                                                  temperature_constraints);
      temperature_constraints.close ();
    }

    setup_stokes_matrix (stokes_partitioning);
    setup_stokes_preconditioner (stokes_partitioning);
    setup_temperature_matrix (temperature_partitioning);

    stokes_rhs.reinit (stokes_partitioning, MPI_COMM_WORLD);
    stokes_rhs_helper.reinit (stokes_partitioning, MPI_COMM_WORLD);
    stokes_solution.reinit (stokes_relevant_partitioning, MPI_COMM_WORLD);
    old_stokes_solution.reinit (stokes_solution);

    temperature_rhs.reinit (temperature_partitioning, MPI_COMM_WORLD);
    temperature_solution.reinit (temperature_relevant_partitioning, MPI_COMM_WORLD);
    old_temperature_solution.reinit (temperature_solution);
    old_old_temperature_solution.reinit (temperature_solution);

    rebuild_stokes_matrix              = true;
    rebuild_stokes_preconditioner      = true;

    computing_timer.exit_section();
  }



  template <int dim>
  void Simulator<dim>::postprocess ()
  {
    computing_timer.enter_section ("Postprocessing");
    pcout << "   Postprocessing:" << std::endl;

    // run all the postprocessing routines and then write
    // the current state of the statistics table to a file
    std::list<std::pair<std::string,std::string> >
    output_list = postprocess_manager.execute (statistics);

    std::ofstream stat_file ((parameters.output_directory+"statistics").c_str());
    statistics.set_scientific("Time (years)", true);
    statistics.set_scientific("Time step size (year)", true);
    statistics.write_text (stat_file,
                           TableHandler::table_with_separate_column_description);

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
    computing_timer.exit_section ();
  }



// Contrary to step-32, we have found that just refining by the temperature
// works well in 2d, but only leads to refinement in the boundary
// layer at the core-mantle boundary in 3d. consequently, we estimate
// the error based both on the temperature and on the velocity; the
// vectors with the resulting error indicators are then both normalized
// to a maximal value of one, and we take the maximum of the two indicators
// to decide whether we want to refine or not. this ensures that we
// also refine into plumes where maybe the temperature gradients aren't
// as strong as in the boundary layer but where nevertheless the gradients
// in the velocity are large
  template <int dim>
  void Simulator<dim>::refine_mesh (const unsigned int max_grid_level)
  {
    computing_timer.enter_section ("Refine mesh structure, part 1");

    Vector<float> estimated_error_per_cell (triangulation.n_active_cells());

    // compute the errors for
    // temperature and stokes solution,
    // then scale them and find the
    // maximum between the two
    {
      Vector<float> estimated_error_per_cell_T (triangulation.n_active_cells());

      KellyErrorEstimator<dim>::estimate (temperature_dof_handler,
                                          QGauss<dim-1>(parameters.temperature_degree+1),
                                          typename FunctionMap<dim>::type(),
                                          temperature_solution,
                                          estimated_error_per_cell_T,
                                          std::vector<bool>(),
                                          0,
                                          0,
                                          triangulation.locally_owned_subdomain());
      estimated_error_per_cell_T /= Utilities::MPI::max (estimated_error_per_cell_T.linfty_norm(),
                                                         MPI_COMM_WORLD);

      Vector<float> estimated_error_per_cell_u (triangulation.n_active_cells());
      std::vector<bool> velocity_mask (dim+1, true);
      velocity_mask[dim] = false;
      KellyErrorEstimator<dim>::estimate (stokes_dof_handler,
                                          QGauss<dim-1>(parameters.stokes_velocity_degree+1),
                                          typename FunctionMap<dim>::type(),
                                          stokes_solution,
                                          estimated_error_per_cell_u,
                                          velocity_mask,
                                          0,
                                          0,
                                          triangulation.locally_owned_subdomain());
      estimated_error_per_cell_u /= Utilities::MPI::max (estimated_error_per_cell_u.linfty_norm(),
                                                         MPI_COMM_WORLD);

      for (unsigned int i=0; i<estimated_error_per_cell.size(); ++i)
        estimated_error_per_cell(i) = std::max (estimated_error_per_cell_T(i),
                                                estimated_error_per_cell_u(i));
    }

    parallel::distributed::GridRefinement::
    refine_and_coarsen_fixed_fraction (triangulation,
                                       estimated_error_per_cell,
                                       parameters.refinement_fraction,
                                       parameters.coarsening_fraction);

    // limit maximum refinement level
    if (triangulation.n_levels() > max_grid_level)
      for (typename Triangulation<dim>::active_cell_iterator
           cell = triangulation.begin_active(max_grid_level);
           cell != triangulation.end(); ++cell)
        cell->clear_refine_flag ();

    std::vector<const TrilinosWrappers::MPI::Vector *> x_temperature (2);
    x_temperature[0] = &temperature_solution;
    x_temperature[1] = &old_temperature_solution;
    std::vector<const TrilinosWrappers::MPI::BlockVector *> x_stokes (2);
    x_stokes[0] = &stokes_solution;
    x_stokes[1] = &old_stokes_solution;

    parallel::distributed::SolutionTransfer<dim,TrilinosWrappers::MPI::Vector>
    temperature_trans(temperature_dof_handler);
    parallel::distributed::SolutionTransfer<dim,TrilinosWrappers::MPI::BlockVector>
    stokes_trans(stokes_dof_handler);

    triangulation.prepare_coarsening_and_refinement();
    temperature_trans.prepare_for_coarsening_and_refinement(x_temperature);
    stokes_trans.prepare_for_coarsening_and_refinement(x_stokes);

    triangulation.execute_coarsening_and_refinement ();
    global_volume = GridTools::volume (triangulation, mapping);
    computing_timer.exit_section();

    setup_dofs ();

    computing_timer.enter_section ("Refine mesh structure, part 2");

    {
      TrilinosWrappers::MPI::Vector
      distributed_temp1 (temperature_rhs);
      TrilinosWrappers::MPI::Vector
      distributed_temp2 (temperature_rhs);

      std::vector<TrilinosWrappers::MPI::Vector *> tmp (2);
      tmp[0] = &(distributed_temp1);
      tmp[1] = &(distributed_temp2);
      temperature_trans.interpolate(tmp);

      temperature_solution     = distributed_temp1;
      old_temperature_solution = distributed_temp2;
    }

    {
      TrilinosWrappers::MPI::BlockVector
      distributed_stokes (stokes_rhs);
      TrilinosWrappers::MPI::BlockVector
      old_distributed_stokes (stokes_rhs);
      std::vector<TrilinosWrappers::MPI::BlockVector *> stokes_tmp (2);
      stokes_tmp[0] = &(distributed_stokes);
      stokes_tmp[1] = &(old_distributed_stokes);

      stokes_trans.interpolate (stokes_tmp);
      stokes_solution     = distributed_stokes;
      old_stokes_solution = old_distributed_stokes;
    }

    computing_timer.exit_section();
  }



  /**
   * This is the main function of the program, containing the overall
   * logic which function is called when.
   */
  template <int dim>
  void Simulator<dim>::run ()
  {
    // if we want to resume a computation from an earlier point
    // then reload it from a snapshot. otherwise do the basic
    // start-up
    if (parameters.resume_computation == true)
      {
        resume_from_snapshot();
      }
    else
      {
        triangulation.refine_global (parameters.initial_global_refinement);
        global_volume = GridTools::volume (triangulation, mapping);

        setup_dofs();
      }

    unsigned int max_refinement_level = parameters.initial_global_refinement +
                                        parameters.initial_adaptive_refinement;
    unsigned int pre_refinement_step = 0;

  start_time_iteration:

    if (parameters.resume_computation == false)
      {
        set_initial_temperature_field ();
        compute_initial_pressure_field ();

        time                      = 0;
        timestep_number           = 0;
        time_step = old_time_step = 0;
      }

    // start the principal loop over time steps:
    do
      {
        pcout << "*** Timestep " << timestep_number
              << ":  t=" << time/year_in_seconds
              << " years"
              << std::endl;

        // set global statistics about this time step
        statistics.add_value("Time step number", timestep_number);
        statistics.add_value("Time (years)", time / year_in_seconds);


        // then do the core work: assemble systems and solve
        assemble_stokes_system ();
        build_stokes_preconditioner ();
        solve ();

        pcout << std::endl;

        // see if we have to start over with a new refinement cycle
        // at the beginning of the simulation
        if ((timestep_number == 0) &&
            (pre_refinement_step < parameters.initial_adaptive_refinement))
          {
            refine_mesh (max_refinement_level);
            ++pre_refinement_step;
            goto start_time_iteration;
          }

        postprocess ();

        // see if this is a time step where additional refinement is requested
        // if so, then loop over as many times as this is necessary
        if ((parameters.additional_refinement_times.size() > 0)
            &&
            (parameters.additional_refinement_times.front () < time+time_step))
          {
            while ((parameters.additional_refinement_times.size() > 0)
                   &&
                   (parameters.additional_refinement_times.front () < time+time_step))
              {
                ++max_refinement_level;
                refine_mesh (max_refinement_level);

                parameters.additional_refinement_times
                .erase (parameters.additional_refinement_times.begin());
              }
          }
        else
          // see if this is a time step where regular refinement is necessary, but only
          // if the previous rule wasn't triggered
          if ((timestep_number > 0)
              &&
              (timestep_number % parameters.adaptive_refinement_interval == 0))
            refine_mesh (max_refinement_level);


        // every 100 time steps output a summary of the current
        // timing information
        if ((timestep_number > 0) && (timestep_number % 100 == 0))
          computing_timer.print_summary ();

        // increment time step by one. then prepare
        // for the next time step by shifting solution vectors
        // by one time step and presetting the current solution based
        // on an extrapolation
        time += time_step;
        ++timestep_number;
        {
          TrilinosWrappers::MPI::BlockVector old_old_stokes_solution;
          old_old_stokes_solution      = old_stokes_solution;
          old_stokes_solution          = stokes_solution;
          old_old_temperature_solution = old_temperature_solution;
          old_temperature_solution     = temperature_solution;
          if (old_time_step > 0)
            {
              stokes_solution.sadd (1.+time_step/old_time_step, -time_step/old_time_step,
                                    old_old_stokes_solution);
              temperature_solution.sadd (1.+time_step/old_time_step,
                                         -time_step/old_time_step,
                                         old_old_temperature_solution);
            }
        }

        // periodically generate snapshots so that we can resume here
        // if the program aborts or is terminated
        if (timestep_number % 50 == 0)
          {
            create_snapshot();
            // matrices will be regenerated after a resume, so do that here too
            // to be consistent. otherwise we would get different results
            // for a restarted computation than for one that ran straight
            // through
            rebuild_stokes_matrix =
              rebuild_stokes_preconditioner = true;
          }
      }
    while (time < parameters.end_time * year_in_seconds);
  }
}



// explicit instantiation of the functions we implement in this file
namespace aspect
{
  template
  class Simulator<deal_II_dimension>;
}
