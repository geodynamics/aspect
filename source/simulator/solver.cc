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
#include <aspect/melt.h>
#include <aspect/simulator/solver/block_stokes_preconditioner.h>
#include <aspect/simulator/solver/stokes_matrix_based.h>
#include <aspect/simulator/solver/stokes_matrix_free.h>
#include <aspect/simulator/solver/stokes_direct.h>
#include <aspect/mesh_deformation/interface.h>

#include <deal.II/base/signaling_nan.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/solver_bicgstab.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/lac/trilinos_vector.h>

namespace aspect
{
  template <int dim>
  double Simulator<dim>::solve_advection (const AdvectionField &advection_field)
  {
    const unsigned int block_idx = advection_field.block_index(introspection);

    const std::string field_name = (advection_field.is_temperature()
                                    ?
                                    "temperature"
                                    :
                                    introspection.name_for_compositional_index(advection_field.compositional_variable) + " composition");

    const double advection_solver_tolerance = (advection_field.is_temperature()) ? (parameters.temperature_solver_tolerance) : (parameters.composition_solver_tolerance);

    const double tolerance = std::max(1e-50,
                                      advection_solver_tolerance*system_rhs.block(block_idx).l2_norm());

    SolverControl solver_control (1000, tolerance);

    solver_control.enable_history_data();

    SolverGMRES<LinearAlgebra::Vector> solver (solver_control,
                                               SolverGMRES<LinearAlgebra::Vector>::AdditionalData(parameters.advection_gmres_restart_length,true));

    // check if matrix and/or RHS are zero
    // note: to avoid a warning, we compare against numeric_limits<double>::min() instead of 0 here
    if (system_rhs.block(block_idx).l2_norm() <= std::numeric_limits<double>::min())
      {
        pcout << "   Skipping " + field_name + " solve because RHS is zero." << std::endl;
        solution.block(block_idx) = 0;

        // signal successful solver and signal residual of zero
        solver_control.check(0, 0.0);
        signals.post_advection_solver(*this,
                                      advection_field.is_temperature(),
                                      advection_field.compositional_variable,
                                      solver_control);

        return 0;
      }

    AssertThrow(system_matrix.block(block_idx,
                                    block_idx).linfty_norm() > std::numeric_limits<double>::min(),
                ExcMessage ("The " + field_name + " equation can not be solved, because the matrix is zero, "
                            "but the right-hand side is nonzero."));

    LinearAlgebra::PreconditionILU preconditioner;
    // first build without diagonal strengthening:
    build_advection_preconditioner(advection_field, preconditioner, 0.);

    computing_timer.enter_subsection(advection_field.is_temperature() ?
                                     "Solve temperature system" :
                                     "Solve composition system");

    if (advection_field.is_temperature())
      {
        pcout << "   Solving temperature system... " << std::flush;
      }
    else
      {
        pcout << "   Solving "
              << introspection.name_for_compositional_index(advection_field.compositional_variable)
              << " system "
              << "... " << std::flush;
      }

    // Create distributed vector (we need all blocks here even though we only
    // solve for the current block) because we only have an AffineConstraints object
    // for the whole system, current_linearization_point contains our initial guess.
    LinearAlgebra::BlockVector distributed_solution (
      introspection.index_sets.system_partitioning,
      mpi_communicator);
    distributed_solution.block(block_idx) = current_linearization_point.block (block_idx);

    // Temporary vector to hold the residual, we don't need a BlockVector here.
    LinearAlgebra::Vector temp (
      introspection.index_sets.system_partitioning[block_idx],
      mpi_communicator);

    current_constraints.set_zero(distributed_solution);

    // Compute the residual before we solve and return this at the end.
    // This is used in the nonlinear solver.
    const double initial_residual = system_matrix.block(block_idx,block_idx).residual
                                    (temp,
                                     distributed_solution.block(block_idx),
                                     system_rhs.block(block_idx));

    // solve the linear system:
    try
      {
        try
          {
            solver.solve (system_matrix.block(block_idx,block_idx),
                          distributed_solution.block(block_idx),
                          system_rhs.block(block_idx),
                          preconditioner);
          }
        catch (const std::exception &exc)
          {
            // Try rebuilding the preconditioner with diagonal strengthening. In general,
            // this increases the number of iterations needed, but helps in rare situations,
            // especially when SUPG is used.
            pcout << "retrying linear solve with different preconditioner..." << std::endl;
            build_advection_preconditioner(advection_field, preconditioner, 1e-5);
            solver.solve (system_matrix.block(block_idx,block_idx),
                          distributed_solution.block(block_idx),
                          system_rhs.block(block_idx),
                          preconditioner);
          }

      }
    // if the solver fails, report the error from processor 0 with some additional
    // information about its location, and throw a quiet exception on all other
    // processors
    catch (const std::exception &exc)
      {
        // signal unsuccessful solver
        signals.post_advection_solver(*this,
                                      advection_field.is_temperature(),
                                      advection_field.compositional_variable,
                                      solver_control);


        Utilities::throw_linear_solver_failure_exception("iterative advection solver",
                                                         "Simulator::solve_advection",
                                                         std::vector<SolverControl> {solver_control},
                                                         exc,
                                                         mpi_communicator,
                                                         parameters.output_directory+"solver_history.txt");
      }

    // signal successful solver
    signals.post_advection_solver(*this,
                                  advection_field.is_temperature(),
                                  advection_field.compositional_variable,
                                  solver_control);

    current_constraints.distribute (distributed_solution);
    solution.block(block_idx) = distributed_solution.block(block_idx);

    // print number of iterations and also record it in the
    // statistics file
    pcout << solver_control.last_step()
          << " iterations." << std::endl;

    if ((advection_field.is_discontinuous(introspection)
         &&
         (
           (advection_field.is_temperature() && parameters.use_limiter_for_discontinuous_temperature_solution)
           ||
           (!advection_field.is_temperature() && parameters.use_limiter_for_discontinuous_composition_solution[advection_field.compositional_variable])
         )))
      {
        apply_limiter_to_dg_solutions(advection_field);

        computing_timer.leave_subsection(advection_field.is_temperature() ?
                                         "Solve temperature system" :
                                         "Solve composition system");

        // by applying the limiter we have modified the solution to no longer
        // satisfy the equation. Therefore the residual is meaningless and cannot
        // converge to zero in nonlinear iterations. Disable residual computation
        // for this field.
        return 0.0;
      }

    computing_timer.leave_subsection(advection_field.is_temperature() ?
                                     "Solve temperature system" :
                                     "Solve composition system");

    return initial_residual;
  }



  template <int dim>
  std::pair<double,double>
  Simulator<dim>::solve_stokes (LinearAlgebra::BlockVector &solution_vector)
  {
    computing_timer.enter_subsection("Solve Stokes system");

    const std::string name = [&]() -> std::string
    {
      if (parameters.stokes_solver_type == Parameters<dim>::StokesSolverType::block_gmg)
        return stokes_matrix_free->name();
      if (parameters.use_direct_stokes_solver)
        return stokes_direct->name();

      return stokes_matrix_based->name();
    }();

    pcout << "   Solving Stokes system (" << name << ")... " << std::flush;

    StokesSolver::SolverOutputs outputs;

    if (parameters.stokes_solver_type == Parameters<dim>::StokesSolverType::block_gmg)
      {
        outputs = stokes_matrix_free->solve(system_matrix,
                                            system_rhs,
                                            assemble_newton_stokes_system,
                                            last_pressure_normalization_adjustment,
                                            solution_vector);
      }
    else if (parameters.use_direct_stokes_solver)
      {
        outputs = stokes_direct->solve(system_matrix,
                                       system_rhs,
                                       assemble_newton_stokes_system,
                                       last_pressure_normalization_adjustment,
                                       solution_vector);
      }
    else
      {
        outputs = stokes_matrix_based->solve(system_matrix,
                                             system_rhs,
                                             assemble_newton_stokes_system,
                                             last_pressure_normalization_adjustment,
                                             solution_vector);
      }

    last_pressure_normalization_adjustment = outputs.pressure_normalization_adjustment;

    // convert melt pressures:
    if (parameters.include_melt_transport)
      melt_handler->compute_melt_variables(system_matrix,solution_vector,system_rhs);

    computing_timer.leave_subsection("Solve Stokes system");

    return {outputs.initial_nonlinear_residual,
            outputs.final_linear_residual
           };
  }

}



// explicit instantiation of the functions we implement in this file
namespace aspect
{
#define INSTANTIATE(dim) \
  template double Simulator<dim>::solve_advection (const AdvectionField &); \
  template std::pair<double,double> Simulator<dim>::solve_stokes (LinearAlgebra::BlockVector &solution_vector);

  ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
}
