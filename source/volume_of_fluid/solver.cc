// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2016-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

#include <aspect/global.h>
#include <aspect/volume_of_fluid/handler.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/solver_cg.h>

#include <deal.II/fe/fe_values.h>

namespace aspect
{
  template <int dim>
  void VolumeOfFluidHandler<dim>::solve_volume_of_fluid_system (const VolumeOfFluidField<dim> &field)
  {
    const unsigned int block_idx = field.volume_fraction.block_index;

    this->get_computing_timer().enter_subsection("Solve volume of fluid system");
    this->get_pcout() << "   Solving volume of fluid system... " << std::flush;

    const double tolerance = std::max(1e-50,
                                      volume_of_fluid_solver_tolerance*sim.system_rhs.block(block_idx).l2_norm());

    SolverControl solver_control (1000, tolerance);
    SolverCG<LinearAlgebra::Vector> solver(solver_control);

    LinearAlgebra::PreconditionJacobi precondition;
    precondition.initialize(sim.system_matrix.block(block_idx, block_idx));

    // Create distributed vector (we need all blocks here even though we only
    // solve for the current block) because only have a AffineConstraints<double>
    // for the whole system, current_linearization_point contains our initial guess.
    LinearAlgebra::BlockVector distributed_solution (
      this->introspection().index_sets.system_partitioning,
      this->get_mpi_communicator());
    distributed_solution.block(block_idx) = sim.current_linearization_point.block (block_idx);

    sim.current_constraints.set_zero(distributed_solution);

    // solve the linear system:
    try
      {
        solver.solve (sim.system_matrix.block(block_idx,block_idx),
                      distributed_solution.block(block_idx),
                      sim.system_rhs.block(block_idx),
                      precondition);
      }
    // if the solver fails, report the error from processor 0 with some additional
    // information about its location, and throw a quiet exception on all other
    // processors
    catch (const std::exception &exc)
      {
        if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
          {
            AssertThrow (false,
                         ExcMessage (std::string("The iterative advection solver "
                                                 "did not converge. It reported the following error:\n\n")
                                     +
                                     exc.what()));
          }
        else
          throw QuietException();
      }

    sim.current_constraints.distribute (distributed_solution);
    sim.solution.block(block_idx) = distributed_solution.block(block_idx);

    // print number of iterations and also record it in the
    // statistics file
    this->get_pcout() << solver_control.last_step()
                      << " iterations." << std::endl;

    // Do not add VolumeOfFluid solver iterations to statistics, duplication due to
    // dimensional splitting results in incorrect line formatting (lines of
    // data split inconsistently with missing values)

    this->get_computing_timer().leave_subsection("Solve volume of fluid system");
  }
}

namespace aspect
{
#define INSTANTIATE(dim) \
  template void VolumeOfFluidHandler<dim>::solve_volume_of_fluid_system (const VolumeOfFluidField<dim> &field);


  ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
}
