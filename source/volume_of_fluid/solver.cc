/*
 Copyright (C) 2016 - 2020 by the authors of the ASPECT code.

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

#include <aspect/global.h>
#include <aspect/volume_of_fluid/handler.h>

#include <deal.II/lac/affine_constraints.h>

#ifdef ASPECT_USE_PETSC
#include <deal.II/lac/solver_cg.h>
#else
#include <deal.II/lac/trilinos_solver.h>
#endif

#include <deal.II/fe/fe_values.h>

namespace aspect
{
  template <int dim>
  void VolumeOfFluidHandler<dim>::solve_volume_of_fluid_system (const VolumeOfFluidField<dim> &field)
  {
    const unsigned int block_idx = field.volume_fraction.block_index;

    TimerOutput::Scope timer (sim.computing_timer, "Solve volume of fluid system");
    this->get_pcout() << "   Solving volume of fluid system... " << std::flush;

    const double tolerance = std::max(1e-50,
                                      volume_of_fluid_solver_tolerance*sim.system_rhs.block(block_idx).l2_norm());

    SolverControl solver_control (1000, tolerance);

#ifdef ASPECT_USE_PETSC
    SolverCG<LinearAlgebra::Vector> solver(solver_control);
    LinearAlgebra::PreconditionJacobi precondition;
    precondition.initialize(sim.system_matrix.block(block_idx, block_idx));
#else
    TrilinosWrappers::SolverCG solver(solver_control);
    TrilinosWrappers::PreconditionJacobi precondition;
    precondition.initialize(sim.system_matrix.block(block_idx, block_idx));
#endif

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
          AssertThrow (false,
                       ExcMessage (std::string("The iterative advection solver "
                                               "did not converge. It reported the following error:\n\n")
                                   +
                                   exc.what()))
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
  }
}

namespace aspect
{
#define INSTANTIATE(dim) \
  template void VolumeOfFluidHandler<dim>::solve_volume_of_fluid_system (const VolumeOfFluidField<dim> &field);


  ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
}
