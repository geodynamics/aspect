/*
  Copyright (C) 2025 - by the authors of the ASPECT code.

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

#include <aspect/simulator/solver/stokes_direct.h>

#include <deal.II/lac/trilinos_solver.h>

namespace aspect
{
  namespace StokesSolver
  {
    template <int dim>
    SolverOutputs
    Direct<dim>::solve(const LinearAlgebra::BlockSparseMatrix &system_matrix,
                       const LinearAlgebra::BlockVector &system_rhs,
                       const bool solve_newton_system,
                       const double last_pressure_normalization_adjustment,
                       LinearAlgebra::BlockVector &solution_vector)
    {

      // In the following, we will operate on a vector that contains only
      // the velocity and pressure DoFs, rather than on the full
      // system. Set such a reduced vector up, without any ghost elements.
      // (Worth noting: for direct solvers, this vector has one block,
      // whereas for the iterative solvers, the result has two blocks.)
      LinearAlgebra::BlockVector distributed_stokes_solution (this->introspection().index_sets.stokes_partitioning,
                                                              this->get_mpi_communicator());

      // We will need the Stokes block indices a lot below, shorten their names
      const unsigned int velocity_block_index = this->introspection().block_indices.velocities;
      const unsigned int pressure_block_index = (this->get_parameters().include_melt_transport) ?
                                                this->introspection().variable("fluid pressure").block_index
                                                : this->introspection().block_indices.pressure;
      (void) velocity_block_index;
      (void) pressure_block_index;


      // Create a view of all constraints that only pertains to the
      // Stokes subset of degrees of freedom. We can then use this later
      // to call constraints.distribute(), constraints.set_zero(), etc.,
      // on those block vectors that only have the Stokes components in
      // them.
      //
      // For the moment, assume that the Stokes degrees are first in the
      // overall vector, so that they form a contiguous range starting
      // at zero. The assertion checks this, but this could easily be
      // generalized if the Stokes block were not starting at zero.
#if DEAL_II_VERSION_GTE(9,6,0)
      {
        Assert (velocity_block_index == 0, ExcNotImplemented());
      }

      IndexSet stokes_dofs (this->get_dof_handler().n_dofs());
      stokes_dofs.add_range (0, distributed_stokes_solution.size());
      const AffineConstraints<double> current_stokes_constraints
        = this->get_current_constraints().get_view (stokes_dofs);
#else
      const AffineConstraints<double> &current_stokes_constraints = this->get_current_constraints();
#endif

      double initial_nonlinear_residual = numbers::signaling_nan<double>();
      double final_linear_residual      = numbers::signaling_nan<double>();


      // We hard-code the blocks down below, so make sure block 0 is indeed
      // the block containing velocity and pressure:
      Assert(distributed_stokes_solution.n_blocks() == 1, ExcInternalError());
      Assert(velocity_block_index == 0, ExcNotImplemented());
      Assert(pressure_block_index == 0
             ||
             (this->get_parameters().include_melt_transport
              && this->introspection().variable("fluid pressure").block_index == 0
              && this->introspection().variable("total pressure").block_index == 0),
             ExcNotImplemented());

      // Clarify that we only use one block for the direct solver
      const unsigned int velocity_and_pressure_block = velocity_block_index;

      // Start with a reasonable guess.
      //
      // While we don't need to set up the initial guess for the direct solver
      // (it will be ignored by the solver anyway), we need this if we are
      // using a nonlinear scheme, because we use this to compute the current
      // nonlinear residual (see initial_residual below).
      solution_vector.block(velocity_and_pressure_block) = this->get_current_linearization_point().block(velocity_and_pressure_block);

      // TODO: if there was an easy way to know if the caller needs the
      // initial residual we could skip all of this stuff.
      distributed_stokes_solution.block(velocity_and_pressure_block) = solution_vector.block(velocity_and_pressure_block);
      this->denormalize_pressure (last_pressure_normalization_adjustment,
                                  distributed_stokes_solution);
      current_stokes_constraints.set_zero (distributed_stokes_solution);

      // Undo the pressure scaling:
      const IndexSet &pressure_idxset = this->get_parameters().include_melt_transport ?
                                        this->introspection().index_sets.locally_owned_melt_pressure_dofs
                                        : this->introspection().index_sets.locally_owned_pressure_dofs;

      for (unsigned int i=0; i< pressure_idxset.n_elements(); ++i)
        {
          types::global_dof_index idx = pressure_idxset.nth_index_in_set(i);

          distributed_stokes_solution(idx) /= this->get_pressure_scaling();
        }
      distributed_stokes_solution.compress(VectorOperation::insert);

      // we need a temporary vector for the residual (even if we don't care about it)
      LinearAlgebra::Vector residual (this->introspection().index_sets.stokes_partitioning[0], this->get_mpi_communicator());

      initial_nonlinear_residual = system_matrix.block(velocity_and_pressure_block,velocity_and_pressure_block).residual(
                                     residual,
                                     distributed_stokes_solution.block(velocity_and_pressure_block),
                                     system_rhs.block(velocity_and_pressure_block));

      SolverControl cn;
      // TODO: can we re-use the direct solver?
      TrilinosWrappers::SolverDirect solver(cn);
      try
        {
          solver.solve(system_matrix.block(velocity_and_pressure_block,velocity_and_pressure_block),
                       distributed_stokes_solution.block(velocity_and_pressure_block),
                       system_rhs.block(velocity_and_pressure_block));

          // if we got here, we have successfully solved the linear system
          // with a direct solver, and the final linear residual should
          // be approximately zero. we could compute it exactly, but
          // this is probably not necessary
          final_linear_residual = 0;
        }
      // if the solver fails, report the error from processor 0 with some additional
      // information about its location, and throw a quiet exception on all other
      // processors
      catch (const std::exception &exc)
        {
          if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
            {
              AssertThrow (false,
                           ExcMessage (std::string("The direct Stokes solver "
                                                   "did not succeed. It reported the following error:\n\n")
                                       +
                                       exc.what()));
            }
          else
            throw QuietException();
        }


      current_stokes_constraints.distribute (distributed_stokes_solution);

      // Now rescale the pressure back to real physical units. Note that we are
      // working on a vector in which all velocities and pressures are in one
      // block (that's the design for block layout in case we're using a direct
      // solver), and so unlike in the "common" case, we can't just scale a
      // whole vector block -- we have to do it element by element.
      {
        const IndexSet &pressure_idxset
          = (this->get_parameters().include_melt_transport ?
             this->introspection().index_sets.locally_owned_melt_pressure_dofs
             : this->introspection().index_sets.locally_owned_pressure_dofs);
        for (const types::global_dof_index i : pressure_idxset)
          distributed_stokes_solution(i) *= this->get_pressure_scaling();

        distributed_stokes_solution.block(velocity_and_pressure_block).compress(VectorOperation::insert);
      }

      // Then copy back the solution from the temporary (non-ghosted) vector
      // into the ghosted one with all solution components. Note that
      // for a direct solver, we have only one block for velocity+pressure,
      // and so only one block needs to be copied.
      solution_vector.block(velocity_and_pressure_block) = distributed_stokes_solution.block(velocity_and_pressure_block);

      // do some cleanup now that we have the solution
      this->remove_nullspace(solution_vector, distributed_stokes_solution);

      SolverOutputs outputs;

      if (solve_newton_system == false)
        outputs.pressure_normalization_adjustment = this->normalize_pressure(solution_vector);

      this->get_pcout() << "done." << std::endl;

      outputs.initial_nonlinear_residual = initial_nonlinear_residual;
      outputs.final_linear_residual = final_linear_residual;

      return outputs;
    }

    template <int dim>
    std::string
    Direct<dim>::name () const
    {
      return "direct";
    }


    // explicit instantiation of the functions we implement in this file
#define INSTANTIATE(dim) \
  template class Direct<dim>;

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
