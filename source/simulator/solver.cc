/*
  Copyright (C) 2011, 2012, 2013 by the authors of the ASPECT code.

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
/*  $Id$  */


#include <aspect/simulator.h>
#include <aspect/global.h>

#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/constraint_matrix.h>

#ifdef USE_PETSC
#include <deal.II/lac/solver_cg.h>
#else
#include <deal.II/lac/trilinos_solver.h>
#endif

#include <deal.II/lac/pointer_matrix.h>


namespace aspect
{
  namespace internal
  {
    using namespace dealii;

    /**
     * Implement multiplication with Stokes part of system matrix. In essence, this
     * object represents a 2x2 block matrix that corresponds to the top left
     * sub-blocks of the entire system matrix (i.e., the Stokes part)
     */
    class StokesBlock
    {
      public:
        /**
         * @brief Constructor
         *
         * @param S The entire system matrix
         */
        StokesBlock (const LinearAlgebra::BlockSparseMatrix  &S)
          : system_matrix(S) {}

        /**
         * Matrix vector product with Stokes block.
         */
        void vmult (LinearAlgebra::BlockVector       &dst,
                    const LinearAlgebra::BlockVector &src) const;

        void Tvmult (LinearAlgebra::BlockVector       &dst,
                     const LinearAlgebra::BlockVector &src) const;

        void vmult_add (LinearAlgebra::BlockVector       &dst,
                        const LinearAlgebra::BlockVector &src) const;

        void Tvmult_add (LinearAlgebra::BlockVector       &dst,
                         const LinearAlgebra::BlockVector &src) const;

        /**
         * Compute the residual with the Stokes block. In a departure from
         * the other functions, the #b variable may actually have more than
         * two blocks so that we can put it a global system_rhs vector. The
         * other vectors need to have 2 blocks only.
         */
        double residual (LinearAlgebra::BlockVector       &dst,
                         const LinearAlgebra::BlockVector &x,
                         const LinearAlgebra::BlockVector &b) const;


      private:

        /**
         * Reference to the system matrix object.
         */
        const LinearAlgebra::BlockSparseMatrix &system_matrix;
    };



    void StokesBlock::vmult (LinearAlgebra::BlockVector       &dst,
                             const LinearAlgebra::BlockVector &src) const
    {
      Assert (src.n_blocks() == 2, ExcInternalError());
      Assert (dst.n_blocks() == 2, ExcInternalError());

      system_matrix.block(0,0).vmult(dst.block(0), src.block(0));
      system_matrix.block(0,1).vmult_add(dst.block(0), src.block(1));

      system_matrix.block(1,0).vmult(dst.block(1), src.block(0));
      system_matrix.block(1,1).vmult_add(dst.block(1), src.block(1));
    }


    void StokesBlock::Tvmult (LinearAlgebra::BlockVector       &dst,
                              const LinearAlgebra::BlockVector &src) const
    {
      Assert (src.n_blocks() == 2, ExcInternalError());
      Assert (dst.n_blocks() == 2, ExcInternalError());

      system_matrix.block(0,0).Tvmult(dst.block(0), src.block(0));
      system_matrix.block(1,0).Tvmult_add(dst.block(0), src.block(1));

      system_matrix.block(0,1).Tvmult(dst.block(1), src.block(0));
      system_matrix.block(1,1).Tvmult_add(dst.block(1), src.block(1));
    }


    void StokesBlock::vmult_add (LinearAlgebra::BlockVector       &dst,
                                 const LinearAlgebra::BlockVector &src) const
    {
      Assert (src.n_blocks() == 2, ExcInternalError());
      Assert (dst.n_blocks() == 2, ExcInternalError());

      system_matrix.block(0,0).vmult_add(dst.block(0), src.block(0));
      system_matrix.block(0,1).vmult_add(dst.block(0), src.block(1));

      system_matrix.block(1,0).vmult_add(dst.block(1), src.block(0));
      system_matrix.block(1,1).vmult_add(dst.block(1), src.block(1));
    }


    void StokesBlock::Tvmult_add (LinearAlgebra::BlockVector       &dst,
                                  const LinearAlgebra::BlockVector &src) const
    {
      Assert (src.n_blocks() == 2, ExcInternalError());
      Assert (dst.n_blocks() == 2, ExcInternalError());

      system_matrix.block(0,0).Tvmult_add(dst.block(0), src.block(0));
      system_matrix.block(1,0).Tvmult_add(dst.block(0), src.block(1));

      system_matrix.block(0,1).Tvmult_add(dst.block(1), src.block(0));
      system_matrix.block(1,1).Tvmult_add(dst.block(1), src.block(1));
    }



    double StokesBlock::residual (LinearAlgebra::BlockVector       &dst,
                                  const LinearAlgebra::BlockVector &x,
                                  const LinearAlgebra::BlockVector &b) const
    {
      Assert (x.n_blocks() == 2, ExcInternalError());
      Assert (dst.n_blocks() == 2, ExcInternalError());

      // compute b-Ax where A is only the top left 2x2 block
      this->vmult (dst, x);
      dst.block(0).sadd (-1, 1, b.block(0));
      dst.block(1).sadd (-1, 1, b.block(1));

      // clear blocks we didn't want to fill
      for (unsigned int block=2; block<dst.n_blocks(); ++block)
        dst.block(block) = 0;

      return dst.l2_norm();
    }


    /**
     * Implement the block Schur preconditioner for the Stokes system.
     */
    template <class PreconditionerA, class PreconditionerMp>
    class BlockSchurPreconditioner : public Subscriptor
    {
      public:
        /**
         * @brief Constructor
         *
         * @param S The entire Stokes matrix
         * @param Spre The matrix whose blocks are used in the definition of
         *     the preconditioning of the Stokes matrix, i.e. containing approximations
         *     of the A and S blocks.
         * @param Mppreconditioner Preconditioner object for the Schur complement,
         *     typically chosen as the mass matrix.
         * @param Apreconditioner Preconditioner object for the matrix A.
         * @param do_solve_A A flag indicating whether we should actually solve with
         *     the matrix $A$, or only apply one preconditioner step with it.
         **/
        BlockSchurPreconditioner (const LinearAlgebra::BlockSparseMatrix  &S,
                                  const LinearAlgebra::BlockSparseMatrix  &Spre,
                                  const PreconditionerMp                     &Mppreconditioner,
                                  const PreconditionerA                      &Apreconditioner,
                                  const bool                                  do_solve_A);

        /**
         * Matrix vector product with this preconditioner object.
         */
        void vmult (LinearAlgebra::BlockVector       &dst,
                    const LinearAlgebra::BlockVector &src) const;

      private:
        /**
         * References to the various matrix object this preconditioner works on.
         */
        const LinearAlgebra::BlockSparseMatrix &stokes_matrix;
        const LinearAlgebra::BlockSparseMatrix &stokes_preconditioner_matrix;
        const PreconditionerMp                    &mp_preconditioner;
        const PreconditionerA                     &a_preconditioner;

        /**
         * Whether to actually invert the $\tilde A$ part of the preconditioner matrix
         * or to just apply a single preconditioner step with it.
         **/
        const bool do_solve_A;
    };


    template <class PreconditionerA, class PreconditionerMp>
    BlockSchurPreconditioner<PreconditionerA, PreconditionerMp>::
    BlockSchurPreconditioner (const LinearAlgebra::BlockSparseMatrix  &S,
                              const LinearAlgebra::BlockSparseMatrix  &Spre,
                              const PreconditionerMp                     &Mppreconditioner,
                              const PreconditionerA                      &Apreconditioner,
                              const bool                                  do_solve_A)
      :
      stokes_matrix     (S),
      stokes_preconditioner_matrix     (Spre),
      mp_preconditioner (Mppreconditioner),
      a_preconditioner  (Apreconditioner),
      do_solve_A        (do_solve_A)
    {}


    template <class PreconditionerA, class PreconditionerMp>
    void
    BlockSchurPreconditioner<PreconditionerA, PreconditionerMp>::
    vmult (LinearAlgebra::BlockVector       &dst,
           const LinearAlgebra::BlockVector &src) const
    {
      LinearAlgebra::Vector utmp(src.block(0));

      {
        SolverControl solver_control(5000, 1e-6 * src.block(1).l2_norm());

#ifdef USE_PETSC
        SolverCG<LinearAlgebra::Vector> solver(solver_control);
#else
        TrilinosWrappers::SolverCG solver(solver_control);
#endif
        // Trilinos reports a breakdown
        // in case src=dst=0, even
        // though it should return
        // convergence without
        // iterating. We simply skip
        // solving in this case.
        if (src.block(1).l2_norm() > 1e-50 || dst.block(1).l2_norm() > 1e-50)
          solver.solve(stokes_preconditioner_matrix.block(1,1),
                       dst.block(1), src.block(1),
                       mp_preconditioner);

        dst.block(1) *= -1.0;
      }

      {
        stokes_matrix.block(0,1).vmult(utmp, dst.block(1)); //B^T
        utmp*=-1.0;
        utmp.add(src.block(0));
      }

      if (do_solve_A == true)
        {
          SolverControl solver_control(5000, utmp.l2_norm()*1e-2);
#ifdef USE_PETSC
          SolverCG<LinearAlgebra::Vector> solver(solver_control);
#else
          TrilinosWrappers::SolverCG solver(solver_control);
#endif
          solver.solve(stokes_matrix.block(0,0), dst.block(0), utmp,
                       a_preconditioner);
        }
      else
        a_preconditioner.vmult (dst.block(0), utmp);
    }

  }

  template <int dim>
  double Simulator<dim>::solve_advection (const TemperatureOrComposition &temperature_or_composition)
  {
    double advection_solver_tolerance = -1;
    unsigned int block_number = temperature_or_composition.block_index(introspection);

    if (temperature_or_composition.is_temperature())
      {
        computing_timer.enter_section ("   Solve temperature system");
        pcout << "   Solving temperature system... " << std::flush;
        advection_solver_tolerance = parameters.temperature_solver_tolerance;
      }
    else
      {
        computing_timer.enter_section ("   Solve composition system");
        pcout << "   Solving composition system "
              << temperature_or_composition.compositional_variable+1
              << "... " << std::flush;
        advection_solver_tolerance = parameters.composition_solver_tolerance;
      }

    const double tolerance = std::max(1e-50,
                                      advection_solver_tolerance*system_rhs.block(block_number).l2_norm());
    SolverControl solver_control (system_matrix.block(block_number, block_number).m(),
                                  tolerance);

    SolverGMRES<LinearAlgebra::Vector>   solver (solver_control,
                                                 SolverGMRES<LinearAlgebra::Vector>::AdditionalData(30,true));

    // Create distributed vector (we need all blocks here even though we only
    // solve for the current block) because only have a ConstraintMatrix
    // for the whole system, current_linearization_point contains our initial guess.
    LinearAlgebra::BlockVector distributed_solution (
      introspection.index_sets.system_partitioning,
      mpi_communicator);
    distributed_solution.block(block_number) = current_linearization_point.block (block_number);

    // Temporary vector to hold the residual, we don't need a BlockVector here.
    LinearAlgebra::Vector temp (
      introspection.index_sets.system_partitioning[block_number],
      mpi_communicator);

    // Compute the residual before we solve and return this at the end.
    // This is used in the nonlinear solver.
    const double initial_residual = system_matrix.block(block_number,block_number).residual
                                    (temp,
                                     distributed_solution.block(block_number),
                                     system_rhs.block(block_number));

    // solve the linear system:
    current_constraints.set_zero(distributed_solution);
    solver.solve (system_matrix.block(block_number,block_number),
                  distributed_solution.block(block_number),
                  system_rhs.block(block_number),
                  (temperature_or_composition.is_temperature()
                   ?
                   *T_preconditioner
                   :
                   *C_preconditioner));

    current_constraints.distribute (distributed_solution);
    solution.block(block_number) = distributed_solution.block(block_number);

    // print number of iterations and also record it in the
    // statistics file
    pcout << solver_control.last_step()
          << " iterations." << std::endl;

    if (temperature_or_composition.is_temperature())
      statistics.add_value("Iterations for temperature solver",
                           solver_control.last_step());
    else
      statistics.add_value("Iterations for composition solver " +
                           Utilities::int_to_string(temperature_or_composition.compositional_variable+1),
                           solver_control.last_step());

    computing_timer.exit_section();

    return initial_residual;
  }




  template <int dim>
  double Simulator<dim>::solve_stokes ()
  {
    computing_timer.enter_section ("   Solve Stokes system");

    pcout << "   Solving Stokes system... " << std::flush;

    const internal::StokesBlock stokes_block(system_matrix);

    // extract Stokes parts of solution vector, without any ghost elements
    LinearAlgebra::BlockVector distributed_stokes_solution (introspection.index_sets.stokes_partitioning, mpi_communicator);

    // create vector with distribution of system_rhs.
    LinearAlgebra::BlockVector remap (introspection.index_sets.stokes_partitioning, mpi_communicator);

    // copy current_linearization_point into it, because its distribution
    // is different.
    remap.block (0) = current_linearization_point.block (0);
    remap.block (1) = current_linearization_point.block (1);

    // before solving we scale the initial solution to the right dimensions
    denormalize_pressure (remap);
    current_constraints.set_zero (remap);
    remap.block (1) /= pressure_scaling;
    // if the model is compressible then we need to adjust the right hand
    // side of the equation to make it compatible with the matrix on the
    // left
    if (material_model->is_compressible ())
      make_pressure_rhs_compatible(system_rhs);

    // (ab)use the distributed solution vector to temporarily put a residual in
    // (we don't care about the residual vector -- all we care about is the
    // value (number) of the initial residual)
    const double initial_residual = stokes_block.residual (distributed_stokes_solution,
                                                           remap,
                                                           system_rhs);

    // then overwrite it again with the current best guess and solve the linear system
    distributed_stokes_solution = remap;

    // extract Stokes parts of rhs vector
    LinearAlgebra::BlockVector distributed_stokes_rhs(introspection.index_sets.stokes_partitioning);

    distributed_stokes_rhs.block(0) = system_rhs.block(0);
    distributed_stokes_rhs.block(1) = system_rhs.block(1);

    PrimitiveVectorMemory< LinearAlgebra::BlockVector > mem;

    // step 1a: try if the simple and fast solver
    // succeeds in 30 steps or less (or whatever the chosen value for the
    // corresponding parameter is).
    const double solver_tolerance = std::max (parameters.linear_stokes_solver_tolerance *
                                              distributed_stokes_rhs.l2_norm(),
                                              1e-12 * initial_residual);
    SolverControl solver_control_cheap (parameters.n_cheap_stokes_solver_steps,
                                        solver_tolerance);
    SolverControl solver_control_expensive (system_matrix.block(0,1).m() +
                                            system_matrix.block(1,0).m(), solver_tolerance);

    try
      {
        // if this cheaper solver is not desired, then simply short-cut
        // the attempt at solving with the cheaper preconditioner
        if (parameters.n_cheap_stokes_solver_steps == 0)
          throw SolverControl::NoConvergence(0,0);

        // otherwise give it a try with a preconditioner that consists
        // of only a single V-cycle
        const internal::BlockSchurPreconditioner<LinearAlgebra::PreconditionAMG,
              LinearAlgebra::PreconditionILU>
              preconditioner (system_matrix, system_preconditioner_matrix,
                              *Mp_preconditioner, *Amg_preconditioner,
                              false);

        SolverFGMRES<LinearAlgebra::BlockVector>
        solver(solver_control_cheap, mem,
               SolverFGMRES<LinearAlgebra::BlockVector>::
               AdditionalData(30, true));
        solver.solve(stokes_block, distributed_stokes_solution,
                     distributed_stokes_rhs, preconditioner);
      }

    // step 1b: take the stronger solver in case
    // the simple solver failed
    catch (SolverControl::NoConvergence)
      {
        const internal::BlockSchurPreconditioner<LinearAlgebra::PreconditionAMG,
              LinearAlgebra::PreconditionILU>
              preconditioner (system_matrix, system_preconditioner_matrix,
                              *Mp_preconditioner, *Amg_preconditioner,
                              true);

        SolverFGMRES<LinearAlgebra::BlockVector>
        solver(solver_control_expensive, mem,
               SolverFGMRES<LinearAlgebra::BlockVector>::
               AdditionalData(50, true));
        solver.solve(stokes_block, distributed_stokes_solution,
                     distributed_stokes_rhs, preconditioner);
      }

    // distribute hanging node and
    // other constraints
    current_constraints.distribute (distributed_stokes_solution);

    // now rescale the pressure back to real physical units
    distributed_stokes_solution.block(1) *= pressure_scaling;

    // then copy back the solution from the temporary (non-ghosted) vector
    // into the ghosted one with all solution components
    solution.block(0) = distributed_stokes_solution.block(0);
    solution.block(1) = distributed_stokes_solution.block(1);

    remove_nullspace(solution, distributed_stokes_solution);

    normalize_pressure(solution);

    // print the number of iterations to screen and record it in the
    // statistics file
    if (solver_control_expensive.last_step() == 0)
      pcout << solver_control_cheap.last_step()  << " iterations.";
    else
      pcout << solver_control_cheap.last_step() << '+'
            << solver_control_expensive.last_step() << " iterations.";
    pcout << std::endl;

    statistics.add_value("Iterations for Stokes solver",
                         solver_control_cheap.last_step() + solver_control_expensive.last_step());

    computing_timer.exit_section();

    return initial_residual;
  }

}





// explicit instantiation of the functions we implement in this file
namespace aspect
{
#define INSTANTIATE(dim) \
  template double Simulator<dim>::solve_advection (const TemperatureOrComposition &); \
  template double Simulator<dim>::solve_stokes ();

  ASPECT_INSTANTIATE(INSTANTIATE)
}
