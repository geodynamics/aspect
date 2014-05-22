/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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


#include <aspect/simulator.h>
#include <aspect/global.h>

#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/constraint_matrix.h>

#ifdef ASPECT_USE_PETSC
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

        unsigned int n_iterations_A() const;
        unsigned int n_iterations_S() const;

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
        mutable unsigned int n_iterations_A_;
        mutable unsigned int n_iterations_S_;
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
      do_solve_A        (do_solve_A),
      n_iterations_A_(0),
      n_iterations_S_(0)
    {}

    template <class PreconditionerA, class PreconditionerMp>
    unsigned int
    BlockSchurPreconditioner<PreconditionerA, PreconditionerMp>::
    n_iterations_A() const
    {
      return n_iterations_A_;
    }
    template <class PreconditionerA, class PreconditionerMp>
    unsigned int
    BlockSchurPreconditioner<PreconditionerA, PreconditionerMp>::
    n_iterations_S() const
    {
      return n_iterations_S_;
    }


    template <class PreconditionerA, class PreconditionerMp>
    void
    BlockSchurPreconditioner<PreconditionerA, PreconditionerMp>::
    vmult (LinearAlgebra::BlockVector       &dst,
           const LinearAlgebra::BlockVector &src) const
    {
      LinearAlgebra::Vector utmp(src.block(0));

      // first solve with the bottom left block, which we have built
      // as a mass matrix with the inverse of the viscosity
      {
        SolverControl solver_control(1000, 1e-6 * src.block(1).l2_norm());

#ifdef ASPECT_USE_PETSC
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
        if (src.block(1).l2_norm() > 1e-50)
          {
            try
              {
                dst.block(1) = 0.0;
                solver.solve(stokes_preconditioner_matrix.block(1,1),
                             dst.block(1), src.block(1),
                             mp_preconditioner);
                n_iterations_S_ += solver_control.last_step();
              }
            // if the solver fails, report the error from processor 0 with some additional
            // information about its location, and throw a quiet exception on all other
            // processors
            catch (const std::exception &exc)
              {
                if (Utilities::MPI::this_mpi_process(src.block(0).get_mpi_communicator()) == 0)
                  AssertThrow (false,
                               ExcMessage (std::string("The iterative (bottom right) solver in BlockSchurPreconditioner::vmult "
                                                       "did not converge. It reported the following error:\n\n")
                                           +
                                           exc.what()))
                  else
                    throw QuietException();
              }
          }

        dst.block(1) *= -1.0;
      }

      // apply the top right block
      {
        stokes_matrix.block(0,1).vmult(utmp, dst.block(1)); //B^T
        utmp*=-1.0;
        utmp.add(src.block(0));
      }

      // now either solve with the top left block (if do_solve_A==true)
      // or just apply one preconditioner sweep (for the first few
      // iterations of our two-stage outer GMRES iteration)
      if (do_solve_A == true)
        {
          SolverControl solver_control(10000, utmp.l2_norm()*1e-2);
#ifdef ASPECT_USE_PETSC
          SolverCG<LinearAlgebra::Vector> solver(solver_control);
#else
          TrilinosWrappers::SolverCG solver(solver_control);
#endif
          try
            {
              dst.block(0) = 0.0;
              solver.solve(stokes_matrix.block(0,0), dst.block(0), utmp,
                           a_preconditioner);
              n_iterations_A_ += solver_control.last_step();
            }
          // if the solver fails, report the error from processor 0 with some additional
          // information about its location, and throw a quiet exception on all other
          // processors
          catch (const std::exception &exc)
            {
              if (Utilities::MPI::this_mpi_process(src.block(0).get_mpi_communicator()) == 0)
                AssertThrow (false,
                             ExcMessage (std::string("The iterative (top left) solver in BlockSchurPreconditioner::vmult "
                                                     "did not converge. It reported the following error:\n\n")
                                         +
                                         exc.what()))
                else
                  throw QuietException();
            }
        }
      else
        {
          a_preconditioner.vmult (dst.block(0), utmp);
          n_iterations_A_ += 1;
        }
    }

  }

  template <int dim>
  double Simulator<dim>::solve_advection (const AdvectionField &advection_field)
  {
    double advection_solver_tolerance = -1;
    unsigned int block_idx = advection_field.block_index(introspection);

    if (advection_field.is_temperature())
      {
        computing_timer.enter_section ("   Solve temperature system");
        pcout << "   Solving temperature system... " << std::flush;
        advection_solver_tolerance = parameters.temperature_solver_tolerance;
      }
    else
      {
        computing_timer.enter_section ("   Solve composition system");
        pcout << "   Solving "
              << introspection.name_for_compositional_index(advection_field.compositional_variable)
              << " system "
              << "... " << std::flush;
        advection_solver_tolerance = parameters.composition_solver_tolerance;
      }

    const double tolerance = std::max(1e-50,
                                      advection_solver_tolerance*system_rhs.block(block_idx).l2_norm());
    SolverControl solver_control (1000, tolerance);

    SolverGMRES<LinearAlgebra::Vector>   solver (solver_control,
                                                 SolverGMRES<LinearAlgebra::Vector>::AdditionalData(30,true));

    // Create distributed vector (we need all blocks here even though we only
    // solve for the current block) because only have a ConstraintMatrix
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
        solver.solve (system_matrix.block(block_idx,block_idx),
                      distributed_solution.block(block_idx),
                      system_rhs.block(block_idx),
                      (advection_field.is_temperature()
                       ?
                       *T_preconditioner
                       :
                       *C_preconditioner));
      }
    // if the solver fails, report the error from processor 0 with some additional
    // information about its location, and throw a quiet exception on all other
    // processors
    catch (const std::exception &exc)
      {
        if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
          AssertThrow (false,
                       ExcMessage (std::string("The iterative advection solver "
                                               "did not converge. It reported the following error:\n\n")
                                   +
                                   exc.what()))
          else
            throw QuietException();
      }

    current_constraints.distribute (distributed_solution);
    solution.block(block_idx) = distributed_solution.block(block_idx);

    // print number of iterations and also record it in the
    // statistics file
    pcout << solver_control.last_step()
          << " iterations." << std::endl;

    if (advection_field.is_temperature())
      statistics.add_value("Iterations for temperature solver",
                           solver_control.last_step());
    else
      statistics.add_value("Iterations for composition solver " +
                           Utilities::int_to_string(advection_field.compositional_variable+1),
                           solver_control.last_step());

    computing_timer.exit_section();

    return initial_residual;
  }




  template <int dim>
  double Simulator<dim>::solve_stokes ()
  {
    computing_timer.enter_section ("   Solve Stokes system");
    pcout << "   Solving Stokes system... " << std::flush;

    if (parameters.use_direct_stokes_solver)
      {
        // We hardcode the blocks down below, so make sure block 0 is indeed
        // the block containing velocity and pressure:
        Assert(introspection.block_indices.velocities == 0, ExcNotImplemented());
        Assert(introspection.block_indices.pressure == 0, ExcNotImplemented());

        LinearAlgebra::BlockVector distributed_stokes_solution (introspection.index_sets.stokes_partitioning, mpi_communicator);

        // While we don't need to set up the initial guess for the direct solver
        // (it will be ignored by the solver anyway), we need this if we are
        // using a nonlinear scheme, because we use this to compute the current
        // nonlinear residual (see initial_residual below).
        // TODO: if there was an easy way to know if the caller needs the
        // initial residual we could skip all of this stuff.
        solution.block(0) = current_linearization_point.block(0);
        denormalize_pressure (solution);
        distributed_stokes_solution.block(0) = solution.block(0);
        current_constraints.set_zero (distributed_stokes_solution);

        // Undo the pressure scaling:
        for (unsigned int i=0; i< introspection.index_sets.locally_owned_pressure_dofs.n_elements(); ++i)
          {
            types::global_dof_index idx = introspection.index_sets.locally_owned_pressure_dofs.nth_index_in_set(i);

            distributed_stokes_solution(idx) /= pressure_scaling;
          }
        distributed_stokes_solution.compress(VectorOperation::insert);

        if (do_pressure_rhs_compatibility_modification)
          make_pressure_rhs_compatible(system_rhs);

        // we need a temporary vector for the residual (even if we don't care about it)
        LinearAlgebra::Vector residual (introspection.index_sets.stokes_partitioning[0], mpi_communicator);

        const double initial_residual = system_matrix.block(0,0).residual(
                                          residual,
                                          distributed_stokes_solution.block(0),
                                          system_rhs.block(0));

        SolverControl cn;
        // TODO: can we re-use the direct solver?
#ifdef ASPECT_USE_PETSC
        PETScWrappers::SparseDirectMUMPS solver(cn, mpi_communicator);
#else
        TrilinosWrappers::SolverDirect solver(cn);
#endif
        try
          {
            solver.solve(system_matrix.block(0,0),
                         distributed_stokes_solution.block(0),
                         system_rhs.block(0));
          }
        // if the solver fails, report the error from processor 0 with some additional
        // information about its location, and throw a quiet exception on all other
        // processors
        catch (const std::exception &exc)
          {
            if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
              AssertThrow (false,
                           ExcMessage (std::string("The direct Stokes solver "
                                                   "did not succeed. It reported the following error:\n\n")
                                       +
                                       exc.what()))
              else
                throw QuietException();
          }


        current_constraints.distribute (distributed_stokes_solution);

        // now rescale the pressure back to real physical units:
        for (unsigned int i=0; i< introspection.index_sets.locally_owned_pressure_dofs.n_elements(); ++i)
          {
            types::global_dof_index idx = introspection.index_sets.locally_owned_pressure_dofs.nth_index_in_set(i);

            distributed_stokes_solution(idx) *= pressure_scaling;
          }
        distributed_stokes_solution.compress(VectorOperation::insert);

        // then copy back the solution from the temporary (non-ghosted) vector
        // into the ghosted one with all solution components
        solution.block(0) = distributed_stokes_solution.block(0);

        remove_nullspace(solution, distributed_stokes_solution);

        normalize_pressure(solution);

        pcout << "done." << std::endl;

        computing_timer.exit_section();

        return initial_residual;
      }


    // Many parts of the solver depend on the block layout (velocity = 0,
    // pressure = 1). For example the remap vector or the StokesBlock matrix
    // wrapper. Let us make sure that this holds (and shorten their names):
    const unsigned int block_vel = introspection.block_indices.velocities;
    const unsigned int block_p = introspection.block_indices.pressure;
    Assert(block_vel == 0, ExcNotImplemented());
    Assert(block_p == 1, ExcNotImplemented());

    const internal::StokesBlock stokes_block(system_matrix);

    // extract Stokes parts of solution vector, without any ghost elements
    LinearAlgebra::BlockVector distributed_stokes_solution (introspection.index_sets.stokes_partitioning, mpi_communicator);

    // create vector with distribution of system_rhs.
    LinearAlgebra::BlockVector remap (introspection.index_sets.stokes_partitioning, mpi_communicator);

    // copy the velocity and pressure from current_linearization_point into
    // the vector remap. We need to do the copy because remap has a different
    // layout than current_linearization_point, which also contains all the
    // other solution variables.
    remap.block (block_vel) = current_linearization_point.block (block_vel);
    remap.block (block_p) = current_linearization_point.block (block_p);

    // before solving we scale the initial solution to the right dimensions
    denormalize_pressure (remap);
    current_constraints.set_zero (remap);
    remap.block (block_p) /= pressure_scaling;
    // if the model is compressible then we need to adjust the right hand
    // side of the equation to make it compatible with the matrix on the
    // left
    if (do_pressure_rhs_compatibility_modification)
      make_pressure_rhs_compatible(system_rhs);

    // (ab)use the distributed solution vector to temporarily put a residual in
    // (we don't care about the residual vector -- all we care about is the
    // value (number) of the initial residual). The initial residual is returned
    // to the caller (for nonlinear computations).
    const double initial_residual = stokes_block.residual (distributed_stokes_solution,
                                                           remap,
                                                           system_rhs);

    // Note: the residual is computed with a zero velocity, effectively computing
    // || B^T p - g ||, which we are going to use for our solver tolerance.
    // We do not use the current velocity for the initial residual because
    // this would not decrease the number of iterations if we had a better
    // initial guess (say using a smaller timestep). But we need to use
    // the pressure instead of only using the norm of the rhs, because we
    // are only interested in the part of the rhs not balanced by the static
    // pressure (the current pressure is a good approximation for the static
    // pressure).
    const double residual_u = system_matrix.block(0,1).residual (distributed_stokes_solution.block(0),
                                                                 remap.block(1),
                                                                 system_rhs.block(0));
    const double residual_p = system_rhs.block(1).l2_norm();
    const double solver_tolerance = parameters.linear_stokes_solver_tolerance *
                                    sqrt(residual_u*residual_u+residual_p*residual_p);

    // Now overwrite the solution vector again with the current best guess
    // to solve the linear system
    distributed_stokes_solution = remap;

    // extract Stokes parts of rhs vector
    LinearAlgebra::BlockVector distributed_stokes_rhs(introspection.index_sets.stokes_partitioning);

    distributed_stokes_rhs.block(block_vel) = system_rhs.block(block_vel);
    distributed_stokes_rhs.block(block_p) = system_rhs.block(block_p);

    PrimitiveVectorMemory< LinearAlgebra::BlockVector > mem;

    // step 1a: try if the simple and fast solver
    // succeeds in 30 steps or less (or whatever the chosen value for the
    // corresponding parameter is).
    SolverControl solver_control_cheap (parameters.n_cheap_stokes_solver_steps,
                                        solver_tolerance);
    SolverControl solver_control_expensive (system_matrix.block(block_vel,block_p).m() +
                                            system_matrix.block(block_p,block_vel).m(), solver_tolerance);

    unsigned int its_A = 0, its_S = 0;
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
        solver.solve (stokes_block,
                      distributed_stokes_solution,
                      distributed_stokes_rhs,
                      preconditioner);

        its_A += preconditioner.n_iterations_A();
        its_S += preconditioner.n_iterations_S();
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

        try
          {
            solver.solve(stokes_block,
                         distributed_stokes_solution,
                         distributed_stokes_rhs,
                         preconditioner);
          }
        // if the solver fails, report the error from processor 0 with some additional
        // information about its location, and throw a quiet exception on all other
        // processors
        catch (const std::exception &exc)
          {
            if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
              AssertThrow (false,
                           ExcMessage (std::string("The iterative Stokes solver "
                                                   "did not converge. It reported the following error:\n\n")
                                       +
                                       exc.what()))
              else
                throw QuietException();
          }


        its_A += preconditioner.n_iterations_A();
        its_S += preconditioner.n_iterations_S();
      }

    // distribute hanging node and
    // other constraints
    current_constraints.distribute (distributed_stokes_solution);

    // now rescale the pressure back to real physical units
    distributed_stokes_solution.block(block_p) *= pressure_scaling;

    // then copy back the solution from the temporary (non-ghosted) vector
    // into the ghosted one with all solution components
    solution.block(block_vel) = distributed_stokes_solution.block(block_vel);
    solution.block(block_p) = distributed_stokes_solution.block(block_p);

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
    statistics.add_value("Velocity iterations in Stokes preconditioner",
                         its_A);
    statistics.add_value("Schur complement iterations in Stokes preconditioner",
                         its_S);

    computing_timer.exit_section();

    return initial_residual;
  }

}





// explicit instantiation of the functions we implement in this file
namespace aspect
{
#define INSTANTIATE(dim) \
  template double Simulator<dim>::solve_advection (const AdvectionField &); \
  template double Simulator<dim>::solve_stokes ();

  ASPECT_INSTANTIATE(INSTANTIATE)
}
