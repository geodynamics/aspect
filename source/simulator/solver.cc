/*
  Copyright (C) 2011 - 2020 by the authors of the ASPECT code.

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
#include <aspect/stokes_matrix_free.h>

#include <deal.II/base/signaling_nan.h>
#include <deal.II/lac/solver_gmres.h>

#ifdef ASPECT_USE_PETSC
#include <deal.II/lac/solver_cg.h>
#else
#include <deal.II/lac/trilinos_solver.h>
#endif

#include <deal.II/fe/fe_values.h>

namespace aspect
{
  namespace internal
  {
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
         * @param A_block_tolerance The tolerance for the CG solver which computes
         *     the inverse of the A block.
         * @param S_block_tolerance The tolerance for the CG solver which computes
         *     the inverse of the S block (Schur complement matrix).
         **/
        BlockSchurPreconditioner (const LinearAlgebra::BlockSparseMatrix  &S,
                                  const LinearAlgebra::BlockSparseMatrix  &Spre,
                                  const PreconditionerMp                     &Mppreconditioner,
                                  const PreconditionerA                      &Apreconditioner,
                                  const bool                                  do_solve_A,
                                  const double                                A_block_tolerance,
                                  const double                                S_block_tolerance);

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
        const double A_block_tolerance;
        const double S_block_tolerance;
    };


    template <class PreconditionerA, class PreconditionerMp>
    BlockSchurPreconditioner<PreconditionerA, PreconditionerMp>::
    BlockSchurPreconditioner (const LinearAlgebra::BlockSparseMatrix  &S,
                              const LinearAlgebra::BlockSparseMatrix  &Spre,
                              const PreconditionerMp                     &Mppreconditioner,
                              const PreconditionerA                      &Apreconditioner,
                              const bool                                  do_solve_A,
                              const double                                A_block_tolerance,
                              const double                                S_block_tolerance)
      :
      stokes_matrix     (S),
      stokes_preconditioner_matrix     (Spre),
      mp_preconditioner (Mppreconditioner),
      a_preconditioner  (Apreconditioner),
      do_solve_A        (do_solve_A),
      n_iterations_A_(0),
      n_iterations_S_(0),
      A_block_tolerance(A_block_tolerance),
      S_block_tolerance(S_block_tolerance)
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
        SolverControl solver_control(1000, src.block(1).l2_norm() * S_block_tolerance);

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
                                                       "did not converge to a tolerance of "
                                                       + Utilities::to_string(solver_control.tolerance()) +
                                                       ". It reported the following error:\n\n")
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
        stokes_matrix.block(0,1).vmult(utmp, dst.block(1)); // B^T or J^{up}
        utmp *= -1.0;
        utmp += src.block(0);
      }

      // now either solve with the top left block (if do_solve_A==true)
      // or just apply one preconditioner sweep (for the first few
      // iterations of our two-stage outer GMRES iteration)
      if (do_solve_A == true)
        {
          SolverControl solver_control(10000, utmp.l2_norm() * A_block_tolerance);
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
                                                     "did not converge to a tolerance of "
                                                     + Utilities::to_string(solver_control.tolerance()) +
                                                     ". It reported the following error:\n\n")
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

  namespace
  {
    void linear_solver_failed(const std::string &solver_name,
                              const std::string &output_filename,
                              const std::vector<SolverControl> &solver_controls,
                              const std::exception &exc)
    {
      // output solver history
      std::ofstream f((output_filename).c_str());

      for (unsigned int i=0; i<solver_controls.size(); ++i)
        {
          if (i>0)
            f << "\n";

          // Only request the solver history if a history has actually been created
          for (unsigned int j=0; j<solver_controls[i].get_history_data().size(); ++j)
            f << j << " " << solver_controls[i].get_history_data()[j] << "\n";
        }

      f.close();

      AssertThrow (false,
                   ExcMessage ("The " + solver_name
                               + " did not converge. It reported the following error:\n\n"
                               +
                               exc.what()
                               + "\n The required residual for convergence is: " + std::to_string(solver_controls.front().tolerance())
                               + ".\n See " + output_filename
                               + " for convergence history."));
    }
  }

  template <int dim>
  double Simulator<dim>::solve_advection (const AdvectionField &advection_field)
  {
    double advection_solver_tolerance = -1;
    unsigned int block_idx = advection_field.block_index(introspection);

    std::string field_name = (advection_field.is_temperature()
                              ?
                              "temperature"
                              :
                              introspection.name_for_compositional_index(advection_field.compositional_variable) + " composition");

    if (advection_field.is_temperature())
      advection_solver_tolerance = parameters.temperature_solver_tolerance;
    else
      advection_solver_tolerance = parameters.composition_solver_tolerance;

    const double tolerance = std::max(1e-50,
                                      advection_solver_tolerance*system_rhs.block(block_idx).l2_norm());

    SolverControl solver_control (1000, tolerance);

    solver_control.enable_history_data();

    SolverGMRES<LinearAlgebra::Vector>   solver (solver_control,
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

    TimerOutput::Scope timer (computing_timer, (advection_field.is_temperature() ?
                                                "Solve temperature system" :
                                                "Solve composition system"));
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
    // solve for the current block) because only have a AffineConstraints<double>
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

        if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
          {
            linear_solver_failed("iterative advection solver",
                                 parameters.output_directory+"solver_history.txt",
                                 std::vector<SolverControl> {solver_control},
                                 exc);
          }
        else
          throw QuietException();
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

    if ((advection_field.is_temperature()
         && parameters.use_discontinuous_temperature_discretization
         && parameters.use_limiter_for_discontinuous_temperature_solution)
        ||
        (!advection_field.is_temperature()
         && parameters.use_discontinuous_composition_discretization
         && parameters.use_limiter_for_discontinuous_composition_solution))
      apply_limiter_to_dg_solutions(advection_field);

    return initial_residual;
  }



  template <int dim>
  std::pair<double,double>
  Simulator<dim>::solve_stokes ()
  {
    TimerOutput::Scope timer (computing_timer, "Solve Stokes system");
    pcout << "   Solving Stokes system... " << std::flush;

    if (parameters.stokes_solver_type == Parameters<dim>::StokesSolverType::block_gmg)
      {
        return stokes_matrix_free->solve();
      }

    // extract Stokes parts of solution vector, without any ghost elements
    LinearAlgebra::BlockVector distributed_stokes_solution (introspection.index_sets.stokes_partitioning, mpi_communicator);

    double initial_nonlinear_residual = numbers::signaling_nan<double>();
    double final_linear_residual      = numbers::signaling_nan<double>();

    if (parameters.use_direct_stokes_solver)
      {
        // We hard-code the blocks down below, so make sure block 0 is indeed
        // the block containing velocity and pressure:
        Assert(introspection.block_indices.velocities == 0, ExcNotImplemented());
        Assert(introspection.block_indices.pressure == 0
               ||
               (parameters.include_melt_transport
                && introspection.variable("fluid pressure").block_index == 0
                && introspection.variable("compaction pressure").block_index == 0),
               ExcNotImplemented());

        // start with a reasonable guess
        solution.block(0) = current_linearization_point.block(0);

        // While we don't need to set up the initial guess for the direct solver
        // (it will be ignored by the solver anyway), we need this if we are
        // using a nonlinear scheme, because we use this to compute the current
        // nonlinear residual (see initial_residual below).

        // TODO: if there was an easy way to know if the caller needs the
        // initial residual we could skip all of this stuff.
        distributed_stokes_solution.block(0) = solution.block(0);
        denormalize_pressure (this->last_pressure_normalization_adjustment,
                              distributed_stokes_solution,
                              solution);
        current_constraints.set_zero (distributed_stokes_solution);

        // Undo the pressure scaling:
        IndexSet &pressure_idxset = parameters.include_melt_transport ?
                                    introspection.index_sets.locally_owned_melt_pressure_dofs
                                    : introspection.index_sets.locally_owned_pressure_dofs;

        for (unsigned int i=0; i< pressure_idxset.n_elements(); ++i)
          {
            types::global_dof_index idx = pressure_idxset.nth_index_in_set(i);

            distributed_stokes_solution(idx) /= pressure_scaling;
          }
        distributed_stokes_solution.compress(VectorOperation::insert);

        // we need a temporary vector for the residual (even if we don't care about it)
        LinearAlgebra::Vector residual (introspection.index_sets.stokes_partitioning[0], mpi_communicator);

        initial_nonlinear_residual = system_matrix.block(0,0).residual(
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
        {
          IndexSet &pressure_idxset = parameters.include_melt_transport ?
                                      introspection.index_sets.locally_owned_melt_pressure_dofs
                                      : introspection.index_sets.locally_owned_pressure_dofs;

          for (unsigned int i=0; i< pressure_idxset.n_elements(); ++i)
            {
              types::global_dof_index idx = pressure_idxset.nth_index_in_set(i);
              distributed_stokes_solution(idx) *= pressure_scaling;
            }
          distributed_stokes_solution.compress(VectorOperation::insert);
        }

        // then copy back the solution from the temporary (non-ghosted) vector
        // into the ghosted one with all solution components
        solution.block(0) = distributed_stokes_solution.block(0);

        pcout << "done." << std::endl;
      }
    else
      {
        // Many parts of the solver depend on the block layout (velocity = 0,
        // pressure = 1). For example the linearized_stokes_initial_guess vector or the StokesBlock matrix
        // wrapper. Let us make sure that this holds (and shorten their names):
        const unsigned int block_vel = introspection.block_indices.velocities;
        const unsigned int block_p = (parameters.include_melt_transport) ?
                                     introspection.variable("fluid pressure").block_index
                                     : introspection.block_indices.pressure;
        Assert(block_vel == 0, ExcNotImplemented());
        Assert(block_p == 1, ExcNotImplemented());
        Assert(!parameters.include_melt_transport
               || introspection.variable("compaction pressure").block_index == 1, ExcNotImplemented());

        const internal::StokesBlock stokes_block(system_matrix);

        // create a completely distributed vector that will be used for
        // the scaled and denormalized solution and later used as a
        // starting guess for the linear solver
        LinearAlgebra::BlockVector linearized_stokes_initial_guess (introspection.index_sets.stokes_partitioning, mpi_communicator);

        // copy the velocity and pressure from current_linearization_point into
        // the vector linearized_stokes_initial_guess. We need to do the copy because
        // linearized_stokes_variables has a different
        // layout than current_linearization_point, which also contains all the
        // other solution variables.
        if (assemble_newton_stokes_system == false)
          {
            linearized_stokes_initial_guess.block (block_vel) = current_linearization_point.block (block_vel);
            linearized_stokes_initial_guess.block (block_p) = current_linearization_point.block (block_p);

            denormalize_pressure (this->last_pressure_normalization_adjustment,
                                  linearized_stokes_initial_guess,
                                  current_linearization_point);
          }
        else
          {
            // The Newton solver solves for updates to variables, for which our best guess is zero when
            // it isn't the first nonlinear iteration. When it is the first nonlinear iteration, we
            // have to assemble the full (non-defect correction) Picard, to get the boundary conditions
            // right in combination with being able to use the initial guess optimally. So we should never
            // end up here when it is the first nonlinear iteration.
            Assert(nonlinear_iteration != 0,
                   ExcMessage ("The Newton solver should not be active in the first nonlinear iteration."));

            linearized_stokes_initial_guess.block (block_vel) = 0;
            linearized_stokes_initial_guess.block (block_p) = 0;
          }

        current_constraints.set_zero (linearized_stokes_initial_guess);
        linearized_stokes_initial_guess.block (block_p) /= pressure_scaling;

        double solver_tolerance = 0;
        if (assemble_newton_stokes_system == false)
          {
            // (ab)use the distributed solution vector to temporarily put a residual in
            // (we don't care about the residual vector -- all we care about is the
            // value (number) of the initial residual). The initial residual is returned
            // to the caller (for nonlinear computations). This value is computed before
            // the solve because we want to compute || A^{k+1} U^k - F^{k+1} ||, which is
            // the nonlinear residual. Because the place where the nonlinear residual is
            // checked against the nonlinear tolerance comes after the solve, the system
            // is solved one time too many in the case of a nonlinear Picard solver.
            initial_nonlinear_residual = stokes_block.residual (distributed_stokes_solution,
                                                                linearized_stokes_initial_guess,
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
                                                                         linearized_stokes_initial_guess.block(1),
                                                                         system_rhs.block(0));
            const double residual_p = system_rhs.block(1).l2_norm();

            solver_tolerance = parameters.linear_stokes_solver_tolerance *
                               std::sqrt(residual_u*residual_u+residual_p*residual_p);
          }
        else
          {
            // if we are solving for the Newton update, then the initial guess of the solution
            // vector is the zero vector, and the starting (nonlinear) residual is simply
            // the norm of the (Newton) right hand side vector
            const double residual_u = system_rhs.block(0).l2_norm();
            const double residual_p = system_rhs.block(1).l2_norm();
            solver_tolerance = parameters.linear_stokes_solver_tolerance *
                               std::sqrt(residual_u*residual_u+residual_p*residual_p);

            // as described in the documentation of the function, the initial
            // nonlinear residual for the Newton method is computed by just
            // taking the norm of the right hand side
            initial_nonlinear_residual = std::sqrt(residual_u*residual_u+residual_p*residual_p);
          }
        // Now overwrite the solution vector again with the current best guess
        // to solve the linear system
        distributed_stokes_solution = linearized_stokes_initial_guess;

        // extract Stokes parts of rhs vector
        LinearAlgebra::BlockVector distributed_stokes_rhs(introspection.index_sets.stokes_partitioning);

        distributed_stokes_rhs.block(block_vel) = system_rhs.block(block_vel);
        distributed_stokes_rhs.block(block_p) = system_rhs.block(block_p);

        PrimitiveVectorMemory< LinearAlgebra::BlockVector > mem;

        // create Solver controls for the cheap and expensive solver phase
        SolverControl solver_control_cheap (parameters.n_cheap_stokes_solver_steps,
                                            solver_tolerance);

        SolverControl solver_control_expensive (parameters.n_expensive_stokes_solver_steps,
                                                solver_tolerance);

        solver_control_cheap.enable_history_data();
        solver_control_expensive.enable_history_data();

        // create a cheap preconditioner that consists of only a single V-cycle
        const internal::BlockSchurPreconditioner<LinearAlgebra::PreconditionAMG,
              LinearAlgebra::PreconditionBase>
              preconditioner_cheap (system_matrix, system_preconditioner_matrix,
                                    *Mp_preconditioner, *Amg_preconditioner,
                                    false,
                                    parameters.linear_solver_A_block_tolerance,
                                    parameters.linear_solver_S_block_tolerance);

        // create an expensive preconditioner that solves for the A block with CG
        const internal::BlockSchurPreconditioner<LinearAlgebra::PreconditionAMG,
              LinearAlgebra::PreconditionBase>
              preconditioner_expensive (system_matrix, system_preconditioner_matrix,
                                        *Mp_preconditioner, *Amg_preconditioner,
                                        true,
                                        parameters.linear_solver_A_block_tolerance,
                                        parameters.linear_solver_S_block_tolerance);

        // step 1a: try if the simple and fast solver
        // succeeds in n_cheap_stokes_solver_steps steps or less.
        try
          {
            // if this cheaper solver is not desired, then simply short-cut
            // the attempt at solving with the cheaper preconditioner
            if (parameters.n_cheap_stokes_solver_steps == 0)
              throw SolverControl::NoConvergence(0,0);

            SolverFGMRES<LinearAlgebra::BlockVector>
            solver(solver_control_cheap, mem,
                   SolverFGMRES<LinearAlgebra::BlockVector>::
                   AdditionalData(parameters.stokes_gmres_restart_length));

            solver.solve (stokes_block,
                          distributed_stokes_solution,
                          distributed_stokes_rhs,
                          preconditioner_cheap);

            final_linear_residual = solver_control_cheap.last_value();
          }

        // step 1b: take the stronger solver in case
        // the simple solver failed and attempt solving
        // it in n_expensive_stokes_solver_steps steps or less.
        catch (const SolverControl::NoConvergence &)
          {
            // use the value defined by the user
            // OR
            // at least a restart length of 100 for melt models
            const unsigned int number_of_temporary_vectors = (parameters.include_melt_transport == false ?
                                                              parameters.stokes_gmres_restart_length :
                                                              std::max(parameters.stokes_gmres_restart_length, 100U));

            SolverFGMRES<LinearAlgebra::BlockVector>
            solver(solver_control_expensive, mem,
                   SolverFGMRES<LinearAlgebra::BlockVector>::
                   AdditionalData(number_of_temporary_vectors));

            try
              {
                AssertThrow (parameters.n_expensive_stokes_solver_steps>0,
                             ExcMessage ("The Stokes solver did not converge in the number of requested cheap iterations and "
                                         "you requested 0 for ``Maximum number of expensive Stokes solver steps''. Aborting."));

                solver.solve(stokes_block,
                             distributed_stokes_solution,
                             distributed_stokes_rhs,
                             preconditioner_expensive);

                final_linear_residual = solver_control_expensive.last_value();
              }
            // if the solver fails, report the error from processor 0 with some additional
            // information about its location, and throw a quiet exception on all other
            // processors
            catch (const std::exception &exc)
              {
                signals.post_stokes_solver(*this,
                                           preconditioner_cheap.n_iterations_S() + preconditioner_expensive.n_iterations_S(),
                                           preconditioner_cheap.n_iterations_A() + preconditioner_expensive.n_iterations_A(),
                                           solver_control_cheap,
                                           solver_control_expensive);

                if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
                  {
                    linear_solver_failed("iterative Stokes solver",
                                         parameters.output_directory+"solver_history.txt",
                                         parameters.n_cheap_stokes_solver_steps > 0 ?
                                         std::vector<SolverControl> {solver_control_cheap, solver_control_expensive} :
                                         std::vector<SolverControl> {solver_control_expensive},
                                         exc);
                  }
                else
                  {
                    throw QuietException();
                  }
              }
          }

        // signal successful solver
        signals.post_stokes_solver(*this,
                                   preconditioner_cheap.n_iterations_S() + preconditioner_expensive.n_iterations_S(),
                                   preconditioner_cheap.n_iterations_A() + preconditioner_expensive.n_iterations_A(),
                                   solver_control_cheap,
                                   solver_control_expensive);

        // distribute hanging node and
        // other constraints
        current_constraints.distribute (distributed_stokes_solution);

        // now rescale the pressure back to real physical units
        distributed_stokes_solution.block(block_p) *= pressure_scaling;

        // then copy back the solution from the temporary (non-ghosted) vector
        // into the ghosted one with all solution components
        solution.block(block_vel) = distributed_stokes_solution.block(block_vel);
        solution.block(block_p) = distributed_stokes_solution.block(block_p);

        // print the number of iterations to screen
        pcout << (solver_control_cheap.last_step() != numbers::invalid_unsigned_int ?
                  solver_control_cheap.last_step():
                  0)
              << '+'
              << (solver_control_expensive.last_step() != numbers::invalid_unsigned_int ?
                  solver_control_expensive.last_step():
                  0)
              << " iterations.";
        pcout << std::endl;
      }


    // do some cleanup now that we have the solution
    remove_nullspace(solution, distributed_stokes_solution);
    if (assemble_newton_stokes_system == false)
      this->last_pressure_normalization_adjustment = normalize_pressure(solution);

    // convert melt pressures:
    if (parameters.include_melt_transport)
      melt_handler->compute_melt_variables(system_matrix,solution,system_rhs);

    return std::pair<double,double>(initial_nonlinear_residual,
                                    final_linear_residual);
  }

}





// explicit instantiation of the functions we implement in this file
namespace aspect
{
#define INSTANTIATE(dim) \
  template double Simulator<dim>::solve_advection (const AdvectionField &); \
  template std::pair<double,double> Simulator<dim>::solve_stokes ();

  ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
}
