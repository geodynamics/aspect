#include <aspect/simulator.h>
#include <aspect/global.h>

#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/pointer_matrix.h>


namespace aspect
{
  namespace internal
  {
    using namespace dealii;

    /**
     * Implement multiplication with Stokes part of system matrix
     */
    class StokesBlock : public PointerMatrixBase<LinearAlgebra::BlockVector>
    {
      public:
        /**
         * @brief Constructor
         *
         * @param S The entire system matrix
         */
        StokesBlock (const LinearAlgebra::BlockSparseMatrix  &S)
          : system_matrix(S) {};

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

        void clear() {};

      private:

        const void *get () const
        {
          return &system_matrix;
        };

        /**
         * References to the system matrix object.
         */
        const LinearAlgebra::BlockSparseMatrix &system_matrix;
    };



    void StokesBlock::vmult (LinearAlgebra::BlockVector       &dst,
                             const LinearAlgebra::BlockVector &src) const
    {
      system_matrix.block(0,0).vmult(dst.block(0), src.block(0));
      system_matrix.block(0,1).vmult_add(dst.block(0), src.block(1));

      system_matrix.block(1,0).vmult(dst.block(1), src.block(0));
      system_matrix.block(1,1).vmult_add(dst.block(1), src.block(1));
    }

    void StokesBlock::Tvmult (LinearAlgebra::BlockVector       &dst,
                              const LinearAlgebra::BlockVector &src) const
    {
      system_matrix.block(0,0).Tvmult(dst.block(0), src.block(0));
      system_matrix.block(1,0).Tvmult_add(dst.block(0), src.block(1));

      system_matrix.block(0,1).Tvmult(dst.block(1), src.block(0));
      system_matrix.block(1,1).Tvmult_add(dst.block(1), src.block(1));
    }

    void StokesBlock::vmult_add (LinearAlgebra::BlockVector       &dst,
                                 const LinearAlgebra::BlockVector &src) const
    {
      system_matrix.block(0,0).vmult_add(dst.block(0), src.block(0));
      system_matrix.block(0,1).vmult_add(dst.block(0), src.block(1));

      system_matrix.block(1,0).vmult_add(dst.block(1), src.block(0));
      system_matrix.block(1,1).vmult_add(dst.block(1), src.block(1));
    }

    void StokesBlock::Tvmult_add (LinearAlgebra::BlockVector       &dst,
                                  const LinearAlgebra::BlockVector &src) const
    {
      system_matrix.block(0,0).Tvmult_add(dst.block(0), src.block(0));
      system_matrix.block(1,0).Tvmult_add(dst.block(0), src.block(1));

      system_matrix.block(0,1).Tvmult_add(dst.block(1), src.block(0));
      system_matrix.block(1,1).Tvmult_add(dst.block(1), src.block(1));
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

        TrilinosWrappers::SolverCG solver(solver_control);

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
          TrilinosWrappers::SolverCG solver(solver_control);
          solver.solve(stokes_matrix.block(0,0), dst.block(0), utmp,
                       a_preconditioner);
        }
      else
        a_preconditioner.vmult (dst.block(0), utmp);
    }

  }



  template <int dim>
  void Simulator<dim>::solve_temperature ()
  {
    computing_timer.enter_section ("   Solve temperature system");
    {
      pcout << "   Solving temperature system... " << std::flush;

      SolverControl solver_control (system_matrix.block(2,2).m(),
                                    1e-12*system_rhs.block(2).l2_norm());

      SolverGMRES<LinearAlgebra::Vector>   solver (solver_control,
                                                   SolverGMRES<LinearAlgebra::Vector>::AdditionalData(30,true));

      LinearAlgebra::BlockVector
      distributed_solution (system_rhs);
      distributed_solution = solution;

      solver.solve (system_matrix.block(2,2), distributed_solution.block(2),
                    system_rhs.block(2), *T_preconditioner);

      current_constraints.distribute (distributed_solution);
      solution.block(2) = distributed_solution.block(2);

      // print number of iterations and also record it in the
      // statistics file
      pcout << solver_control.last_step()
            << " iterations." << std::endl;

      statistics.add_value("Iterations for temperature solver",
                           solver_control.last_step());
    }
    computing_timer.exit_section();
  }



  template <int dim>
  void Simulator<dim>::solve_stokes ()
  {
    computing_timer.enter_section ("   Solve Stokes system");

    pcout << "   Solving Stokes system... " << std::flush;

    // extract Stokes parts of solution vector, without any ghost elements
    LinearAlgebra::BlockVector distributed_stokes_solution;
    distributed_stokes_solution.reinit(system_rhs);
    distributed_stokes_solution.block(0) = solution.block(0);
    distributed_stokes_solution.block(1) = solution.block(1);

    // before solving we scale the initial solution to the right dimensions
    distributed_stokes_solution.block(1) /= pressure_scaling;

    // TODO: why do we only set the constraints for the pressure? It should not hurt
    // do that also for the velocity. Note that we should implement constraint_matrix.set_zero()
    // in deal.II and call it instead of this manual code. (TH)
    const unsigned int
    start = (distributed_stokes_solution.block(0).size() +
             distributed_stokes_solution.block(1).local_range().first),
            end   = (distributed_stokes_solution.block(0).size() +
                     distributed_stokes_solution.block(1).local_range().second);
    for (unsigned int i=start; i<end; ++i)
      if (current_constraints.is_constrained (i))
        distributed_stokes_solution(i) = 0;

    // if the model is compressible then we need to adjust the right hand
    // side of the equation to make it compatible with the matrix on the
    // left
    if (material_model->is_compressible ())
      make_pressure_rhs_compatible(system_rhs);

    // extract Stokes parts of rhs vector
    LinearAlgebra::BlockVector distributed_stokes_rhs;
    distributed_stokes_rhs.reinit(system_rhs);
    distributed_stokes_rhs.block(0) = system_rhs.block(0);
    distributed_stokes_rhs.block(1) = system_rhs.block(1);

    PrimitiveVectorMemory< LinearAlgebra::BlockVector > mem;

    const internal::StokesBlock stokes_block(system_matrix);

    // step 1a: try if the simple and fast solver
    // succeeds in 30 steps or less.
    const double solver_tolerance = 1e-7 * distributed_stokes_rhs.l2_norm();
    SolverControl solver_control_cheap (30, solver_tolerance);
    SolverControl solver_control_expensive (system_matrix.block(0,1).m() +
                                            system_matrix.block(1,0).m(), solver_tolerance);

    try
      {
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

    // then copy back the solution from the temporary (non-ghosted) vector
    // into the ghosted one with all solution components
    solution.block(0) = distributed_stokes_solution.block(0);
    solution.block(1) = distributed_stokes_solution.block(1);

    // now rescale the pressure back to real physical units
    solution.block(1) *= pressure_scaling;

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
  }

}





// explicit instantiation of the functions we implement in this file
namespace aspect
{
  template void Simulator<deal_II_dimension>::solve_temperature ();
  template void Simulator<deal_II_dimension>::solve_stokes ();
}
