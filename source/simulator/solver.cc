/* $Id$ */
/* Author: Martin Kronbichler, Uppsala University,
           Wolfgang Bangerth, Texas A&M University,
     Timo Heister, University of Goettingen, 2008-2011 */
/*                                                                */
/*    Copyright (C) 2008, 2009, 2010, 2011, 2012 by the deal.II authors */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */

#include <aspect/simulator.h>
#include <aspect/global.h>

#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/trilinos_solver.h>


namespace aspect
{
  namespace internal
  {
    using namespace dealii;

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
        BlockSchurPreconditioner (const TrilinosWrappers::BlockSparseMatrix  &S,
                                  const TrilinosWrappers::BlockSparseMatrix  &Spre,
                                  const PreconditionerMp                     &Mppreconditioner,
                                  const PreconditionerA                      &Apreconditioner,
                                  const bool                                  do_solve_A);

        /**
         * Matrix vector product with this preconditioner object.
         */
        void vmult (TrilinosWrappers::MPI::BlockVector       &dst,
                    const TrilinosWrappers::MPI::BlockVector &src) const;

      private:
        /**
         * References to the various matrix object this preconditioner works on.
         */
        const TrilinosWrappers::BlockSparseMatrix &stokes_matrix;
        const TrilinosWrappers::BlockSparseMatrix &stokes_preconditioner_matrix;
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
    BlockSchurPreconditioner (const TrilinosWrappers::BlockSparseMatrix  &S,
                              const TrilinosWrappers::BlockSparseMatrix  &Spre,
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
    vmult (TrilinosWrappers::MPI::BlockVector       &dst,
           const TrilinosWrappers::MPI::BlockVector &src) const
    {
      TrilinosWrappers::MPI::Vector utmp(src.block(0));

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
  void Simulator<dim>::solve ()
  {
    computing_timer.enter_section ("   Solve Stokes system");

    // STEP 1: solve the Stokes system
    {
      pcout << "   Solving Stokes system... " << std::flush;

      TrilinosWrappers::MPI::BlockVector
      distributed_stokes_solution (stokes_rhs);
      distributed_stokes_solution = stokes_solution;

      // before solving we scale the initial solution to the right dimensions
      distributed_stokes_solution.block(1) /= pressure_scaling;

      const unsigned int
      start = (distributed_stokes_solution.block(0).size() +
               distributed_stokes_solution.block(1).local_range().first),
              end   = (distributed_stokes_solution.block(0).size() +
                       distributed_stokes_solution.block(1).local_range().second);
      for (unsigned int i=start; i<end; ++i)
        if (stokes_constraints.is_constrained (i))
          distributed_stokes_solution(i) = 0;

      // if the model is compressible then we need to adjust the right hand
      // side of the equation to make it compatible with the matrix on the
      // left
      if (material_model->is_compressible ())
        make_pressure_rhs_compatible(stokes_rhs);

      PrimitiveVectorMemory< TrilinosWrappers::MPI::BlockVector > mem;

      // step 1a: try if the simple and fast solver
      // succeeds in 30 steps or less.
      const double solver_tolerance = 1e-7 * stokes_rhs.l2_norm();
      SolverControl solver_control_cheap (30, solver_tolerance);
      SolverControl solver_control_expensive (stokes_matrix.m(), solver_tolerance);

      try
        {
          const internal::BlockSchurPreconditioner<TrilinosWrappers::PreconditionAMG,
                TrilinosWrappers::PreconditionILU>
                preconditioner (stokes_matrix, stokes_preconditioner_matrix,
                                *Mp_preconditioner, *Amg_preconditioner,
                                false);

          SolverFGMRES<TrilinosWrappers::MPI::BlockVector>
          solver(solver_control_cheap, mem,
                 SolverFGMRES<TrilinosWrappers::MPI::BlockVector>::
                 AdditionalData(30, true));
          solver.solve(stokes_matrix, distributed_stokes_solution, stokes_rhs,
                       preconditioner);
        }

      // step 1b: take the stronger solver in case
      // the simple solver failed
      catch (SolverControl::NoConvergence)
        {
          const internal::BlockSchurPreconditioner<TrilinosWrappers::PreconditionAMG,
                TrilinosWrappers::PreconditionILU>
                preconditioner (stokes_matrix, stokes_preconditioner_matrix,
                                *Mp_preconditioner, *Amg_preconditioner,
                                true);

          SolverFGMRES<TrilinosWrappers::MPI::BlockVector>
          solver(solver_control_expensive, mem,
                 SolverFGMRES<TrilinosWrappers::MPI::BlockVector>::
                 AdditionalData(50, true));
          solver.solve(stokes_matrix, distributed_stokes_solution, stokes_rhs,
                       preconditioner);
        }


      stokes_constraints.distribute (distributed_stokes_solution);

      // now rescale the pressure back to real physical units
      distributed_stokes_solution.block(1) *= pressure_scaling;

      stokes_solution = distributed_stokes_solution;

      normalize_pressure(stokes_solution);

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
    }
    computing_timer.exit_section();


    // STEP 2: Update the time step size
    {
      old_time_step = time_step;
      time_step = compute_time_step();

      if (parameters.convert_to_years == true)
        statistics.add_value("Time step size (years)", time_step / year_in_seconds);
      else
        statistics.add_value("Time step size (seconds)", time_step);

      temperature_solution = old_temperature_solution;
    }

    // STEP 3: Set up the temperature linear system
    assemble_temperature_system ();

    // STEP 4: Solve the temperature system
    computing_timer.enter_section ("   Solve temperature system");
    {
      pcout << "   Solving temperature system... " << std::flush;

      SolverControl solver_control (temperature_matrix.m(),
                                    1e-12*temperature_rhs.l2_norm());
      SolverGMRES<TrilinosWrappers::MPI::Vector>   solver (solver_control,
                                                           SolverGMRES<TrilinosWrappers::MPI::Vector>::AdditionalData(30,true));

      TrilinosWrappers::MPI::Vector
      distributed_temperature_solution (temperature_rhs);
      distributed_temperature_solution = temperature_solution;

      solver.solve (temperature_matrix, distributed_temperature_solution,
                    temperature_rhs, *T_preconditioner);

      temperature_constraints.distribute (distributed_temperature_solution);
      temperature_solution = distributed_temperature_solution;

      // print number of iterations and also record it in the
      // statistics file
      pcout << solver_control.last_step()
            << " iterations." << std::endl;

      statistics.add_value("Iterations for temperature solver",
                           solver_control.last_step());
    }
    computing_timer.exit_section();
  }
}


// explicit instantiation of the functions we implement in this file
namespace aspect
{
  template void Simulator<deal_II_dimension>::solve ();
}
