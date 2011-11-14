//-------------------------------------------------------------
//    $Id: simulator.h 232 2011-10-19 13:30:15Z bangerth $
//
//    Copyright (C) 2011 by the authors of the ASPECT code
//
//-------------------------------------------------------------
#ifndef __aspect__solver_h
#define __aspect__solver_h

#include <deal.II/lac/trilinos_block_vector.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>

namespace LinearSolvers
{
  using namespace dealii;

  template <class PreconditionerA, class PreconditionerMp>
  class BlockSchurPreconditioner : public Subscriptor
  {
    public:
      BlockSchurPreconditioner (const TrilinosWrappers::BlockSparseMatrix  &S,
                                const TrilinosWrappers::BlockSparseMatrix  &Spre,
                                const PreconditionerMp                     &Mppreconditioner,
                                const PreconditionerA                      &Apreconditioner,
                                const bool                                  do_solve_A)
        :
        stokes_matrix     (&S),
        stokes_preconditioner_matrix     (&Spre),
        mp_preconditioner (Mppreconditioner),
        a_preconditioner  (Apreconditioner),
        do_solve_A        (do_solve_A)
      {}

      void vmult (TrilinosWrappers::MPI::BlockVector       &dst,
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
            solver.solve(stokes_preconditioner_matrix->block(1,1),
                         dst.block(1), src.block(1),
                         mp_preconditioner);

          dst.block(1) *= -1.0;
        }

        {
          stokes_matrix->block(0,1).vmult(utmp, dst.block(1)); //B^T
          utmp*=-1.0;
          utmp.add(src.block(0));
        }

        if (do_solve_A == true)
          {
            SolverControl solver_control(5000, utmp.l2_norm()*1e-2);
            TrilinosWrappers::SolverCG solver(solver_control);
            solver.solve(stokes_matrix->block(0,0), dst.block(0), utmp,
                         a_preconditioner);
          }
        else
          a_preconditioner.vmult (dst.block(0), utmp);
      }

    private:
      const SmartPointer<const TrilinosWrappers::BlockSparseMatrix> stokes_matrix;
      const SmartPointer<const TrilinosWrappers::BlockSparseMatrix> stokes_preconditioner_matrix;
      const PreconditionerMp &mp_preconditioner;
      const PreconditionerA  &a_preconditioner;
      const bool do_solve_A;
  };
}


#endif
