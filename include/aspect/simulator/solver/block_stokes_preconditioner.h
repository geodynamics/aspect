/*
  Copyright (C) 2011 - 2025 by the authors of the ASPECT code.

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


#ifndef aspect_block_stokes_preconditioner_h
#define aspect_block_stokes_preconditioner_h

#include <deal.II/lac/solver_bicgstab.h>
#include <deal.II/lac/solver_cg.h>


#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/la_parallel_block_vector.h>

namespace aspect
{

  namespace internal
  {
    /**
      * This class is used in the implementation of the right preconditioner
      * as an approximation for the inverse of the velocity (A) block.
      * This operator can either just apply the preconditioner (AMG)
      * or perform an inner CG / BiCGStab solve with the same preconditioner.
      */
    template <class PreconditionerA, class VectorType, class ABlockType>
    class InverseVelocityBlock
    {
      public:
        /**
         * Constructor.
         * @param matrix The matrix that contains A (from the system matrix)
         * @param preconditioner The preconditioner to be used
         * @param do_solve_A A flag indicating whether we should actually solve with
         *     the matrix $A$, or only apply one preconditioner step with it.
         * @param A_block_is_symmetric A flag indicating whether the matrix $A$ is symmetric.
         * @param solver_tolerance The tolerance for the CG solver which computes
         *     the inverse of the A block.
         */
        InverseVelocityBlock(const ABlockType &matrix,
                             const PreconditionerA &preconditioner,
                             const bool do_solve_A,
                             const bool A_block_is_symmetric,
                             const double solver_tolerance);

        void vmult(VectorType &dst,
                   const VectorType &src) const;

        unsigned int n_iterations() const;

      private:
        mutable unsigned int n_iterations_;
        const ABlockType &matrix;
        const PreconditionerA &preconditioner;
        const bool do_solve_A;
        const bool A_block_is_symmetric;
        const double solver_tolerance;
    };



    template <class PreconditionerA,class VectorType, class ABlockType>
    InverseVelocityBlock<PreconditionerA,VectorType,ABlockType>::InverseVelocityBlock(
      const ABlockType &matrix,
      const PreconditionerA &preconditioner,
      const bool do_solve_A,
      const bool A_block_is_symmetric,
      const double solver_tolerance)
      : n_iterations_ (0),
        matrix (matrix),
        preconditioner (preconditioner),
        do_solve_A (do_solve_A),
        A_block_is_symmetric (A_block_is_symmetric),
        solver_tolerance (solver_tolerance)
    {}



    /**
    * Implements the vmult for InverseVelocityBlock. This applies the action of A^{-1} by either
    * performing a solve with A or using a preconditioner sweep.
    */
    template <class PreconditionerA, class VectorType, class ABlockType>
    void InverseVelocityBlock<PreconditionerA,VectorType,ABlockType>::vmult(VectorType &dst,
                                                                            const VectorType &src) const
    {

      // Either solve with the top left block
      // or just apply one preconditioner sweep (for the first few
      // iterations of our two-stage outer GMRES iteration)
      if (do_solve_A == true)
        {
          SolverControl solver_control(10000, src.l2_norm() * solver_tolerance);
          PrimitiveVectorMemory<VectorType> mem;

          try
            {
              dst = 0.0;

              if (A_block_is_symmetric)
                {
                  SolverCG<VectorType> solver(solver_control, mem);
                  solver.solve(matrix, dst, src, preconditioner);
                }
              else
                {
                  // Use BiCGStab for non-symmetric matrices.
                  // BiCGStab can also solve indefinite systems if necessary.
                  // Do not compute the exact residual, as this
                  // is more expensive, and we only need an approximate solution.
                  SolverBicgstab<VectorType>
                  solver(solver_control,
                         mem,
                         typename SolverBicgstab<VectorType>::AdditionalData(/*exact_residual=*/ false));
                  solver.solve(matrix, dst, src, preconditioner);
                }
              n_iterations_ += solver_control.last_step();
            }
          catch (const std::exception &exc)
            {
              // if the solver fails, report the error from processor 0 with some additional
              // information about its location, and throw a quiet exception on all other
              // processors
              Utilities::throw_linear_solver_failure_exception("iterative (top left) solver",
                                                               "BlockSchurPreconditioner::vmult",
                                                               std::vector<SolverControl> {solver_control},
                                                               exc,
                                                               src.get_mpi_communicator());
            }
        }
      else
        {
          preconditioner.vmult (dst, src);
          n_iterations_ += 1;
        }
    }



    template <class PreconditionerA, class VectorType, class ABlockType>
    unsigned int InverseVelocityBlock<PreconditionerA, VectorType, ABlockType>::n_iterations() const
    {
      return n_iterations_;
    }

    /**
     * Implement the block Schur preconditioner
     * (A B^T; 0 S)^{-1}.
     */
    template <class AInvOperator, class SInvOperator, class BTOperator,  class VectorType>
    class BlockSchurPreconditioner : public
#if DEAL_II_VERSION_GTE(9,7,0)
      EnableObserverPointer
#else
      Subscriptor
#endif

    {
      public:
        /**
         * @brief Constructor
         * @param A_inverse_operator Approximation of the inverse of the velocity block.
         * @param S_inverse_operator Approximation for the inverse Schur complement.
         * @param BT_operator Operator for the B^T block of the Stokes system.
         */
        BlockSchurPreconditioner (
          const AInvOperator                         &A_inverse_operator,
          const SInvOperator                         &S_inverse_operator,
          const BTOperator                           &BT_operator);

        /**
         * Matrix vector product with this preconditioner object.
         */
        void vmult (VectorType       &dst,
                    const VectorType &src) const;

      private:
        /**
         * References to the various operators this preconditioner works with.
         */

        const AInvOperator                     &A_inverse_operator;
        const SInvOperator                     &S_inverse_operator;
        const BTOperator                       &BT_operator;
        mutable VectorType                      tmp;
    };


    template <class AInvOperator, class SInvOperator, class BTOperator,  class VectorType>
    BlockSchurPreconditioner<AInvOperator, SInvOperator, BTOperator, VectorType>::
    BlockSchurPreconditioner (
      const AInvOperator                         &A_inverse_operator,
      const SInvOperator                         &S_inverse_operator,
      const BTOperator                           &BT_operator)
      :
      A_inverse_operator (A_inverse_operator),
      S_inverse_operator (S_inverse_operator),
      BT_operator        (BT_operator)
    {}



    template <class AInvOperator, class SInvOperator, class BTOperator, class VectorType>
    void
    BlockSchurPreconditioner<AInvOperator, SInvOperator, BTOperator, VectorType>::
    vmult (VectorType       &dst,
           const VectorType &src) const
    {
      if (tmp.size() == 0)
        tmp.reinit(src);


      dst=0.0;

      // first apply the Schur Complement inverse operator.
      {
        S_inverse_operator.vmult(dst.block(1),src.block(1));
        dst.block(1) *= -1.0;
      }


      // apply the top right block: B^T or J^{up}
      {
        // For matrix-free usage, the BT operator needs to operate on the whole vector,
        // but the matrix-based version is a SparseMatrix::vmult() that requires passing the
        // individual blocks.
        if constexpr (std::is_same_v<VectorType, dealii::LinearAlgebra::distributed::BlockVector<double>>)
          BT_operator.vmult(tmp, dst);

        else
          BT_operator.vmult(tmp.block(0), dst.block(1));

        tmp.block(0) *= -1.0;
        tmp.block(0) += src.block(0);
      }

      A_inverse_operator.vmult(dst.block(0), tmp.block(0));
    }



    template<class OperatorType, class StokesMatrixType, class SchurComplementMatrixType, class VectorType>
    class SchurApproximation
    {
      public:
        /**
         * @brief Constructor
         * @param schur_preconditioner Preconditioner for the Schur Complement.
         * @param stokes_matrix Stokes system. This is not necessarily a matrix.
         * @param Schur_complement_block Operator for the Schur complement.
         * @param do_solve_Schur_complement Flag that determines whether to do a full solve with the Schur complement or a v-cycle.
         * @param Schur_complement_tolerance Tolerance in case a full solve for the Schur complement is used.
         */
        SchurApproximation(const OperatorType &schur_preconditioner,
                           const StokesMatrixType &stokes_matrix,
                           const SchurComplementMatrixType &Schur_complement_block,
                           const bool do_solve_Schur_complement,
                           const double Schur_complement_tolerance);
        void vmult( VectorType &dst, const VectorType &src) const;
        unsigned int n_iterations() const;

      private:
        const OperatorType &schur_preconditioner;
        const StokesMatrixType &stokes_matrix;
        const SchurComplementMatrixType &Schur_complement_block;
        const bool do_solve_Schur_complement;
        const double Schur_complement_tolerance;
        mutable unsigned int n_iterations_Schur_complement_;


    };



    template <class OperatorType, class StokesMatrixType, class SchurComplementMatrixType, class VectorType>
    SchurApproximation<OperatorType, StokesMatrixType, SchurComplementMatrixType, VectorType>::SchurApproximation(const OperatorType &schur_preconditioner,
        const StokesMatrixType &stokes_matrix,
        const SchurComplementMatrixType &Schur_complement_block,
        const bool do_solve_Schur_complement,
        const double Schur_complement_tolerance)
      :
      schur_preconditioner(schur_preconditioner),
      stokes_matrix(stokes_matrix),
      Schur_complement_block(Schur_complement_block),
      do_solve_Schur_complement(do_solve_Schur_complement),
      Schur_complement_tolerance(Schur_complement_tolerance),
      n_iterations_Schur_complement_(0)
    {}


    template<class OperatorType, class StokesMatrixType, class SchurComplementMatrixType, class VectorType>
    void SchurApproximation<OperatorType, StokesMatrixType, SchurComplementMatrixType, VectorType>::vmult(VectorType &dst,
        const VectorType &src) const
    {
      if (do_solve_Schur_complement)
        {

          SolverControl solver_control(100, src.l2_norm() * Schur_complement_tolerance,true);

          SolverCG<VectorType> solver(solver_control);
          // Trilinos reports a breakdown in case src=dst=0, even
          // though it should return convergence without
          // iterating. We simply skip solving in this case.
          if (src.l2_norm() > 1e-50)
            {
              try
                {
                  // explicitly zero out because GMRES does not guarantee that dst is zeroed out
                  dst = 0.0;

                  solver.solve(Schur_complement_block,
                               dst, src,
                               schur_preconditioner);
                  n_iterations_Schur_complement_ += solver_control.last_step();
                }
              // if the solver fails, report the error from processor 0 with some additional
              // information about its location, and throw a quiet exception on all other
              // processors
              catch (const std::exception &exc)
                {
                  Utilities::throw_linear_solver_failure_exception("iterative (bottom right) solver",
                                                                   "BlockSchurPreconditioner::vmult",
                                                                   std::vector<SolverControl> {solver_control},
                                                                   exc,
                                                                   src.get_mpi_communicator());
                }
            }
        }
      else
        {
          dst = 0.0;
          schur_preconditioner.vmult(dst,src);
          n_iterations_Schur_complement_ += 1;
        }
    }
    template<class OperatorType, class StokesMatrixType, class SchurComplementMatrixType, class VectorType>
    unsigned int SchurApproximation<OperatorType, StokesMatrixType, SchurComplementMatrixType, VectorType>::n_iterations() const
    {
      return n_iterations_Schur_complement_;
    }


  }
}



#endif
