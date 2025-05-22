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

#ifndef _aspect_simulator_solver_stokes_matrix_free_h
#define _aspect_simulator_solver_stokes_matrix_free_h

#include <aspect/global.h>
#include <aspect/parameters.h>
#include <aspect/simulator/solver/interface.h>
#include <aspect/utilities.h>

#include <deal.II/multigrid/mg_transfer_matrix_free.h>
#include <deal.II/multigrid/mg_transfer_global_coarsening.templates.h>

#include <deal.II/lac/solver_bicgstab.h>

/**
 * Typedef for the number type for the multigrid operators. Can be either float or double.
 */
using GMGNumberType = double;

namespace aspect
{
  /**
   * Base class for the matrix free GMG solver for the Stokes system. The
   * actual implementation is found inside StokesMatrixFreeHandlerLocalSmoothingImplementation.
   */
  template <int dim>
  class StokesMatrixFreeHandler: public StokesSolver::Interface<dim>
  {
    public:
      /**
       * Allocates and sets up the members of the StokesMatrixFreeHandler. This
       * is called by Simulator<dim>::setup_dofs()
       */
      virtual void setup_dofs()=0;

      /**
       * Perform various tasks to update the linear system to solve
       * for. Note that we are not assembling a matrix (as this is a
       * matrix-free algorithm), but we are evaluating the material
       * model and storing the information necessary for a later call
       * to solve().
       */
      virtual void assemble()=0;

      /**
       * Computes and sets the diagonal for both the mass matrix
       * operator and the A-block operators on each level for the
       * purpose of smoothing inside the multigrid v-cycle.
       */
      virtual void build_preconditioner()=0;

      /**
       * Declare parameters necessary for this solver.
       */
      static
      void declare_parameters (ParameterHandler &prm);

      /**
       * Return a reference to the DoFHandler that is used for velocity in
       * the block GMG solver.
       */
      virtual const DoFHandler<dim> &
      get_dof_handler_v () const = 0;

      /**
       * Return a reference to the DoFHandler that is used for pressure in
       * the block GMG solver.
       */
      virtual const DoFHandler<dim> &
      get_dof_handler_p () const = 0;

      /**
       * Return a reference to the DoFHandler that is used for the coefficient
       * projection in the block GMG solver.
       */
      virtual const DoFHandler<dim> &
      get_dof_handler_projection () const = 0;

      /**
       * Return a pointer to the object that describes the velocity DoF
       * constraints for the block GMG Stokes solver.
       */
      virtual const AffineConstraints<double> &
      get_constraints_v () const = 0;

      /**
       * Return a pointer to the object that describes the pressure DoF
       * constraints for the block GMG Stokes solver.
       */
      virtual const AffineConstraints<double> &
      get_constraints_p () const = 0;

      /**
       * Return a pointer to the MGTransfer object used for the A block
       * of the block GMG Stokes solver.
       */
      virtual const MGTransferMF<dim,GMGNumberType> &
      get_mg_transfer_A () const = 0;

      /**
       * Return a pointer to the MGTransfer object used for the Schur
       * complement block of the block GMG Stokes solver.
       */
      virtual const MGTransferMF<dim,GMGNumberType> &
      get_mg_transfer_S () const = 0;

      /**
       * Return the memory consumption in bytes that are used to store
       * equation data like viscosity to be able to apply the operators.
       */
      virtual std::size_t get_cell_data_memory_consumption() const = 0;
  };

  /**
   * Create an instance of the StokesMatrixFreeHandler based on the input parameters.
   */
  template <int dim>
  std::unique_ptr<StokesMatrixFreeHandler<dim>> create_matrix_free_solver(Simulator<dim> &simulator, const Parameters<dim> &parameters);


  namespace internal
  {

    /**
     * Implement the block Schur preconditioner for the Stokes system.
     */
    template <class StokesMatrixType, class ABlockMatrixType, class SchurComplementMatrixType,
              class ABlockPreconditionerType, class SchurComplementPreconditionerType>
    class BlockSchurGMGPreconditioner : public Subscriptor
    {
      public:
        /**
        * @brief Constructor
        *
        * @param Stokes_matrix The entire Stokes matrix
        * @param A_block The A block of the Stokes matrix
        * @param Schur_complement_block The matrix which describes the Schur complement approximation
        * @param A_block_preconditioner Preconditioner object for the matrix A.
        * @param Schur_complement_preconditioner Preconditioner object for the Schur complement.
        * @param do_solve_A A flag indicating whether we should actually solve with
        *     the matrix $A_block$, or only apply one preconditioner step with it.
        * @param do_solve_Schur_complement A flag indicating whether we should actually solve with
        *     the matrix $Schur_complement_block$, or only apply one preconditioner step with it.
        * @param A_block_is_symmetric A flag indicating whether the A block is symmetric.
        * @param A_block_tolerance The tolerance for the CG solver which computes
        *     the inverse of the A block.
        * @param Schur_complement_tolerance The tolerance for the CG solver which computes
        *     the inverse of the Schur complement block (Schur complement approximation matrix).
        */
        BlockSchurGMGPreconditioner (const StokesMatrixType                  &Stokes_matrix,
                                     const ABlockMatrixType                  &A_block,
                                     const SchurComplementMatrixType         &Schur_complement_block,
                                     const ABlockPreconditionerType          &A_block_preconditioner,
                                     const SchurComplementPreconditionerType &Schur_complement_preconditioner,
                                     const bool                               do_solve_A,
                                     const bool                               do_solve_Schur_complement,
                                     const bool                               A_block_is_symmetric,
                                     const double                             A_block_tolerance,
                                     const double                             Schur_complement_tolerance);

        /**
        * Matrix vector product with this preconditioner object.
        */
        void vmult (dealii::LinearAlgebra::distributed::BlockVector<double>       &dst,
                    const dealii::LinearAlgebra::distributed::BlockVector<double> &src) const;

        unsigned int n_iterations_A_block() const;
        unsigned int n_iterations_Schur_complement() const;


      private:
        /**
        * References to the various matrix object this preconditioner works on.
        */
        const StokesMatrixType                  &stokes_matrix;
        const ABlockMatrixType                  &A_block;
        const SchurComplementMatrixType         &Schur_complement_block;
        const ABlockPreconditionerType          &A_block_preconditioner;
        const SchurComplementPreconditionerType &Schur_complement_preconditioner;

        /**
        * Whether to actually invert the $\tilde M$ or $\tilde A$ of the preconditioner matrix
        * or to just apply a single preconditioner step with it.
        */
        const bool                                                      do_solve_A;
        const bool                                                      do_solve_Schur_complement;
        const bool                                                      A_block_is_symmetric;
        mutable unsigned int                                            n_iterations_A_;
        mutable unsigned int                                            n_iterations_Schur_complement_;
        const double                                                    A_block_tolerance;
        const double                                                    Schur_complement_tolerance;
        mutable dealii::LinearAlgebra::distributed::BlockVector<double> utmp;
    };

    template <class StokesMatrixType, class ABlockMatrixType, class SchurComplementMatrixType,
              class ABlockPreconditionerType, class SchurComplementPreconditionerType>
    BlockSchurGMGPreconditioner<StokesMatrixType, ABlockMatrixType, SchurComplementMatrixType,
                                ABlockPreconditionerType, SchurComplementPreconditionerType>::
                                BlockSchurGMGPreconditioner (const StokesMatrixType                  &Stokes_matrix,
                                                             const ABlockMatrixType                  &A_block,
                                                             const SchurComplementMatrixType         &Schur_complement_block,
                                                             const ABlockPreconditionerType          &A_block_preconditioner,
                                                             const SchurComplementPreconditionerType &Schur_complement_preconditioner,
                                                             const bool                               do_solve_A,
                                                             const bool                               do_solve_Schur_complement,
                                                             const bool                               A_block_symmetric,
                                                             const double                             A_block_tolerance,
                                                             const double                             Schur_complement_tolerance)
                                  :
                                  stokes_matrix                   (Stokes_matrix),
                                  A_block                         (A_block),
                                  Schur_complement_block          (Schur_complement_block),
                                  A_block_preconditioner          (A_block_preconditioner),
                                  Schur_complement_preconditioner (Schur_complement_preconditioner),
                                  do_solve_A                      (do_solve_A),
                                  do_solve_Schur_complement       (do_solve_Schur_complement),
                                  A_block_is_symmetric            (A_block_symmetric),
                                  n_iterations_A_                 (0),
                                  n_iterations_Schur_complement_  (0),
                                  A_block_tolerance               (A_block_tolerance),
                                  Schur_complement_tolerance      (Schur_complement_tolerance)
    {}

    template <class StokesMatrixType, class ABlockMatrixType, class SchurComplementMatrixType,
              class ABlockPreconditionerType, class SchurComplementPreconditionerType>
    unsigned int
    BlockSchurGMGPreconditioner<StokesMatrixType, ABlockMatrixType, SchurComplementMatrixType,
                                ABlockPreconditionerType, SchurComplementPreconditionerType>::
                                n_iterations_A_block() const
    {
      return n_iterations_A_;
    }

    template <class StokesMatrixType, class ABlockMatrixType, class SchurComplementMatrixType,
              class ABlockPreconditionerType, class SchurComplementPreconditionerType>
    unsigned int
    BlockSchurGMGPreconditioner<StokesMatrixType, ABlockMatrixType, SchurComplementMatrixType,
                                ABlockPreconditionerType, SchurComplementPreconditionerType>::
                                n_iterations_Schur_complement() const
    {
      return n_iterations_Schur_complement_;
    }

    template <class StokesMatrixType, class ABlockMatrixType, class SchurComplementMatrixType,
              class ABlockPreconditionerType, class SchurComplementPreconditionerType>
    void
    BlockSchurGMGPreconditioner<StokesMatrixType, ABlockMatrixType, SchurComplementMatrixType,
                                ABlockPreconditionerType, SchurComplementPreconditionerType>::
                                vmult (dealii::LinearAlgebra::distributed::BlockVector<double>       &dst,
                                       const dealii::LinearAlgebra::distributed::BlockVector<double>  &src) const
    {
      if (utmp.size()==0)
        utmp.reinit(src);

      // This needs to be done explicitly, as GMRES does not
      // initialize the data of the vector dst before calling
      // us. Otherwise we might use random data as our initial guess.
      dst = 0.0;

      // either solve with the Schur complement matrix (if do_solve_Schur_complement==true)
      // or just apply one preconditioner sweep (for the first few
      // iterations of our two-stage outer GMRES iteration)
      if (do_solve_Schur_complement)
        {
          // first solve with the bottom right block, which we have built
          // as a mass matrix with the inverse of the viscosity
          SolverControl solver_control(100, src.block(1).l2_norm() * Schur_complement_tolerance,true);

          SolverCG<dealii::LinearAlgebra::distributed::Vector<double>> solver(solver_control);
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
                  solver.solve(Schur_complement_block,
                               dst.block(1), src.block(1),
                               Schur_complement_preconditioner);
                  n_iterations_Schur_complement_ += solver_control.last_step();
                }
              // if the solver fails, report the error from processor 0 with some additional
              // information about its location, and throw a quiet exception on all other
              // processors
              catch (const std::exception &exc)
                {
                  Utilities::throw_linear_solver_failure_exception("iterative (bottom right) solver",
                                                                   "BlockSchurGMGPreconditioner::vmult",
                                                                   std::vector<SolverControl> {solver_control},
                                                                   exc,
                                                                   src.block(0).get_mpi_communicator());
                }
            }
        }
      else
        {
          Schur_complement_preconditioner.vmult(dst.block(1),src.block(1));
          n_iterations_Schur_complement_ += 1;
        }

      dst.block(1) *= -1.0;

      {
        // the u-block of dst only contains zeros
        stokes_matrix.vmult(utmp, dst); // B^T
        utmp.block(0) *= -1.0;
        utmp.block(0) += src.block(0);
      }

      // now either solve with the top left block (if do_solve_A==true)
      // or just apply one preconditioner sweep (for the first few
      // iterations of our two-stage outer GMRES iteration)
      if (do_solve_A == true)
        {
          SolverControl solver_control(1000, utmp.block(0).l2_norm() * A_block_tolerance);
          PrimitiveVectorMemory<dealii::LinearAlgebra::distributed::Vector<double>> mem;

          try
            {
              if (A_block_is_symmetric)
                {
                  SolverCG<dealii::LinearAlgebra::distributed::Vector<double>> solver(solver_control,mem);
                  solver.solve(A_block, dst.block(0), utmp.block(0),
                               A_block_preconditioner);
                }
              else
                {
                  // Use BiCGStab for non-symmetric matrices.
                  // BiCGStab can also solve indefinite systems if necessary.
                  // Do not compute the exact residual, as this
                  // is more expensive, and we only need an approximate solution.
                  SolverBicgstab<dealii::LinearAlgebra::distributed::Vector<double>>
                  solver(solver_control,
                         mem,
                         SolverBicgstab<dealii::LinearAlgebra::distributed::Vector<double>>::AdditionalData(/*exact_residual=*/ false));
                  solver.solve(A_block, dst.block(0), utmp.block(0),
                               A_block_preconditioner);
                }

              n_iterations_A_ += solver_control.last_step();
            }
          // if the solver fails, report the error from processor 0 with some additional
          // information about its location, and throw a quiet exception on all other
          // processors
          catch (const std::exception &exc)
            {
              Utilities::throw_linear_solver_failure_exception("iterative (top left) solver",
                                                               "BlockSchurGMGPreconditioner::vmult",
                                                               std::vector<SolverControl> {solver_control},
                                                               exc,
                                                               src.block(0).get_mpi_communicator());
            }
        }
      else
        {
          A_block_preconditioner.vmult (dst.block(0), utmp.block(0));
          n_iterations_A_ += 1;
        }
    }
  }
}

#endif
