/*
  Copyright (C) 2018 - 2024 by the authors of the ASPECT code.

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

#ifndef _aspect_simulator_stokes_matrix_free_operators_h
#define _aspect_simulator_stokes_matrix_free_operators_h

#include <aspect/global.h>
#include <aspect/simulator/solver/interface.h>
#include <aspect/simulator/solver/stokes_matrix_free.h>
#include <aspect/simulator.h>

#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>
#include <deal.II/matrix_free/fe_evaluation.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/multigrid/multigrid.h>
#include <deal.II/multigrid/mg_transfer_matrix_free.h>
#include <deal.II/multigrid/mg_transfer_global_coarsening.templates.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_matrix.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/la_parallel_block_vector.h>

namespace aspect
{
  namespace internal
  {
    /**
     * Matrix-free operators must use deal.II defined vectors, while the rest of the ASPECT
     * software is based on Trilinos vectors. Here we define functions which copy between the
     * vector types.
     */
    namespace ChangeVectorTypes
    {
      void import(TrilinosWrappers::MPI::Vector &out,
                  const dealii::LinearAlgebra::ReadWriteVector<double> &rwv,
                  const VectorOperation::values                 operation);

      void copy(TrilinosWrappers::MPI::Vector &out,
                const dealii::LinearAlgebra::distributed::Vector<double> &in);

      void copy(dealii::LinearAlgebra::distributed::Vector<double> &out,
                const TrilinosWrappers::MPI::Vector &in);

      void copy(TrilinosWrappers::MPI::BlockVector &out,
                const dealii::LinearAlgebra::distributed::BlockVector<double> &in);

      void copy(dealii::LinearAlgebra::distributed::BlockVector<double> &out,
                const TrilinosWrappers::MPI::BlockVector &in);
    }
  }

  /**
   * This namespace contains all matrix-free operators used in the Stokes solver.
   */
  namespace MatrixFreeStokesOperators
  {

    /**
     * This struct stores the data for the current linear operator that is required to perform
     * matrix-vector products.
     *
     * The members of type Table<2, VectorizedArray<X>> contain values
     * of type X, grouped by cell batch using the VectorizedArray. The
     * table is indexed by the index of the cell batch and quadrature
     * point index.  In other words, you can access the value by
     * <tt>table(cell_batch_index, q_index)[cell_index]</tt>
     */
    template <int dim, typename number>
    struct OperatorCellData
    {
      /**
       * Information on the compressibility of the flow.
       */
      bool is_compressible;

      /**
       * Pressure scaling constant.
       */
      double pressure_scaling;

      /**
       * If true, Newton terms are part of the operator.
       */
      bool enable_newton_derivatives;

      /**
       * Symmetrize the Newton system when it's true (i.e., the
       * stabilization is symmetric or SPD).
       */
      bool symmetrize_newton_system;

      /**
       * If true, apply the stabilization on free surface faces.
       */
      bool apply_stabilization_free_surface_faces;

      /**
       * Table which stores viscosity values for each cell.
       *
       * If the second dimension is of size 1, the viscosity is
       * assumed to be constant per cell.
       */
      Table<2, VectorizedArray<number>> viscosity;

      /**
       * Table which stores the strain rate for each cell to be used
       * for the Newton terms.
       */
      Table<2, SymmetricTensor<2, dim, VectorizedArray<number>>> strain_rate_table;

      /**
       * Table which stores the product of the following three
       * variables: viscosity derivative with respect to pressure,
       * the Newton derivative scaling factor, and the averaging weight.
       */
      Table<2, VectorizedArray<number>> newton_factor_wrt_pressure_table;

      /**
       * Table which stores the product of the following four
       * variables: viscosity derivative with respect to strain rate,
       * newton derivative scaling factor, alpha, and the averaging
       * weight. Here alpha is the spd factor when the stabilization
       * is PD or SPD, otherwise, it is 1.
       */
      Table<2, SymmetricTensor<2, dim, VectorizedArray<number>>>
      newton_factor_wrt_strain_rate_table;

      /**
       * Table which stores the product of the pressure perturbation
       * and the normalized gravity. The size is n_face_boundary * n_face_q_points,
       * but only those on the free surface are computed and stored.
       */
      Table<2, Tensor<1, dim, VectorizedArray<number>>> free_surface_stabilization_term_table;

      /**
       * Boundary indicators of those boundaries with a free surface.
       */
      std::set<types::boundary_id> free_surface_boundary_indicators;

      /**
       * Determine an estimate for the memory consumption (in bytes) of this
       * object.
       */
      std::size_t
      memory_consumption() const;

      /**
       * Reset the object and free all memory
       */
      void clear();
    };

    /**
     * Operator for the entire Stokes block.
     */
    template <int dim, int degree_v, typename number>
    class StokesOperator
      : public MatrixFreeOperators::Base<dim, dealii::LinearAlgebra::distributed::BlockVector<number>>
    {
      public:

        /**
         * Constructor.
         */
        StokesOperator ();

        /**
         * Reset object.
         */
        void clear () override;

        /**
         * Pass in a reference to the problem data.
         */
        void set_cell_data (const OperatorCellData<dim,number> &data);

        /**
         * Computes the diagonal of the matrix. Since matrix-free operators have not access
         * to matrix elements, we must apply the matrix-free operator to the unit vectors to
         * recover the diagonal.
         */
        void compute_diagonal () override;

      private:

        /**
         * Performs the application of the matrix-free operator. This function is called by
         * vmult() functions MatrixFreeOperators::Base.
         */
        void apply_add (dealii::LinearAlgebra::distributed::BlockVector<number> &dst,
                        const dealii::LinearAlgebra::distributed::BlockVector<number> &src) const override;

        /**
         * Defines the application of the cell matrix.
         */
        void local_apply (const dealii::MatrixFree<dim, number> &data,
                          dealii::LinearAlgebra::distributed::BlockVector<number> &dst,
                          const dealii::LinearAlgebra::distributed::BlockVector<number> &src,
                          const std::pair<unsigned int, unsigned int> &cell_range) const;

        /**
         * This function doesn't do anything, it's created to use the matrixfree loop.
         */
        void local_apply_face (const dealii::MatrixFree<dim, number> &data,
                               dealii::LinearAlgebra::distributed::BlockVector<number> &dst,
                               const dealii::LinearAlgebra::distributed::BlockVector<number> &src,
                               const std::pair<unsigned int, unsigned int> &face_range) const;

        /**
         * Apply the stabilization on free surface faces.
         */
        void local_apply_boundary_face (const dealii::MatrixFree<dim, number> &data,
                                        dealii::LinearAlgebra::distributed::BlockVector<number> &dst,
                                        const dealii::LinearAlgebra::distributed::BlockVector<number> &src,
                                        const std::pair<unsigned int, unsigned int> &face_range) const;

        /**
         * A pointer to the current cell data that contains viscosity and other required parameters per cell.
         */
        const OperatorCellData<dim,number> *cell_data;
    };



    /**
    * Operator for the B^T block.
    */
    template <int dim, int degree_v, typename number>
    class BTBlockOperator
      : public MatrixFreeOperators::Base<dim, dealii::LinearAlgebra::distributed::BlockVector<number>>
    {
      public:

        /**
         * Constructor.
         */
        BTBlockOperator ();

        /**
         * Reset object.
         */
        void clear () override;

        /**
         * Pass in a reference to the problem data.
         */
        void set_cell_data (const OperatorCellData<dim,number> &data);

        /**
         * Computes the diagonal of the matrix. Since matrix-free operators have not access
         * to matrix elements, we must apply the matrix-free operator to the unit vectors to
         * recover the diagonal.
         */
        void compute_diagonal () override;

      private:

        /**
         * Performs the application of the matrix-free operator. This function is called by
         * vmult() functions MatrixFreeOperators::Base.
         */
        void apply_add (dealii::LinearAlgebra::distributed::BlockVector<number> &dst,
                        const dealii::LinearAlgebra::distributed::BlockVector<number> &src) const override;

        /**
         * Defines the application of the cell matrix.
         */
        void local_apply (const dealii::MatrixFree<dim, number> &data,
                          dealii::LinearAlgebra::distributed::BlockVector<number> &dst,
                          const dealii::LinearAlgebra::distributed::BlockVector<number> &src,
                          const std::pair<unsigned int, unsigned int> &cell_range) const;

        /**
         * This function doesn't do anything, it's created to use the matrixfree loop.
         */
        void local_apply_face (const dealii::MatrixFree<dim, number> &data,
                               dealii::LinearAlgebra::distributed::BlockVector<number> &dst,
                               const dealii::LinearAlgebra::distributed::BlockVector<number> &src,
                               const std::pair<unsigned int, unsigned int> &face_range) const;

        /**
         * A pointer to the current cell data that contains viscosity and other required parameters per cell.
         */
        const OperatorCellData<dim,number> *cell_data;
    };



    /**
     * Operator for the pressure mass matrix used in the block preconditioner
     */
    template <int dim, int degree_p, typename number>
    class MassMatrixOperator
      : public MatrixFreeOperators::Base<dim, dealii::LinearAlgebra::distributed::Vector<number>>
    {
      public:

        /**
         * Constructor
         */
        MassMatrixOperator ();

        /**
         * Reset the object.
         */
        void clear () override;

        /**
         * Pass in a reference to the problem data.
         */
        void set_cell_data (const OperatorCellData<dim,number> &data);

        /**
         * Computes the diagonal of the matrix. Since matrix-free operators have not access
         * to matrix elements, we must apply the matrix-free operator to the unit vectors to
         * recover the diagonal.
         */
        void compute_diagonal () override;

      private:

        /**
         * Performs the application of the matrix-free operator. This function is called by
         * vmult() functions MatrixFreeOperators::Base.
         */
        void apply_add (dealii::LinearAlgebra::distributed::Vector<number> &dst,
                        const dealii::LinearAlgebra::distributed::Vector<number> &src) const override;

        /**
         * Defines the application of the cell matrix.
         */
        void local_apply (const dealii::MatrixFree<dim, number> &data,
                          dealii::LinearAlgebra::distributed::Vector<number> &dst,
                          const dealii::LinearAlgebra::distributed::Vector<number> &src,
                          const std::pair<unsigned int, unsigned int> &cell_range) const;

        /**
         * This function contains the inner-most operation done on a single cell
         */
        void inner_cell_operation(FEEvaluation<dim,
                                  degree_p,
                                  degree_p+2,
                                  1,
                                  number> &pressure) const;

        /**
         * A pointer to the current cell data that contains viscosity and other required parameters per cell.
         */
        const OperatorCellData<dim,number> *cell_data;
    };

    /**
     * Operator for the A block of the Stokes matrix. The same class is used for both
     * active and level mesh operators.
     */
    template <int dim, int degree_v, typename number>
    class ABlockOperator
      : public MatrixFreeOperators::Base<dim, dealii::LinearAlgebra::distributed::Vector<number>>
    {
      public:

        /**
         * Constructor
         */
        ABlockOperator ();

        /**
         * Reset the operator.
         */
        void clear () override;

        /**
         * Pass in a reference to the problem data.
         */
        void set_cell_data (const OperatorCellData<dim,number> &data);

        /**
         * Computes the diagonal of the matrix. Since matrix-free operators have not access
         * to matrix elements, we must apply the matrix-free operator to the unit vectors to
         * recover the diagonal.
         */
        void compute_diagonal () override;

        /**
         * Manually set the diagonal inside the matrix-free object. This function is needed
         * when using tangential constraints as the function compute_diagonal() cannot handle
         * non-Dirichlet boundary conditions.
         */
        void set_diagonal (const dealii::LinearAlgebra::distributed::Vector<number> &diag);

      private:
        /**
         * Defines the inner-most operator on a single cell batch with
         * the loop over quadrature points.
         */
        void inner_cell_operation(FEEvaluation<dim,
                                  degree_v,
                                  degree_v+1,
                                  dim,
                                  number> &velocity) const;

        /**
         * Defines the operation on a single cell batch including
         * load/store and calls inner_cell_operation().
         */
        void cell_operation(FEEvaluation<dim,
                            degree_v,
                            degree_v+1,
                            dim,
                            number> &velocity) const;

        /**
         * Performs the application of the matrix-free operator. This
         * function is called by vmult() functions
         * MatrixFreeOperators::Base.
         */
        void apply_add (dealii::LinearAlgebra::distributed::Vector<number> &dst,
                        const dealii::LinearAlgebra::distributed::Vector<number> &src) const override;

        /**
         * Defines the application of the cell matrix.
         */
        void local_apply (const dealii::MatrixFree<dim, number> &data,
                          dealii::LinearAlgebra::distributed::Vector<number> &dst,
                          const dealii::LinearAlgebra::distributed::Vector<number> &src,
                          const std::pair<unsigned int, unsigned int> &cell_range) const;

        /**
         * A pointer to the current cell data that contains viscosity and other required parameters per cell.
         */
        const OperatorCellData<dim,number> *cell_data;
    };
  }

}

#endif
