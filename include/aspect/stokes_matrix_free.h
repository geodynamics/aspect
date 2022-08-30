/*
  Copyright (C) 2018 - 2022 by the authors of the ASPECT code.

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


#ifndef _aspect_stokes_matrix_free_h
#define _aspect_stokes_matrix_free_h

#include <aspect/global.h>

#include <aspect/simulator.h>

#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>
#include <deal.II/matrix_free/fe_evaluation.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/multigrid/multigrid.h>
#include <deal.II/multigrid/mg_transfer_matrix_free.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_matrix.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/la_parallel_block_vector.h>

/**
 * Typedef for the number type for the multigrid operators. Can be either float or double.
 */
using GMGNumberType = double;

namespace aspect
{
  using namespace dealii;

  namespace internal
  {
    namespace TangentialBoundaryFunctions
    {

      template <int dim>
      void compute_no_normal_flux_constraints_box (const DoFHandler<dim>    &dof,
                                                   const types::boundary_id  bid,
                                                   const unsigned int first_vector_component,
                                                   MGConstrainedDoFs         &mg_constrained_dofs);
      template <int dim>
      void
      add_constraint(const std::array<types::global_dof_index,dim> &dof_indices,
                     const Tensor<1, dim> &constraining_vector,
                     AffineConstraints<double> &constraints,
                     const double inhomogeneity = 0);

      template <int dim, int spacedim>
      void compute_no_normal_flux_constraints_shell(const DoFHandler<dim,spacedim> &dof_handler,
                                                    const MGConstrainedDoFs        &mg_constrained_dofs,
                                                    const Mapping<dim> &mapping,
                                                    const unsigned int level,
                                                    const unsigned int first_vector_component,
                                                    const std::set<types::boundary_id> &boundary_ids,
                                                    AffineConstraints<double> &constraints);
    }

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
     * This struct stores the data for the current linear operator that is requried to perform
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
       * Table which stores the product of the viscosity derivative
       * with respect to pressure and the Newton derivative scaling
       * factor alpha.
       */
      Table<2, VectorizedArray<number>> newton_factor_wrt_pressure_table;

      /**
       * Table which stores the product of the following three
       * variables: viscosity derivative with respect to strain rate,
       * newton derivative scaling factor, and alpha. Here alpha is
       * the spd factor when the stabilization is PD or SPD,
       * otherwise, it is 1.
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
         * Computes the diagonal contribution from a cell matrix.
         */
        void local_compute_diagonal (const MatrixFree<dim,number>                     &data,
                                     dealii::LinearAlgebra::distributed::Vector<number>  &dst,
                                     const unsigned int                               &dummy,
                                     const std::pair<unsigned int,unsigned int>       &cell_range) const;

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
         * Computes the diagonal contribution from a cell matrix.
         */
        void local_compute_diagonal (const MatrixFree<dim,number>                     &data,
                                     dealii::LinearAlgebra::distributed::Vector<number>  &dst,
                                     const unsigned int                               &dummy,
                                     const std::pair<unsigned int,unsigned int>       &cell_range) const;

        /**
         * A pointer to the current cell data that contains viscosity and other required parameters per cell.
         */
        const OperatorCellData<dim,number> *cell_data;
    };
  }

  /**
   * Base class for the matrix free GMG solver for the Stokes system. The
   * actual implementation is found inside StokesMatrixFreeHandlerImplementation below.
   */
  template<int dim>
  class StokesMatrixFreeHandler
  {
    public:
      /**
       * virtual Destructor.
       */
      virtual ~StokesMatrixFreeHandler() = default;

      /**
       * Solves the Stokes linear system matrix-free. This is called
       * by Simulator<dim>::solve_stokes().
       */
      virtual std::pair<double,double> solve()=0;

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
       * Declare parameters.
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
      virtual const MGTransferMatrixFree<dim,GMGNumberType> &
      get_mg_transfer_A () const = 0;

      /**
       * Return a pointer to the MGTransfer object used for the Schur
       * complement block of the block GMG Stokes solver.
       */
      virtual const MGTransferMatrixFree<dim,GMGNumberType> &
      get_mg_transfer_S () const = 0;

      /**
       * Return the memory consumption in bytes that are used to store
       * equation data like viscosity to be able to apply the operators.
       */
      virtual std::size_t get_cell_data_memory_consumption() const = 0;
  };

  /**
   * Main class of the Matrix-free method. Here are all the functions for
   * setup, assembly and solving the Stokes system.
   *
   * We need to derive from StokesMatrixFreeHandler to be able to introduce a
   * second template argument for the degree of the Stokes finite
   * element. This way, the main simulator does not need to know about the
   * degree by using a pointer to the base class and we can pick the desired
   * velocity degree at runtime.
   */
  template<int dim, int velocity_degree>
  class StokesMatrixFreeHandlerImplementation: public StokesMatrixFreeHandler<dim>
  {
    public:
      /**
       * Initialize this class, allowing it to read in
       * relevant parameters as well as giving it a reference to the
       * Simulator that owns it, since it needs to make fairly extensive
       * changes to the internals of the simulator.
       */
      StokesMatrixFreeHandlerImplementation(Simulator<dim> &, ParameterHandler &prm);

      /**
       * Destructor.
       */
      ~StokesMatrixFreeHandlerImplementation() override = default;

      /**
       * Solves the Stokes linear system matrix-free. This is called
       * by Simulator<dim>::solve_stokes().
       */
      std::pair<double,double> solve() override;

      /**
       * Allocates and sets up the members of the StokesMatrixFreeHandler. This
       * is called by Simulator<dim>::setup_dofs()
       */
      void setup_dofs() override;

      /**
       * Perform various tasks to update the linear system to solve
       * for. Note that we are not assembling a matrix (as this is a
       * matrix-free algorithm), but we are evaluating the material
       * model and storing the information necessary for a later call
       * to solve().
       */
      virtual void assemble() override;

      /**
       * Computes and sets the diagonal for both the mass matrix operator and the A-block
       * operators on each level for the purpose of smoothing inside the multigrid v-cycle.
       */
      void build_preconditioner() override;

      /**
       * Declare parameters. (No actual parameters at the moment).
       */
      static
      void declare_parameters (ParameterHandler &prm);

      /**
       * Return a reference to the DoFHandler that is used for velocity in
       * the block GMG solver.
       */
      const DoFHandler<dim> &
      get_dof_handler_v () const override;

      /**
       * Return a reference to the DoFHandler that is used for pressure in
       * the block GMG solver.
       */
      const DoFHandler<dim> &
      get_dof_handler_p () const override;

      /**
       * Return a reference to the DoFHandler that is used for the coefficient
       * projection in the block GMG solver.
       */
      const DoFHandler<dim> &
      get_dof_handler_projection () const override;

      /**
       * Return a pointer to the object that describes the velocity DoF
       * constraints for the block GMG Stokes solver.
       */
      const AffineConstraints<double> &
      get_constraints_v () const override;

      /**
       * Return a pointer to the object that describes the pressure DoF
       * constraints for the block GMG Stokes solver.
       */
      const AffineConstraints<double> &
      get_constraints_p () const override;

      /**
       * Return a pointer to the MGTransfer object used for the A block
       * of the block GMG Stokes solver.
       */
      const MGTransferMatrixFree<dim,GMGNumberType> &
      get_mg_transfer_A () const override;

      /**
       * Return a pointer to the MGTransfer object used for the Schur
       * complement block of the block GMG Stokes solver.
       */
      const MGTransferMatrixFree<dim,GMGNumberType> &
      get_mg_transfer_S () const override;


      /**
       * Return the memory consumption in bytes that are used to store
       * equation data like viscosity to be able to apply the operators.
       */
      std::size_t get_cell_data_memory_consumption() const override;

    private:
      /**
       * Parse parameters.
       */
      void parse_parameters (ParameterHandler &prm);

      /**
       * Evaluate the MaterialModel to query information like the viscosity and
       * project this viscosity to the multigrid hierarchy. Also queries
       * other parameters like pressure scaling.
       */
      void evaluate_material_model();

      /**
       * Add correction to system RHS for non-zero boundary condition. See description in
       * StokesMatrixFreeHandler::correct_stokes_rhs() for more information.
       */
      void correct_stokes_rhs();


      Simulator<dim> &sim;

      bool print_details;

      /**
       * If true, it will time the key components of this matrix-free implementation, such as
       * vmult of different matrices, solver IDR with the cheap preconditioner, etc.
       */
      bool do_timings;

      /**
       * The max/min of the evaluated viscosities.
       */
      double minimum_viscosity;
      double maximum_viscosity;

      DoFHandler<dim> dof_handler_v;
      DoFHandler<dim> dof_handler_p;
      DoFHandler<dim> dof_handler_projection;

      FESystem<dim> fe_v;
      FESystem<dim> fe_p;
      FESystem<dim> fe_projection;

      /**
       * Store the data for the Stokes operator (viscosity, etc.) for the active cells.
       */
      MatrixFreeStokesOperators::OperatorCellData<dim, GMGNumberType> active_cell_data;

      /**
       * Store the data for the Stokes operator (viscosity, etc.) for each multigrid level.
       */
      MGLevelObject<MatrixFreeStokesOperators::OperatorCellData<dim, GMGNumberType>> level_cell_data;

      // This variable is needed only in the setup in both evaluate_material_model()
      // and build_preconditioner(). It will be deleted after the last use.
      MGLevelObject<dealii::LinearAlgebra::distributed::Vector<GMGNumberType>> level_viscosity_vector;

      using StokesMatrixType = MatrixFreeStokesOperators::StokesOperator<dim,velocity_degree,double>;
      using SchurComplementMatrixType = MatrixFreeStokesOperators::MassMatrixOperator<dim,velocity_degree-1,double>;
      using ABlockMatrixType = MatrixFreeStokesOperators::ABlockOperator<dim,velocity_degree,double>;

      using GMGSchurComplementMatrixType = MatrixFreeStokesOperators::MassMatrixOperator<dim,velocity_degree-1,GMGNumberType>;
      using GMGABlockMatrixType = MatrixFreeStokesOperators::ABlockOperator<dim,velocity_degree,GMGNumberType>;

      StokesMatrixType stokes_matrix;
      ABlockMatrixType A_block_matrix;
      SchurComplementMatrixType Schur_complement_block_matrix;

      AffineConstraints<double> constraints_v;
      AffineConstraints<double> constraints_p;

      MGLevelObject<GMGABlockMatrixType> mg_matrices_A_block;
      MGLevelObject<GMGSchurComplementMatrixType> mg_matrices_Schur_complement;

      MGConstrainedDoFs mg_constrained_dofs_A_block;
      MGConstrainedDoFs mg_constrained_dofs_Schur_complement;
      MGConstrainedDoFs mg_constrained_dofs_projection;

      MGTransferMatrixFree<dim,GMGNumberType> mg_transfer_A_block;
      MGTransferMatrixFree<dim,GMGNumberType> mg_transfer_Schur_complement;

      std::vector<std::shared_ptr<MatrixFree<dim,double>>> matrix_free_objects;
  };
}


#endif
