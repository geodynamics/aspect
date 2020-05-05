/*
  Copyright (C) 2018 - 2019 by the authors of the ASPECT code.

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

namespace aspect
{
  using namespace dealii;

  /**
   * This namespace contains all matrix-free operators used in the Stokes solver.
   */
  namespace MatrixFreeStokesOperators
  {
    /**
     * Operator for the entire Stokes block.
     */
    template <int dim, int degree_v, typename number>
    class StokesOperator
      : public MatrixFreeOperators::Base<dim, dealii::LinearAlgebra::distributed::BlockVector<number> >
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
         * Fills in the viscosity table, sets the value for the pressure scaling constant,
         * and gives information regarding compressibility.
         */
        void fill_cell_data (const Table<1, VectorizedArray<number>> &viscosity_table,
                             const double pressure_scaling,
                             const bool is_compressible);

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
         * Table which stores a viscosity value for each cell.
         */
        const Table<1, VectorizedArray<number>> *viscosity;

        /**
         * Pressure scaling constant.
         */
        double pressure_scaling;

        /**
          * Information on the compressibility of the flow.
          */
        bool is_compressible;
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
         * Fills in the viscosity table and sets the value for the pressure scaling constant. The input
         * @p is_mg_level_data describes whether the viscosity values are defined for a multigrid level
         * matrix or for the active level matrix.
         */
        void fill_cell_data (const Table<1, VectorizedArray<number>> &viscosity_table,
                             const double pressure_scaling);


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
         * Table which stores a viscosity value for each cell.
         */
        const Table<1, VectorizedArray<number>> *viscosity;

        /**
         * Pressure scaling constant.
         */
        double pressure_scaling;
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
         * Fills in the viscosity table and gives information regarding compressibility. The input
         * @p is_mg_level_data describes whether the viscosity values are defined for a multigrid level
         * matrix or for the active level matrix.
         */
        void fill_cell_data (const Table<1, VectorizedArray<number>> &viscosity_table,
                             const bool is_compressible);

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
         * Table which stores a viscosity value for each cell.
         */
        const Table<1, VectorizedArray<number>> *viscosity;

        /**
          * Information on the compressibility of the flow.
          */
        bool is_compressible;

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
       * Evalute the MaterialModel to query for the viscosity on the active cells,
       * project this viscosity to the multigrid hierarchy, and cache the information
       * for later usage. Also sets pressure scaling and information regarding the
       * compressiblity of the flow.
       */
      virtual void evaluate_material_model()=0;

      /**
       * Add correction to system RHS for non-zero boundary condition.
       */
      virtual void correct_stokes_rhs()=0;

      /**
       * Computes and sets the diagonal for both the mass matrix operator and the A-block
       * operators on each level for the purpose of smoothing inside the multigrid v-cycle.
       */
      virtual void build_preconditioner()=0;

      /**
       * Declare parameters.
       */
      static
      void declare_parameters (ParameterHandler &prm);
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
       * Evalute the MaterialModel to query for the viscosity on the active cells,
       * project this viscosity to the multigrid hierarchy, and cache the information
       * for later usage. Also sets pressure scaling and information regarding the
       * compressiblity of the flow.
       */
      void evaluate_material_model() override;

      /**
       * Add correction to system RHS for non-zero boundary condition.
       */
      void correct_stokes_rhs() override;

      /**
       * Computes and sets the diagonal for both the mass matrix operator and the A-block
       * operators on each level for the purpose of smoothing inside the multigrid v-cycle.
       */
      void build_preconditioner() override;

      /**
       * Get the workload imbalance of the distribution
       * of the level hierarchy.
       */
      double get_workload_imbalance();

      /**
       * Declare parameters. (No actual parameters at the moment).
       */
      static
      void declare_parameters (ParameterHandler &prm);

    private:
      /**
       * Parse parameters. (No actual parameters at the moment).
       */
      void parse_parameters (ParameterHandler &prm);


      Simulator<dim> &sim;

      DoFHandler<dim> dof_handler_v;
      DoFHandler<dim> dof_handler_p;
      DoFHandler<dim> dof_handler_projection;

      FESystem<dim> fe_v;
      FESystem<dim> fe_p;
      FESystem<dim> fe_projection;

      Table<1, VectorizedArray<double>> active_viscosity_table;
      MGLevelObject<Table<1, VectorizedArray<double>>> level_viscosity_tables;

      // This variable is needed only in the setup in both evaluate_material_model()
      // and build_preconditioner(). It will be deleted after the last use.
      MGLevelObject<dealii::LinearAlgebra::distributed::Vector<double> > level_viscosity_vector;

      typedef MatrixFreeStokesOperators::StokesOperator<dim,velocity_degree,double> StokesMatrixType;
      typedef MatrixFreeStokesOperators::MassMatrixOperator<dim,velocity_degree-1,double> SchurComplementMatrixType;
      typedef MatrixFreeStokesOperators::ABlockOperator<dim,velocity_degree,double> ABlockMatrixType;

      StokesMatrixType stokes_matrix;
      ABlockMatrixType A_block_matrix;
      SchurComplementMatrixType Schur_complement_block_matrix;

      ConstraintMatrix constraints_v;
      ConstraintMatrix constraints_p;

      MGLevelObject<ABlockMatrixType> mg_matrices_A_block;
      MGLevelObject<SchurComplementMatrixType> mg_matrices_Schur_complement;

      MGConstrainedDoFs mg_constrained_dofs_A_block;
      MGConstrainedDoFs mg_constrained_dofs_Schur_complement;
      MGConstrainedDoFs mg_constrained_dofs_projection;

      MGTransferMatrixFree<dim,double> mg_transfer_A_block;
      MGTransferMatrixFree<dim,double> mg_transfer_Schur_complement;
  };
}


#endif
