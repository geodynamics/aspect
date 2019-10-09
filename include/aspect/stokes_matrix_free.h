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
         * Reset the viscosity table.
         */
        void clear ();

        /**
         * Fills in the viscosity table, sets the value for the pressure scaling constant,
         * and gives information regarding compressibility.
         */
        void fill_cell_data(const dealii::LinearAlgebra::distributed::Vector<number> &viscosity_values,
                            const double pressure_scaling,
                            const Triangulation<dim> &tria,
                            const DoFHandler<dim> &dof_handler_for_projection,
                            const bool is_compressible);

        /**
         * Returns the viscosity table.
         */
        const Table<2, VectorizedArray<number> > &
        get_viscosity_x_2_table();

        /**
         * Computes the diagonal of the matrix. Since matrix-free operators have not access
         * to matrix elements, we must apply the matrix-free operator to the unit vectors to
         * recover the diagonal.
         */
        virtual void compute_diagonal ();

      private:

        /**
         * Performs the application of the matrix-free operator. This function is called by
         * vmult() functions MatrixFreeOperators::Base.
         */
        virtual void apply_add (dealii::LinearAlgebra::distributed::BlockVector<number> &dst,
                                const dealii::LinearAlgebra::distributed::BlockVector<number> &src) const;

        /**
         * Defines the application of the cell matrix.
         */
        void local_apply (const dealii::MatrixFree<dim, number> &data,
                          dealii::LinearAlgebra::distributed::BlockVector<number> &dst,
                          const dealii::LinearAlgebra::distributed::BlockVector<number> &src,
                          const std::pair<unsigned int, unsigned int> &cell_range) const;

        /**
         * Table which stores the viscosity on each quadrature point.
         */
        Table<2, VectorizedArray<number> > viscosity_x_2;

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
         * Reset the viscosity table.
         */
        void clear ();

        /**
         * Fills in the viscosity table and sets the value for the pressure scaling constant.
         */
        void fill_cell_data (const dealii::LinearAlgebra::distributed::Vector<number> &viscosity_values,
                             const double pressure_scaling,
                             const Triangulation<dim> &tria,
                             const DoFHandler<dim> &dof_handler_for_projection);


        /**
         * Computes the diagonal of the matrix. Since matrix-free operators have not access
         * to matrix elements, we must apply the matrix-free operator to the unit vectors to
         * recover the diagonal.
         */
        virtual void compute_diagonal ();

      private:

        /**
         * Performs the application of the matrix-free operator. This function is called by
         * vmult() functions MatrixFreeOperators::Base.
         */
        virtual void apply_add (dealii::LinearAlgebra::distributed::Vector<number> &dst,
                                const dealii::LinearAlgebra::distributed::Vector<number> &src) const;

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
         * Table which stores the viscosity on each quadrature point.
         */
        Table<2, VectorizedArray<number> > one_over_viscosity;

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
         * Reset the viscosity table.
         */
        void clear ();

        /**
         * Fills in the viscosity table and gives information regarding compressibility.
         */
        void fill_cell_data(const dealii::LinearAlgebra::distributed::Vector<number> &viscosity_values,
                            const Triangulation<dim> &tria,
                            const DoFHandler<dim> &dof_handler_for_projection,
                            const bool for_mg,
                            const bool is_compressible);

        /**
         * Computes the diagonal of the matrix. Since matrix-free operators have not access
         * to matrix elements, we must apply the matrix-free operator to the unit vectors to
         * recover the diagonal.
         */
        virtual void compute_diagonal ();

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
        virtual void apply_add (dealii::LinearAlgebra::distributed::Vector<number> &dst,
                                const dealii::LinearAlgebra::distributed::Vector<number> &src) const;

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
         * Table which stores the viscosity on each quadrature point.
         */
        Table<2, VectorizedArray<number> > viscosity_x_2;

        /**
          * Information on the compressibility of the flow.
          */
        bool is_compressible;

    };
  }

  /**
   * Main class of the Matrix-free method. Here are all the functions for setup, assembly and solving
   * the Stokes system.
   */
  template<int dim>
  class StokesMatrixFreeHandler
  {
    public:
      /**
       * Initialize this class, allowing it to read in
       * relevant parameters as well as giving it a reference to the
       * Simulator that owns it, since it needs to make fairly extensive
       * changes to the internals of the simulator.
       */
      StokesMatrixFreeHandler(Simulator<dim> &, ParameterHandler &prm);

      /**
       * Destructor.
       */
      ~StokesMatrixFreeHandler();

      /**
       * Solves the Stokes linear system matrix-free. This is called
       * by Simulator<dim>::solve_stokes().
       */
      std::pair<double,double> solve();

      /**
       * Allocates and sets up the members of the StokesMatrixFreeHandler. This
       * is called by Simulator<dim>::setup_dofs()
       */
      void setup_dofs();

      /**
       * Evaluate the material model and update internal data structures before the
       * actual solve().
       */
      void build_preconditioner();

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

      /**
       * Parse parameters. (No actual parameters at the moment).
       */
      void parse_parameters (ParameterHandler &prm);

    private:

      /**
       * Evalute the MaterialModel to query for the viscosity on the active cells,
       * project this viscosity to the multigrid hierarchy, and cache the information
       * for later usage. Also sets pressure scaling and information regarding the
       * compressiblity of the flow.
       */
      void evaluate_material_model();

      /**
       * Add correction to system RHS for non-zero boundary condition.
       */
      void correct_stokes_rhs();

      /**
       * Computes and sets the diagonal for the A-block operators on each level for
       * the purpose of smoothing inside the multigrid v-cycle.
       */
      void compute_A_block_diagonals();

      Simulator<dim> &sim;

      DoFHandler<dim> dof_handler_v;
      DoFHandler<dim> dof_handler_p;
      DoFHandler<dim> dof_handler_projection;

      FESystem<dim> stokes_fe;
      FESystem<dim> fe_v;
      FESystem<dim> fe_p;
      FESystem<dim> fe_projection;

      // TODO: velocity degree not only 2, Choosing quadrature degree?
      typedef MatrixFreeStokesOperators::StokesOperator<dim,2,double> StokesMatrixType;
      typedef MatrixFreeStokesOperators::MassMatrixOperator<dim,1,double> MassMatrixType;
      typedef MatrixFreeStokesOperators::ABlockOperator<dim,2,double> ABlockMatrixType;

      StokesMatrixType stokes_matrix;
      ABlockMatrixType velocity_matrix;
      MassMatrixType mass_matrix;

      ConstraintMatrix constraints_v;
      ConstraintMatrix constraints_p;
      ConstraintMatrix constraints_projection;

      MGLevelObject<ABlockMatrixType> mg_matrices;
      MGConstrainedDoFs              mg_constrained_dofs;
      MGConstrainedDoFs mg_constrained_dofs_projection;

      dealii::LinearAlgebra::distributed::Vector<double> active_coef_dof_vec;
      MGLevelObject<dealii::LinearAlgebra::distributed::Vector<double> > level_coef_dof_vec;

      MGTransferMatrixFree<dim,double> mg_transfer;
  };
}


#endif
