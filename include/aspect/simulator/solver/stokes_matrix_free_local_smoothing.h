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

#ifndef _aspect_simulator_solver_stokes_matrix_free_local_smoothing_h
#define _aspect_simulator_solver_stokes_matrix_free_local_smoothing_h

#include <aspect/global.h>
#include <aspect/simulator/solver/stokes_matrix_free.h>
#include <aspect/simulator/solver/matrix_free_operators.h>

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

namespace aspect
{
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
  template <int dim, int velocity_degree>
  class StokesMatrixFreeHandlerLocalSmoothingImplementation: public StokesMatrixFreeHandler<dim>
  {
    public:
      /**
       * Constructor. Give it a reference to the
       * Simulator that owns it, since it needs to make fairly extensive
       * changes to the internals of the simulator.
       */
      StokesMatrixFreeHandlerLocalSmoothingImplementation(Simulator<dim> &simulator,
                                                          const Parameters<dim> &parameters);

      /**
       * Destructor.
       */
      ~StokesMatrixFreeHandlerLocalSmoothingImplementation() override = default;

      /**
       * Initialize the matrix-free solver.
       */
      void initialize() override;

      /**
       * Return the name of the solver for screen output.
       */
      std::string name() const override;

      /**
       * Solves the Stokes linear system using the matrix-free
       * solver.
       *
       * @param system_matrix The system matrix. Note that we do not actually
       * use this matrix for this matrix free solver.
       * @param system_rhs The right hand side vector of the system.
       * @param solve_newton_system A flag indicating whether the system to be
       * solved is the normal linear system or the Newton system. If the Newton
       * system is solved, some operations have to change, e.g. the residual
       * is computed differently.
       * @param last_pressure_normalization_adjustment The amount by which the
       * pressure was adjusted to satisfy the chosen pressure normalization. This
       * information is used to undo the normalization before the solve.
       * @param solution_vector The solution vector that will be
       * updated with the new solution. This vector is expected to have the
       * block structure of the full solution vector, and its velocity and
       * pressure blocks will be updated with the new solution.
       *
       * @return A structure that contains information about the solver, like
       * the initial and final residual.
       */
      StokesSolver::SolverOutputs
      solve(const LinearAlgebra::BlockSparseMatrix &system_matrix,
            const LinearAlgebra::BlockVector &system_rhs,
            const bool solve_newton_system,
            const double last_pressure_normalization_adjustment,
            LinearAlgebra::BlockVector &solution_vector) override;

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
      void assemble() override;

      /**
       * Computes and sets the diagonal for both the mass matrix operator and the A-block
       * operators on each level for the purpose of smoothing inside the multigrid v-cycle.
       */
      void build_preconditioner() override;

      /**
       * Declare parameters.
       */
      static
      void declare_parameters (ParameterHandler &prm);

      /**
       * Parse parameters.
       */
      void parse_parameters (ParameterHandler &prm) override;

      /**
       * Return memory consumption in bytes for all DoFHandler objects.
       */
      std::size_t get_dof_handler_memory_consumption() const override;

      /**
       * Return memory consumption in bytes for all transfer objects.
       */
      std::size_t get_mg_transfer_memory_consumption() const override;

      /**
       * Return memory consumption in bytes for all transfer objects.
       */
      std::size_t get_constraint_memory_consumption() const override;

      /**
       * Return the memory consumption in bytes that are used to store
       * equation data like viscosity to be able to apply the operators.
       */
      std::size_t get_cell_data_memory_consumption() const override;

    private:
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

      using StokesMatrixType = MatrixFreeStokesOperators::StokesOperator<dim,velocity_degree,double>;
      using SchurComplementMatrixType = MatrixFreeStokesOperators::MassMatrixOperator<dim,velocity_degree-1,double>;
      using ABlockMatrixType = MatrixFreeStokesOperators::ABlockOperator<dim,velocity_degree,double>;
      using BTBlockOperatorType = MatrixFreeStokesOperators::BTBlockOperator<dim,velocity_degree,double>;
      using GMGSchurComplementMatrixType = MatrixFreeStokesOperators::MassMatrixOperator<dim,velocity_degree-1,GMGNumberType>;
      using GMGABlockMatrixType = MatrixFreeStokesOperators::ABlockOperator<dim,velocity_degree,GMGNumberType>;

      StokesMatrixType stokes_matrix;
      ABlockMatrixType A_block_matrix;
      BTBlockOperatorType BT_block;
      SchurComplementMatrixType Schur_complement_block_matrix;

      AffineConstraints<double> constraints_v;
      AffineConstraints<double> constraints_p;

      MGLevelObject<GMGABlockMatrixType> mg_matrices_A_block;
      MGLevelObject<GMGSchurComplementMatrixType> mg_matrices_Schur_complement;

      MGConstrainedDoFs mg_constrained_dofs_A_block;
      MGConstrainedDoFs mg_constrained_dofs_Schur_complement;
      MGConstrainedDoFs mg_constrained_dofs_projection;

      MGTransferMF<dim,GMGNumberType> mg_transfer_A_block;
      MGTransferMF<dim,GMGNumberType> mg_transfer_Schur_complement;

      std::vector<std::shared_ptr<MatrixFree<dim,double>>> matrix_free_objects;
  };
}

#endif
