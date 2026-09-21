// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

#ifndef _aspect_simulator_solver_stokes_matrix_free_h
#define _aspect_simulator_solver_stokes_matrix_free_h

#include <aspect/global.h>
#include <aspect/parameters.h>
#include <aspect/simulator/solver/interface.h>

namespace aspect
{
  /**
   * Typedef for the number type for the multigrid operators. Can be either float or double.
   */
  using GMGNumberType = double;

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
       * Set up the triangulations used by the multigrid hierarchy.
       */
      virtual void setup_multigrid_hierarchy();

      /**
       * Return the triangulations used by the multigrid hierarchy.
       */
      virtual const std::vector<std::shared_ptr<const Triangulation<dim, dim>>> &
      get_multigrid_triangulations() const = 0;

      /**
       * Perform various tasks to update the linear system to solve
       * for. Note that we are not assembling a matrix (as this is a
       * matrix-free algorithm), but we are evaluating the material
       * model and storing the information necessary for a later call
       * to solve(). The function can modify @p system_rhs to account
       * for boundary conditions as is necessary for matrix-free methods.
       */
      virtual void assemble(LinearAlgebra::BlockVector &system_rhs)=0;

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
       * Return memory consumption in bytes for all DoFHandler objects.
       */
      virtual std::size_t get_dof_handler_memory_consumption() const = 0;

      /**
       * Return memory consumption in bytes for all transfer objects.
       */
      virtual std::size_t get_mg_transfer_memory_consumption() const = 0;

      /**
       * Return memory consumption in bytes for all transfer objects.
       */
      virtual std::size_t get_constraint_memory_consumption() const = 0;

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



}

#endif
