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

#include <aspect/simulator/solver/matrix_free_operators.h>
#include <aspect/global.h>
#include <aspect/parameters.h>
#include <aspect/simulator/solver/interface.h>
#include <aspect/utilities.h>

#include <deal.II/multigrid/mg_transfer_matrix_free.h>
#include <deal.II/multigrid/mg_transfer_global_coarsening.templates.h>

#include <deal.II/lac/solver_bicgstab.h>

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
