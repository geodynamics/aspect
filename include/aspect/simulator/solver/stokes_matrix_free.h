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

}

#endif
