// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2025 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

#ifndef _aspect_simulator_solver_stokes_direct_h
#define _aspect_simulator_solver_stokes_direct_h

#include <aspect/global.h>
#include <aspect/simulator/solver/interface.h>


namespace aspect
{
  namespace StokesSolver
  {
    /**
     * A direct solver for the Stokes equations using the Trilinos Amesos package.
     * The solver used is Amesos_Klu.
     */
    template <int dim>
    class Direct: public Interface<dim>
    {
      public:
        /**
         * Solves the linear system.
         *
         * @param system_matrix The system matrix. Note that even if the matrix
         * is not assembled (e.g. for matrix free solvers), a reference to the
         * system matrix will be provided to the solver.
         * @param system_rhs The right hand side vector of the system.
         * @param solve_newton_system A flag indicating whether the system to be
         * solved is the normal linear system or the Newton system. If the Newton
         * system is solved, some operations have to change, e.g. the residual
         * is computed differently.
         * @param last_pressure_normalization_adjustment The amount by which the
         * pressure was adjusted to satisfy the chosen pressure normalization. This
         * information is used to undo the normalization before the solve.
         * @param solution_vector The existing solution vector that will be
         * updated with the new solution. This vector is expected to have the
         * block structure of the full solution vector, and the blocks that
         * are to be solved in this solver will be overwritten.
         *
         * @return A structure that contains information about the solver, like
         * the initial and final residual.
         */
        SolverOutputs solve(const LinearAlgebra::BlockSparseMatrix &system_matrix,
                            const LinearAlgebra::BlockVector &system_rhs,
                            const bool solve_newton_system,
                            const double last_pressure_normalization_adjustment,
                            LinearAlgebra::BlockVector &solution_vector) override;

        /**
         * Return the name of the solver for screen output.
         */
        std::string name() const override;
    };
  }
}

#endif
