/*
  Copyright (C) 2025 - by the authors of the ASPECT code.

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

#ifndef _aspect_simulator_solver_interface_h
#define _aspect_simulator_solver_interface_h

#include <aspect/global.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace Solver
  {
    /**
     * Base class for ASPECT solvers.
     */
    template <int dim>
    class Interface: public SimulatorAccess<dim>, public Plugins::InterfaceBase
    {
      public:
        /**
         * Solves the linear system.
         *
         * @param system_matrix The system matrix. Note that even if the matrix
         * is not assembled (e.g. for matrix free solvers), a reference to the
         * system matrix will be provided to the solver.
         * @param system_rhs The right hand side vector of the system.
         * @param solution_vector The existing solution vector that will be
         * updated with the new solution. This vector is expected to have the
         * block structure of the full solution vector, and the blocks that
         * are to be solved in this solver will be overwritten.
         *
         * @return A pair of the initial and final residual of the linear solver.
         */
        virtual
        std::pair<double,double> solve(const LinearAlgebra::BlockSparseMatrix &system_matrix,
                                       const LinearAlgebra::BlockVector &system_rhs,
                                       LinearAlgebra::BlockVector &solution_vector) = 0;

        /**
         * Return the name of the solver for screen output.
         */
        virtual
        std::string name() const = 0;
    };
  }
}

#endif
