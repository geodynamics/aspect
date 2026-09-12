/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

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

#ifndef _aspect_adjoint_state_h
#define _aspect_adjoint_state_h

#include <aspect/global.h>

#include <string>

namespace aspect
{
  namespace Adjoint
  {
    /**
     * A data structure that stores the value and right hand side contribution
     * of one objective before the adjoint linear solve.
     *
     * @ingroup Adjoint
     */
    template <int dim>
    struct ObjectiveResult
    {
      /**
       * Objective functional name used in output and repositories.
       */
      std::string objective_name;

      /**
       * Scalar objective value.
       */
      double value = 0.0;

      /**
       * Adjoint right hand side contribution associated with this objective.
       */
      LinearAlgebra::BlockVector rhs;
    };

    /**
     * A data structure that stores the adjoint solution associated with a
     * single objective.
     *
     * @ingroup Adjoint
     */
    template <int dim>
    struct AdjointState
    {
      /**
       * Objective functional name associated with this adjoint solution.
       */
      std::string objective_name;

      /**
       * Solution of the adjoint Stokes system.
       */
      LinearAlgebra::BlockVector solution;
    };

  }
}

#endif
