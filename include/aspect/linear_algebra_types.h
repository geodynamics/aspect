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


#ifndef _aspect_linear_algebra_types_h
#define _aspect_linear_algebra_types_h

#include <aspect/global.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS

#include <deal.II/lac/generic_linear_algebra.h>

DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

namespace aspect
{
  /**
   * A namespace that contains typedefs for classes used in the linear algebra
   * description.
   */
  namespace LinearAlgebra
  {
    /**
     * Typedef for the base class for all preconditioners.
     */
    using PreconditionBase = dealii::TrilinosWrappers::PreconditionBase;

    /**
     * Typedef for the AMG preconditioner type used for the top left block of
     * the Stokes matrix.
     */
    using PreconditionAMG = dealii::TrilinosWrappers::PreconditionAMG;

    /**
     * Typedef for the Incomplete LU decomposition preconditioner used for
     * other blocks of the system matrix.
     */
    using PreconditionILU = dealii::TrilinosWrappers::PreconditionILU;

    /**
     * Typedef for the Jacobi preconditioner used for free surface velocity
     * projection.
     */
    using PreconditionJacobi = dealii::TrilinosWrappers::PreconditionJacobi;

    /**
     * Typedef for the block compressed sparsity pattern type.
     */
    using BlockDynamicSparsityPattern = dealii::TrilinosWrappers::BlockSparsityPattern;

    /**
     * Typedef for the compressed sparsity pattern type.
     */
    using DynamicSparsityPattern = dealii::TrilinosWrappers::SparsityPattern;
  }
}



#endif
