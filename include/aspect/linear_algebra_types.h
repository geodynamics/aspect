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


#ifndef _aspect_linear_algebra_types_h
#define _aspect_linear_algebra_types_h

#include <aspect/global.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS

#ifdef ASPECT_USE_TPETRA
#include <deal.II/base/memory_space.h>
#include <deal.II/lac/trilinos_tpetra_precondition.h>
#include <deal.II/lac/trilinos_tpetra_sparsity_pattern.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#else
#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/lac/trilinos_precondition.h>
#endif

DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

namespace aspect
{
  /**
   * A namespace that contains typedefs for classes used in the linear algebra
   * description.
   */
  namespace LinearAlgebra
  {
#ifdef ASPECT_USE_TPETRA
    /**
     * Typedef for the base class for all preconditioners.
     */
    using PreconditionBase = dealii::LinearAlgebra::TpetraWrappers::PreconditionBase<double>;

    /**
     * Typedef for the AMG preconditioner type used for the top left block of
     * the Stokes matrix.
     */
    using PreconditionAMG = dealii::LinearAlgebra::TpetraWrappers::PreconditionAMGMueLu<double>;

    /**
     * Typedef for the Incomplete LU decomposition preconditioner used for
     * other blocks of the system matrix.
     */
    using PreconditionILU = dealii::LinearAlgebra::TpetraWrappers::PreconditionILU<double>;

    /**
     * Typedef for the Jacobi preconditioner used for free surface velocity
     * projection.
     */
    using PreconditionJacobi = dealii::LinearAlgebra::TpetraWrappers::PreconditionJacobi<double>;
#else
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
#endif

    /**
     * Typedef for the block compressed sparsity pattern type.
     */
    using BlockDynamicSparsityPattern = dealii::BlockDynamicSparsityPattern;

    /**
     * Typedef for the compressed sparsity pattern type.
     */
    using DynamicSparsityPattern = dealii::DynamicSparsityPattern;
  }
}



#endif
