/*
  Copyright (C) 2011, 2012 by the authors of the ASPECT code.

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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/
/*  $Id$  */


#ifndef __aspect__global_h
#define __aspect__global_h

// uncomment this to use PETSc for linear algebra
//#define USE_PETSC

#ifdef USE_PETSC
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/petsc_block_sparse_matrix.h>
#include <deal.II/lac/petsc_precondition.h>
#else
#include <deal.II/lac/trilinos_block_vector.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_precondition.h>
#endif

#include <deal.II/lac/generic_linear_algebra.h>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>



namespace aspect
{
  /**
   * A variable whose value denotes the number of seconds in one year.
   */
  extern const double year_in_seconds;

  /**
   * A variable that denotes whether we should periodically
   * output statistics about memory consumption, run times, etc
   * via the Simulator::output_statistics() function or other
   * means.
   */
  extern const bool output_parallel_statistics;


  /**
   * A typedef that denotes the BOOST stream type for reading data
   * during serialization. The type chosen here is a binary archive
   * which we subsequently will have to un-compress.
   */
  typedef boost::archive::binary_iarchive iarchive;

  /**
   * A typedef that denotes the BOOST stream type for writing data
   * during serialization. The type chosen here is a binary archive
   * which we compress before writing it into a file.
   */
  typedef boost::archive::binary_oarchive oarchive;

  /**
   * A namespace that contains typedefs for classes used in
   * the linear algebra description.
   */
  namespace LinearAlgebra
  {
    using namespace dealii;

#ifdef USE_PETSC
    /**
     * Typedef for the vector type used.
     */
    typedef PETScWrappers::MPI::Vector Vector;

    /**
     * Typedef for the type used to describe vectors that
     * consist of multiple blocks.
     */
    typedef PETScWrappers::MPI::BlockVector BlockVector;

    /**
     * Typedef for the sparse matrix type used.
     */
    typedef PETScWrappers::MPI::SparseMatrix SparseMatrix;

    /**
     * Typedef for the type used to describe sparse matrices that
     * consist of multiple blocks.
     */
    typedef PETScWrappers::MPI::BlockSparseMatrix BlockSparseMatrix;

    /**
     * Typedef for the AMG preconditioner type used for the
     * top left block of the Stokes matrix.
     */
    typedef PETScWrappers::PreconditionBoomerAMG PreconditionAMG;

    /**
     * Typedef for the Incomplete Cholesky preconditioner used
     * for other blocks of the system matrix.
     */
    typedef PETScWrappers::PreconditionICC PreconditionIC;

    /**
     * Typedef for the Incomplete LU decomposition preconditioner used
     * for other blocks of the system matrix. Note that PETSc does not
     * support a communicating ILU, so we use Jacobi here.
     */
    typedef PETScWrappers::PreconditionBlockJacobi PreconditionILU;

    typedef LinearAlgebraPETSc::MPI::CompressedBlockSparsityPattern CompressedBlockSparsityPattern;

#else
    /**
     * Typedef for the vector type used.
     */
    typedef TrilinosWrappers::MPI::Vector Vector;

    /**
     * Typedef for the type used to describe vectors that
     * consist of multiple blocks.
     */
    typedef TrilinosWrappers::MPI::BlockVector BlockVector;

    /**
     * Typedef for the sparse matrix type used.
     */
    typedef TrilinosWrappers::SparseMatrix SparseMatrix;

    /**
     * Typedef for the type used to describe sparse matrices that
     * consist of multiple blocks.
     */
    typedef TrilinosWrappers::BlockSparseMatrix BlockSparseMatrix;

    /**
     * Typedef for the AMG preconditioner type used for the
     * top left block of the Stokes matrix.
     */
    typedef TrilinosWrappers::PreconditionAMG PreconditionAMG;

    /**
     * Typedef for the Incomplete Cholesky preconditioner used
     * for other blocks of the system matrix.
     */
    typedef TrilinosWrappers::PreconditionIC PreconditionIC;

    /**
     * Typedef for the Incomplete LU decomposition preconditioner used
     * for other blocks of the system matrix.
     */
    typedef TrilinosWrappers::PreconditionILU PreconditionILU;
#endif
  }
}


/**
 * A macro that is used in instantiating the ASPECT classes and functions
 * for both 2d and 3d. Call this macro with the name of another macro that
 * when called with a single integer argument instantiates the respective
 * classes in the given space dimension.
 */
#define ASPECT_INSTANTIATE(INSTANTIATIONS) \
  INSTANTIATIONS(2) \
  INSTANTIATIONS(3)

#endif
