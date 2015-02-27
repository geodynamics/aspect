/*
  Copyright (C) 2011, 2012, 2014 by the authors of the ASPECT code.

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


#ifndef __aspect__global_h
#define __aspect__global_h

#ifdef ASPECT_USE_PETSC
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

#include <deal.II/base/mpi.h>
#include <deal.II/base/multithread_info.h>


namespace aspect
{
  /**
   * A variable whose value denotes the number of seconds in one year.
   */
  extern const double year_in_seconds;

  /**
   * A variable that denotes whether we should periodically output statistics
   * about memory consumption, run times, etc via the
   * Simulator::output_statistics() function or other means.
   */
  extern const bool output_parallel_statistics;


  /**
   * A typedef that denotes the BOOST stream type for reading data during
   * serialization. The type chosen here is a binary archive which we
   * subsequently will have to un-compress.
   */
  typedef boost::archive::binary_iarchive iarchive;

  /**
   * A typedef that denotes the BOOST stream type for writing data during
   * serialization. The type chosen here is a binary archive which we compress
   * before writing it into a file.
   */
  typedef boost::archive::binary_oarchive oarchive;

  /**
   * A class we throw in exceptions in parallel jobs and that we can silently
   * treat in main(). We do this, for example, in read_parameters() where each
   * processor would otherwise throw the same exception and every processor
   * would produce a tangle of output that is impenetrable in large parallel
   * jobs. The same situation happens if a linear solver fails. Rather, we
   * make processor 0 throw the real exception and every other processor
   * converts the exception it wants to throw to an object of the current type
   * -- which is caught in main() but doesn't produce any output (because
   * processor 0 will already produce the output).
   */
  class QuietException {};


  /**
   * A namespace that contains typedefs for classes used in the linear algebra
   * description.
   */
  namespace LinearAlgebra
  {
#ifdef ASPECT_USE_PETSC
    /**
     * Typedef for the vector type used.
     */
    typedef dealii::PETScWrappers::MPI::Vector Vector;

    /**
     * Typedef for the type used to describe vectors that consist of multiple
     * blocks.
     */
    typedef dealii::PETScWrappers::MPI::BlockVector BlockVector;

    /**
     * Typedef for the sparse matrix type used.
     */
    typedef dealii::PETScWrappers::MPI::SparseMatrix SparseMatrix;

    /**
     * Typedef for the type used to describe sparse matrices that consist of
     * multiple blocks.
     */
    typedef dealii::PETScWrappers::MPI::BlockSparseMatrix BlockSparseMatrix;

    /**
     * Typedef for the AMG preconditioner type used for the top left block of
     * the Stokes matrix.
     */
    typedef dealii::PETScWrappers::PreconditionBoomerAMG PreconditionAMG;

    /**
     * Typedef for the Incomplete Cholesky preconditioner used for other
     * blocks of the system matrix.
     */
    typedef dealii::PETScWrappers::PreconditionICC PreconditionIC;

    /**
     * Typedef for the Incomplete LU decomposition preconditioner used for
     * other blocks of the system matrix. Note that PETSc does not support a
     * communicating ILU, so we use Jacobi here.
     */
    typedef dealii::PETScWrappers::PreconditionBlockJacobi PreconditionILU;

    /**
     * Typedef for the Jacobi preconditioner used for free surface velocity
     * projection.
     */
    typedef dealii::PETScWrappers::PreconditionJacobi PreconditionJacobi;

    /**
     * Typedef for the block compressed sparsity pattern type.
     */
    typedef dealii::BlockCompressedSimpleSparsityPattern BlockCompressedSparsityPattern;

#else
    /**
     * Typedef for the vector type used.
     */
    typedef dealii::TrilinosWrappers::MPI::Vector Vector;

    /**
     * Typedef for the type used to describe vectors that consist of multiple
     * blocks.
     */
    typedef dealii::TrilinosWrappers::MPI::BlockVector BlockVector;

    /**
     * Typedef for the sparse matrix type used.
     */
    typedef dealii::TrilinosWrappers::SparseMatrix SparseMatrix;

    /**
     * Typedef for the type used to describe sparse matrices that consist of
     * multiple blocks.
     */
    typedef dealii::TrilinosWrappers::BlockSparseMatrix BlockSparseMatrix;

    /**
     * Typedef for the AMG preconditioner type used for the top left block of
     * the Stokes matrix.
     */
    typedef dealii::TrilinosWrappers::PreconditionAMG PreconditionAMG;

    /**
     * Typedef for the Incomplete Cholesky preconditioner used for other
     * blocks of the system matrix.
     */
    typedef dealii::TrilinosWrappers::PreconditionIC PreconditionIC;

    /**
     * Typedef for the Incomplete LU decomposition preconditioner used for
     * other blocks of the system matrix.
     */
    typedef dealii::TrilinosWrappers::PreconditionILU PreconditionILU;

    /**
     * Typedef for the Jacobi preconditioner used for free surface velocity
     * projection.
     */
    typedef dealii::TrilinosWrappers::PreconditionJacobi PreconditionJacobi;

    /**
     * Typedef for the block compressed sparsity pattern type.
     */
    typedef dealii::TrilinosWrappers::BlockSparsityPattern BlockCompressedSparsityPattern;

#endif
  }
}


template < class Stream>
void print_aspect_header(Stream &stream)
{
  const int n_tasks = dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  stream << "-----------------------------------------------------------------------------\n"
         << "-- This is ASPECT, the Advanced Solver for Problems in Earth's ConvecTion.\n"
         << "--     . version 1.3.pre\n" //VERSION-INFO. Do not edit by hand.
#ifdef DEBUG
         << "--     . running in DEBUG mode\n"
#else
         << "--     . running in OPTIMIZED mode\n"
#endif
         << "--     . running with " << n_tasks << " MPI process" << (n_tasks == 1 ? "\n" : "es\n");
  const int n_threads =
#if DEAL_II_VERSION_GTE(8,3,0)
      dealii::MultithreadInfo::n_threads();
#else
      dealii::multithread_info.n_threads();
#endif
  if (n_threads>1)
    stream << "--     . using " << n_threads << " threads " << (n_tasks == 1 ? "\n" : "each\n");
#ifdef ASPECT_USE_PETSC
  stream << "--     . using PETSc\n";
#else
  stream << "--     . using Trilinos\n";
#endif
  stream << "-----------------------------------------------------------------------------\n"
         << std::endl;
}


/**
 * A macro that is used in instantiating the ASPECT classes and functions for
 * both 2d and 3d. Call this macro with the name of another macro that when
 * called with a single integer argument instantiates the respective classes
 * in the given space dimension.
 */
#define ASPECT_INSTANTIATE(INSTANTIATIONS) \
  INSTANTIATIONS(2) \
  INSTANTIATIONS(3)

#endif
