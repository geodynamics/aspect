/*
  Copyright (C) 2011 - 2019 by the authors of the ASPECT code.

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


#ifndef _aspect_global_h
#define _aspect_global_h

#include <aspect/config.h>

#include <deal.II/base/mpi.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS

#include <deal.II/lac/generic_linear_algebra.h>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

DEAL_II_ENABLE_EXTRA_DIAGNOSTICS


#include <aspect/compat.h>

namespace aspect
{
  /**
   * The following are a set of global constants which may be used by ASPECT:
   * (for sources of data and values used by ASPECT, see source/global.cc)
   */
  namespace constants
  {
    /**
     * Number of seconds in a year [s]
     */
    constexpr double year_in_seconds = 60*60*24*365.2425;

    /**
     * Zero degrees Celsius to Kelvin [K]
     */
    constexpr double celsius_to_kelvin = 273.15;

    /**
     * Gas constant (also known as R) [J K^-1 mol^-1]
     */
    constexpr double gas_constant = 8.3144621;
    /**
     * Avogadro's constant [mol^-1]
     */
    constexpr double avogadro = 6.02214129e23;
    /**
     * Gravitational constant [m^3 kg^-1 s^-2]
     */
    constexpr double big_g = 6.67384e-11;

    /**
     * Constants for Earth:
     */
    namespace earth
    {

      /**
       * Masses are taken from Yoder (1995)
       */
      namespace masses
      {
        /**
         * Planet mass [kg]
         */
        constexpr double planet = 5.9736e24;
        /**
         * Mass of the whole core [kg]
         */
        constexpr double core = 1.932e24;
        /**
         * Mass of the mantle [kg]
         */
        constexpr double mantle = 4.043e24;
      }

      /**
       * Earth structure radii taken from the IASP91 model
       */
      namespace iasp91_radii
      {
        /**
         * Inner core radius [m], equivalent of 5150 km depth
         */
        constexpr double inner_core = 1.2171e6;
        /**
         * Inner core radius [m], equivalent of 2889 km depth
         */
        constexpr double core = 3.482e6;
        /**
         * Lower mantle radius [m], equivalent of 660 km depth
         */
        constexpr double lower_mantle = 5.711e6;
        /**
         * Radius [m], equivalent of 5150 km depth
         */
        constexpr double planet = 6.371e6;
      }

      /**
       * Gravity values taken from the PREM (Dziewonski and Anderson, 1981)
       */
      namespace prem_gravity
      {
        /**
         * Inner core boundary gravity [ms^-2]
         */
        constexpr double icb = 4.4002;
        /**
         * Core-mantle boundary gravity [ms^-2]
         */
        constexpr double cmb = 10.6823;
        /**
         * Upper-lower mantle boundary gravity [ms^-2]
         */
        constexpr double ulmb = 10.0143;
        /**
         * Surface gravity [ms^-2]
         */
        constexpr double surface = 9.8156;
      }

      /**
       * NIST "Standard gravity" (average gravitational acceleration at surface
       * [ms^-2]
       */
      constexpr double surface_gravity = 9.80665;
    }

    /**
     * Constants for Mars
     */
    namespace mars
    {
      /**
       * Mars structure radii
       */
      namespace radii
      {
        /**
         * Planetary radius from Seidermann et al., 2007 [m]
         */
        constexpr double planet = 3.3895e6;
        /**
         * Core radius from Rivoldini et al., 2011 [m]
         */
        constexpr double core = 1.794e6;
      }
      /**
       * Surface gravity from Lodders et al., 1998 [ms^-2]
       */
      constexpr double surface_gravity = 3.711;
    }
  }

  /**
   * Number of seconds in a year [s] (deprecated)
   */
  using constants::year_in_seconds;


  /**
   * A typedef that denotes the BOOST stream type for reading data during
   * serialization. The type chosen here is a binary archive which we
   * subsequently will have to un-compress.
   */
  using iarchive = boost::archive::binary_iarchive;

  /**
   * A typedef that denotes the BOOST stream type for writing data during
   * serialization. The type chosen here is a binary archive which we compress
   * before writing it into a file.
   */
  using oarchive = boost::archive::binary_oarchive;

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
    using Vector = dealii::PETScWrappers::MPI::Vector;

    /**
     * Typedef for the type used to describe vectors that consist of multiple
     * blocks.
     */
    using BlockVector = dealii::PETScWrappers::MPI::BlockVector;

    /**
     * Typedef for the sparse matrix type used.
     */
    using SparseMatrix = dealii::PETScWrappers::MPI::SparseMatrix;

    /**
     * Typedef for the type used to describe sparse matrices that consist of
     * multiple blocks.
     */
    using BlockSparseMatrix = dealii::PETScWrappers::MPI::BlockSparseMatrix;

    /**
     * Typedef for the base class for all preconditioners.
     */
    using PreconditionBase = dealii::PETScWrappers::PreconditionerBase;

    /**
     * Typedef for the AMG preconditioner type used for the top left block of
     * the Stokes matrix.
     */
    using PreconditionAMG = dealii::PETScWrappers::PreconditionBoomerAMG;

    /**
     * Typedef for the Incomplete Cholesky preconditioner used for other
     * blocks of the system matrix.
     */
    using PreconditionIC = dealii::PETScWrappers::PreconditionICC;

    /**
     * Typedef for the Incomplete LU decomposition preconditioner used for
     * other blocks of the system matrix. Note that PETSc does not support a
     * communicating ILU, so we use Jacobi here.
     */
    using PreconditionILU = dealii::PETScWrappers::PreconditionBlockJacobi;

    /**
     * Typedef for the Jacobi preconditioner used for free surface velocity
     * projection.
     */
    using PreconditionJacobi = dealii::PETScWrappers::PreconditionJacobi;

    /**
     * Typedef for the block compressed sparsity pattern type.
     */
    using BlockDynamicSparsityPattern = dealii::BlockDynamicSparsityPattern;

    /**
     * Typedef for the compressed sparsity pattern type.
     */
    using DynamicSparsityPattern = dealii::DynamicSparsityPattern;
#else
    /**
     * Typedef for the vector type used.
     */
    using Vector = dealii::TrilinosWrappers::MPI::Vector;

    /**
     * Typedef for the type used to describe vectors that consist of multiple
     * blocks.
     */
    using BlockVector = dealii::TrilinosWrappers::MPI::BlockVector;

    /**
     * Typedef for the sparse matrix type used.
     */
    using SparseMatrix = dealii::TrilinosWrappers::SparseMatrix;

    /**
     * Typedef for the type used to describe sparse matrices that consist of
     * multiple blocks.
     */
    using BlockSparseMatrix = dealii::TrilinosWrappers::BlockSparseMatrix;

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
     * Typedef for the Incomplete Cholesky preconditioner used for other
     * blocks of the system matrix.
     */
    using PreconditionIC = dealii::TrilinosWrappers::PreconditionIC;

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
#endif
  }
}


/**
 * Print a header into the given stream that will be written both to screen
 * and to the log file and that provides basic information about what is
 * running, with how many processes, and using which linear algebra library.
 */
template <class Stream>
void print_aspect_header(Stream &stream);

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
