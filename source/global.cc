/*
  Copyright (C) 2011 - 2020 by the authors of the ASPECT code.

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


#include <aspect/global.h>
#include <aspect/revision.h>

#include <deal.II/base/multithread_info.h>
#include <deal.II/base/revision.h>
#include <deal.II/base/vectorization.h>

#include <cstring>



template <class Stream>
void print_aspect_header(Stream &stream)
{
  const int n_tasks = dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  stream << "-----------------------------------------------------------------------------\n"
         << "-- This is ASPECT, the Advanced Solver for Problems in Earth's ConvecTion.\n"
         << "--     . version " << ASPECT_PACKAGE_VERSION;
  if (strcmp(ASPECT_GIT_BRANCH,"") != 0)
    stream << " (" << ASPECT_GIT_BRANCH << ", " << ASPECT_GIT_SHORTREV << ")\n";
  else
    stream << "\n";

  stream << "--     . using deal.II " << DEAL_II_PACKAGE_VERSION;
  if (strcmp(DEAL_II_GIT_BRANCH,"") != 0)
    stream << " (" << DEAL_II_GIT_BRANCH << ", " << DEAL_II_GIT_SHORTREV << ")";
  stream << "\n";
  stream << "--     .       with "
#ifdef DEAL_II_WITH_64BIT_INDICES
         << "64"
#else
         << "32"
#endif
         << " bit indices and vectorization level ";
  const unsigned int n_vect_bits =
#if DEAL_II_VERSION_GTE(9,2,0)
    dealii::VectorizedArray<double>::size() * 8 * sizeof(double);
#else
    dealii::VectorizedArray<double>::n_array_elements * 8 * sizeof(double);
#endif

  stream << DEAL_II_COMPILER_VECTORIZATION_LEVEL
         << " (" << n_vect_bits << " bits)\n";

#ifdef ASPECT_USE_PETSC
  stream << "--     . using PETSc "
         << PETSC_VERSION_MAJOR    << '.'
         << PETSC_VERSION_MINOR    << '.'
         << PETSC_VERSION_SUBMINOR << '\n';
#else
  stream << "--     . using Trilinos "
         << DEAL_II_TRILINOS_VERSION_MAJOR    << '.'
         << DEAL_II_TRILINOS_VERSION_MINOR    << '.'
         << DEAL_II_TRILINOS_VERSION_SUBMINOR << '\n';
#endif
  stream << "--     . using p4est "
         << DEAL_II_P4EST_VERSION_MAJOR << '.'
         << DEAL_II_P4EST_VERSION_MINOR << '.'
         << DEAL_II_P4EST_VERSION_SUBMINOR << '\n';

#ifdef DEBUG
  stream << "--     . running in DEBUG mode\n"
#else
  stream << "--     . running in OPTIMIZED mode\n"
#endif
         << "--     . running with " << n_tasks << " MPI process" << (n_tasks == 1 ? "\n" : "es\n");

  const int n_threads =
    dealii::MultithreadInfo::n_threads();
  if (n_threads>1)
    stream << "--     . using " << n_threads << " threads " << (n_tasks == 1 ? "\n" : "each\n");

  stream << "-----------------------------------------------------------------------------\n"
         << std::endl;
}

template void print_aspect_header<std::ostream> (std::ostream &stream);
template void print_aspect_header<std::ofstream> (std::ofstream &stream);
