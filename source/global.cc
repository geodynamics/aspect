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


#include <aspect/global.h>
#include <aspect/revision.h>

#include <deal.II/base/multithread_info.h>
#include <deal.II/base/revision.h>
#include <deal.II/base/vectorization.h>

#ifdef ASPECT_WITH_WORLD_BUILDER
#include <world_builder/config.h>
#endif

#include <cstring>
#include <fstream>

#include <mpi.h>


namespace aspect
{
  template <class Stream>
  void print_aspect_header(Stream &stream)
  {
    const int n_tasks = dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

    stream << "-----------------------------------------------------------------------------\n"
           << "--                             This is ASPECT                              --\n"
           << "-- The Advanced Solver for Planetary Evolution, Convection, and Tectonics. --\n"
           << "-----------------------------------------------------------------------------\n"
           << "--     . version " << ASPECT_PACKAGE_VERSION;
    if (std::strcmp(ASPECT_GIT_BRANCH,"") != 0)
      stream << " (" << ASPECT_GIT_BRANCH << ", " << ASPECT_GIT_SHORTREV << ")\n";
    else
      stream << "\n";

    stream << "--     . using deal.II " << DEAL_II_PACKAGE_VERSION;
    if (std::strcmp(DEAL_II_GIT_BRANCH,"") != 0)
      stream << " (" << DEAL_II_GIT_BRANCH << ", " << DEAL_II_GIT_SHORTREV << ')';
    stream << "\n";
    stream << "--     .       with "
#ifdef DEAL_II_WITH_64BIT_INDICES
           << "64"
#else
           << "32"
#endif
           << " bit indices\n";

    const unsigned int n_vectorization_doubles = dealii::VectorizedArray<double>::size();
    const unsigned int n_vectorization_bits = 8 * sizeof(double) * n_vectorization_doubles;
    const std::string vectorization_level = dealii::Utilities::System::get_current_vectorization_level();

    stream << "--     .       with vectorization level ";
    stream << DEAL_II_COMPILER_VECTORIZATION_LEVEL
           << " (" << vectorization_level << ", "
           << n_vectorization_doubles << " doubles, "
           << n_vectorization_bits << " bits)\n";

    stream << "--     . using Trilinos "
           << DEAL_II_TRILINOS_VERSION_MAJOR    << '.'
           << DEAL_II_TRILINOS_VERSION_MINOR    << '.'
           << DEAL_II_TRILINOS_VERSION_SUBMINOR << '\n';
    stream << "--     . using p4est "
           << DEAL_II_P4EST_VERSION_MAJOR << '.'
           << DEAL_II_P4EST_VERSION_MINOR << '.'
           << DEAL_II_P4EST_VERSION_SUBMINOR << '\n';

#ifdef ASPECT_WITH_WORLD_BUILDER
    stream << "--     . using Geodynamic World Builder "
           << WORLD_BUILDER_VERSION_MAJOR << '.'
           << WORLD_BUILDER_VERSION_MINOR << '.'
           << WORLD_BUILDER_VERSION_PATCH;

    if (WorldBuilder::Version::GIT_SHA1 != "")
      stream << " (" << WorldBuilder::Version::GIT_BRANCH << ", " << WorldBuilder::Version::GIT_SHA1.substr(0,9) << ')';

    stream << '\n';
#endif

#ifdef ASPECT_WITH_FASTSCAPE
    stream << "--     . using FastScape "
           << FASTSCAPE_VERSION_MAJOR << '.'
           << FASTSCAPE_VERSION_MINOR << '.'
           << FASTSCAPE_VERSION_SUBMINOR << '\n';
#endif

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

}
