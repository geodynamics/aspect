// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2022-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

#include <aspect/global.h>

#include <deal.II/base/mpi.h>

#include <iostream>
#include <thread>
#include <chrono>
#include <cstdlib>
#include <string>

/*
 * Launch the following function when this plugin is created. Copy checkpoint
 * files into the correct place to resume model.
 */
int f()
{
  // Wait for file system operations of the test checkpoint_07_enable_free_surface_create to finish
  std::this_thread::sleep_for(std::chrono::seconds(15));

  if (dealii::Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) != 0)
    return 42;

  std::cout << "* Copying checkpoint files." << std::endl;

  std::string command;

  command = ("mkdir -p output-checkpoint_07_enable_free_surface_resume; "
             "mv output-checkpoint_07_enable_free_surface_create/restart output-checkpoint_07_enable_free_surface_resume/");
  std::cout << "Executing the following command:\n"
            << command
            << std::endl;
  const int ret = system (command.c_str());
  if (ret!=0)
    {
      std::cout << "system() returned error " << ret << std::endl;
      exit(1);
    }

  std::cout << "* Finished copying files. Now resuming model." << std::endl;

  return 42;
}

// run this function by initializing a global variable by it
int i = f();
