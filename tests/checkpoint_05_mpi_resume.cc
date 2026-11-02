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
  // Wait for file system operations of the test checkpoint_05_mpi_create to finish
  std::this_thread::sleep_for(std::chrono::seconds(15));

  if (dealii::Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    {
      std::cout << "* Copying checkpoint files." << std::endl;

      std::string command;

      command = ("mkdir -p output-checkpoint_05_mpi_resume; "
                 "cp -r output-checkpoint_05_mpi_create/restart/ output-checkpoint_05_mpi_resume/");
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
    }

  // make sure the cp command above finished on all ranks:
  MPI_Barrier(MPI_COMM_WORLD);

  return 42;
}

// run this function by initializing a global variable by it
int i = f();
