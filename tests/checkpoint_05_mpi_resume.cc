/*
  Copyright (C) 2022 by the authors of the ASPECT code.

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

#include <aspect/simulator.h>

#include <iostream>
#include <thread>
#include <chrono>

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
