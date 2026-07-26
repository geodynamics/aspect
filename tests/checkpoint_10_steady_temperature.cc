/*
  Copyright (C) 2026 by the authors of the ASPECT code.

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
  along with ASPECT; if not, write to the Free Software Foundation,
  Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#include <aspect/simulator.h>

#include <cstdlib>
#include <iostream>


namespace
{
  void
  run_command (const std::string &command)
  {
    const int ret = std::system(command.c_str());
    if (ret != 0)
      {
        std::cout << "system() returned error " << ret << std::endl;
        std::exit(1);
      }
  }


}


/*
 * Run a steady-temperature calculation to completion, then resume from the
 * checkpoint made before it terminates. The steady-state time window spans
 * that checkpoint. Without restoring SteadyTemperature::time_temperature,
 * the resumed calculation needs to rebuild its history and terminates later.
 */
int f()
{
  std::cout << "* running continuous reference:" << std::endl;

  run_command ("cd output-checkpoint_10_steady_temperature ; "
               "(cat " ASPECT_SOURCE_DIR "/tests/checkpoint_10_steady_temperature.prm "
               " ; "
               " echo 'set Output directory = output_continuous.tmp' "
               " ; "
               " rm -rf output_continuous.tmp "
               ") "
               "| ../../aspect -- > /dev/null");

  run_command ("cd output-checkpoint_10_steady_temperature ; "
               "rm -rf output_restart.tmp ; mkdir output_restart.tmp ; "
               "cp -r output_continuous.tmp/restart output_restart.tmp/ ; "
               "echo 1 > output_restart.tmp/restart/last_good_checkpoint.txt");

  std::cout << "* now resuming from the pre-termination checkpoint:" << std::endl;

  run_command ("cd output-checkpoint_10_steady_temperature ; "
               "(cat " ASPECT_SOURCE_DIR "/tests/checkpoint_10_steady_temperature.prm "
               " ; "
               " echo 'set Output directory = output_restart.tmp' "
               " ; "
               " echo 'set Resume computation = true' "
               ") "
               "| ../../aspect -- > /dev/null");

  std::cout << "* restarted run log:" << std::endl;
  run_command ("cat output-checkpoint_10_steady_temperature/"
               "output_restart.tmp/log.txt");

  std::cout << "* final continuous and restart statistics:" << std::endl;

  run_command ("cd output-checkpoint_10_steady_temperature ; "
               "tail -n 1 output_continuous.tmp/statistics ; "
               "tail -n 1 output_restart.tmp/statistics");

  std::cout << "* checkpoint resume preserved steady-temperature termination history "
            << "if the last two lines are equal."
            << std::endl;

  std::exit(0);
  return 42;
}


// Run this function by initializing a global variable by it.
int i = f();
