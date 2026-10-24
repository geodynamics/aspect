// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------



#include <aspect/global.h>

#include <cstdlib>
#include <iostream>
#include <string>


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
 * Run a steady-heat-flux calculation to completion, then resume from the
 * checkpoint made before it terminates. The steady-state time window spans
 * that checkpoint. Without restoring SteadyHeatFlux::time_heat_flux, the
 * resumed calculation needs to rebuild its history and terminates later.
 */
int f()
{
  std::cout << "* running continuous reference:" << std::endl;

  run_command ("cd output-checkpoint_11_steady_heat_flux ; "
               "(cat " ASPECT_SOURCE_DIR "/tests/checkpoint_11_steady_heat_flux.prm "
               " ; "
               " echo 'set Output directory = output_continuous.tmp' "
               " ; "
               " rm -rf output_continuous.tmp "
               ") "
               "| ../../aspect -- > /dev/null");

  run_command ("cd output-checkpoint_11_steady_heat_flux ; "
               "rm -rf output_restart.tmp ; mkdir output_restart.tmp ; "
               "cp -r output_continuous.tmp/restart output_restart.tmp/ ; "
               "echo 1 > output_restart.tmp/restart/last_good_checkpoint.txt");

  std::cout << "* now resuming from the pre-termination checkpoint:" << std::endl;

  run_command ("cd output-checkpoint_11_steady_heat_flux ; "
               "(cat " ASPECT_SOURCE_DIR "/tests/checkpoint_11_steady_heat_flux.prm "
               " ; "
               " echo 'set Output directory = output_restart.tmp' "
               " ; "
               " echo 'set Resume computation = true' "
               ") "
               "| ../../aspect -- > /dev/null");

  std::cout << "* restarted run log:" << std::endl;
  run_command ("cat output-checkpoint_11_steady_heat_flux/"
               "output_restart.tmp/log.txt");

  std::cout << "* final continuous and restart statistics:" << std::endl;

  run_command ("cd output-checkpoint_11_steady_heat_flux ; "
               "tail -n 1 output_continuous.tmp/statistics ; "
               "tail -n 1 output_restart.tmp/statistics");

  std::cout << "* checkpoint resume preserved steady-heat-flux termination history "
            << "if the last two lines are equal."
            << std::endl;

  std::exit(0);
  return 42;
}


// Run this function by initializing a global variable by it.
int i = f();
