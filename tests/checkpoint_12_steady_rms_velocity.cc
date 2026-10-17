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
 * Run a steady-RMS-velocity calculation to completion, then resume from the
 * checkpoint made before it terminates. The steady-state time window spans
 * that checkpoint. Without restoring SteadyRMSVelocity::time_rmsvel, the
 * resumed calculation needs to rebuild its history and terminates later.
 */
int f()
{
  std::cout << "* running continuous reference:" << std::endl;

  run_command ("cd output-checkpoint_12_steady_rms_velocity ; "
               "rm -rf output_steady_rms_continuous.tmp ; "
               "(cat " ASPECT_SOURCE_DIR "/tests/checkpoint_12_steady_rms_velocity.prm "
               " ; "
               " echo 'set Output directory = output_steady_rms_continuous.tmp' "
               ") "
               "| ../../aspect -- > /dev/null");

  run_command ("cd output-checkpoint_12_steady_rms_velocity ; "
               "rm -rf output_steady_rms_restart.tmp ; mkdir output_steady_rms_restart.tmp ; "
               "cp -r output_steady_rms_continuous.tmp/restart output_steady_rms_restart.tmp/ ; "
               "echo 1 > output_steady_rms_restart.tmp/restart/last_good_checkpoint.txt");

  std::cout << "* now resuming from the pre-termination checkpoint:" << std::endl;

  run_command ("cd output-checkpoint_12_steady_rms_velocity ; "
               "(cat " ASPECT_SOURCE_DIR "/tests/checkpoint_12_steady_rms_velocity.prm "
               " ; "
               " echo 'set Output directory = output_steady_rms_restart.tmp' "
               " ; "
               " echo 'set Resume computation = true' "
               ") "
               "| ../../aspect -- > /dev/null");

  std::cout << "* restarted run log:" << std::endl;
  run_command ("cat output-checkpoint_12_steady_rms_velocity/"
               "output_steady_rms_restart.tmp/log.txt");

  std::cout << "* final continuous and restart statistics:" << std::endl;

  run_command ("cd output-checkpoint_12_steady_rms_velocity ; "
               "tail -n 1 output_steady_rms_continuous.tmp/statistics ; "
               "tail -n 1 output_steady_rms_restart.tmp/statistics");

  std::cout << "* checkpoint resume preserved steady-RMS-velocity termination history "
            << "if the last two lines are equal."
            << std::endl;

  std::exit(0);
  return 42;
}


// Run this function by initializing a global variable by it.
int i = f();
