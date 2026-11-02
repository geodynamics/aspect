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
    const int ret = std::system (command.c_str());
    if (ret != 0)
      {
        std::cout << "system() returned error " << ret << std::endl;
        std::exit (1);
      }
  }
}



/*
 * Launch ASPECT three times when this plugin is created. The parameter file
 * uses a compositional field advected by particles, but does not enable the
 * particles postprocessor. The test compares a continuous run with a run that
 * is resumed from a checkpoint to make sure particle state is preserved even
 * without particle output.
 */
int f()
{
  std::cout << "* running continuous reference:" << std::endl;

  run_command ("cd output-checkpoint_09_particles_no_postprocessor ; "
               "(cat " ASPECT_SOURCE_DIR "/tests/checkpoint_09_particles_no_postprocessor.prm "
               " ; "
               " echo 'set Output directory = output_continuous.tmp' "
               " ; "
               " rm -rf output_continuous.tmp "
               ") "
               "| ../../aspect -- > /dev/null");

  run_command ("cd output-checkpoint_09_particles_no_postprocessor ; "
               "test ! -d output_continuous.tmp/particles");

  std::cout << "* running restart setup:" << std::endl;

  run_command ("cd output-checkpoint_09_particles_no_postprocessor ; "
               "(cat " ASPECT_SOURCE_DIR "/tests/checkpoint_09_particles_no_postprocessor.prm "
               " ; "
               " echo 'set Output directory = output_restart_start.tmp' "
               " ; "
               " echo 'subsection Termination criteria' "
               " ; "
               " echo '  set End step = 3' "
               " ; "
               " echo 'end' "
               " ; "
               " rm -rf output_restart_start.tmp "
               ") "
               "| ../../aspect -- > /dev/null");

  run_command ("cd output-checkpoint_09_particles_no_postprocessor ; "
               "test ! -d output_restart_start.tmp/particles ; "
               "rm -rf output_restart_resume.tmp ; mkdir output_restart_resume.tmp ; "
               "cp -r output_restart_start.tmp/restart output_restart_resume.tmp/");

  std::cout << "* now resuming:" << std::endl;

  run_command ("cd output-checkpoint_09_particles_no_postprocessor ; "
               "(cat " ASPECT_SOURCE_DIR "/tests/checkpoint_09_particles_no_postprocessor.prm "
               " ; "
               " echo 'set Output directory = output_restart_resume.tmp' "
               " ; "
               " echo 'set Resume computation = true' "
               ") "
               "| ../../aspect -- > /dev/null");

  run_command ("cd output-checkpoint_09_particles_no_postprocessor ; "
               "test ! -d output_restart_resume.tmp/particles");

  std::cout << "* comparing restarted run with continuous run:" << std::endl;

  run_command ("cd output-checkpoint_09_particles_no_postprocessor ; "
               "tail -n 1 output_continuous.tmp/statistics > continuous.last ; "
               "tail -n 1 output_restart_resume.tmp/statistics > restart.last ; "
               "diff -u continuous.last restart.last");

  std::cout << "* checkpoint resume preserved particle state without particle output." << std::endl;

  std::exit (0);
  return 42;
}



// run this function by initializing a global variable by it
int i = f();
