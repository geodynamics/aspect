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
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/

#include <aspect/simulator.h>

#include <cstdlib>
#include <iostream>

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
 * Launch ASPECT twice when this plugin is created. The parameter file uses a
 * compositional field advected by particles, but does not enable the particles
 * postprocessor. The test succeeds if a checkpoint written by the first run can
 * be resumed by the second run.
 */
int f()
{
  std::cout << "* starting from beginning:" << std::endl;

  run_command ("cd output-checkpoint_09_particles_no_postprocessor ; "
               "(cat " ASPECT_SOURCE_DIR "/tests/checkpoint_09_particles_no_postprocessor.prm "
               " ; "
               " echo 'set Output directory = output1.tmp' "
               " ; "
               " rm -rf output1.tmp ; mkdir output1.tmp "
               ") "
               "| ../../aspect -- > /dev/null");

  run_command ("cd output-checkpoint_09_particles_no_postprocessor ; "
               "test ! -d output1.tmp/particles ; "
               "rm -rf output2.tmp ; mkdir output2.tmp ; "
               "mv output1.tmp/restart output2.tmp/");

  std::cout << "* now resuming:" << std::endl;

  run_command ("cd output-checkpoint_09_particles_no_postprocessor ; "
               "(cat " ASPECT_SOURCE_DIR "/tests/checkpoint_09_particles_no_postprocessor.prm "
               " ; "
               " echo 'set Output directory = output2.tmp' "
               " ; "
               " echo 'set Resume computation = true' "
               ") "
               "| ../../aspect -- > /dev/null");

  run_command ("cd output-checkpoint_09_particles_no_postprocessor ; "
               "test ! -d output2.tmp/particles");

  std::cout << "* checkpoint resume succeeded without particle output." << std::endl;

  std::exit (0);
  return 42;
}



// run this function by initializing a global variable by it
int i = f();
