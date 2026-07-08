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
#include <iostream>

/*
 * Launch the following function when this plugin is created. Run ASPECT twice
 * to test FastScape checkpoint/resume of the basement and silt fraction, then
 * terminate the outer ASPECT run.
 *
 * The first run goes to the End time while writing checkpoints. The second run
 * resumes from the last checkpoint and continues to the same End time. Because
 * FastScape is deterministic on the root process, the resumed run must
 * reproduce the uninterrupted run; comparing the two topography outputs and
 * statistics detects any corruption of the restarted basement/silt state.
 */
int f()
{
  std::cout << "* starting from beginning:" << std::endl;

  // call ASPECT with "--" and pipe an existing input file into it.
  int ret;
  std::string command;

  command = ("cd output-fastscape_restart_basement_silt ; "
             "(cat " ASPECT_SOURCE_DIR "/tests/fastscape_restart_basement_silt.prm "
             " ; "
             " echo 'set Output directory = output1.tmp' "
             " ; "
             " rm -rf output1.tmp ; mkdir output1.tmp "
             ") "
             "| ../../aspect -- > /dev/null");
  std::cout << "Executing the following command:\n"
            << command
            << std::endl;
  ret = system (command.c_str());
  if (ret!=0)
    std::cout << "system() returned error " << ret << std::endl;

  command = ("cd output-fastscape_restart_basement_silt ; "
             " rm -rf output2.tmp ; mkdir output2.tmp ; "
             " mv output1.tmp/restart* output2.tmp/");
  std::cout << "Executing the following command:\n"
            << command
            << std::endl;
  ret = system (command.c_str());
  if (ret!=0)
    std::cout << "system() returned error " << ret << std::endl;


  std::cout << "* now resuming:" << std::endl;
  command = ("cd output-fastscape_restart_basement_silt ; "
             "(cat " ASPECT_SOURCE_DIR "/tests/fastscape_restart_basement_silt.prm "
             " ; "
             " echo 'set Output directory = output2.tmp' "
             " ; "
             " echo 'set Resume computation = true' "
             ") "
             "| ../../aspect -- > /dev/null");
  std::cout << "Executing the following command:\n"
            << command
            << std::endl;
  ret = system (command.c_str());
  if (ret!=0)
    std::cout << "system() returned error " << ret << std::endl;

  std::cout << "* now comparing:" << std::endl;

  // Copy the final topography of both the uninterrupted (1) and the resumed (2)
  // run as references for the testsuite to diff against.
  ret = system ("cd output-fastscape_restart_basement_silt ; "
                "cp output1.tmp/topography/topography.00006 topography.000061;"
                "cp output2.tmp/topography/topography.00006 topography.000062;"
                "cp output1.tmp/statistics statistics1;"
                "cp output2.tmp/statistics statistics2;"
                "");
  if (ret!=0)
    std::cout << "system() returned error " << ret << std::endl;

  // Directly assert the restart invariant: a correct restart reproduces the
  // uninterrupted run, so the two final topographies must be identical. This
  // check is environment-independent and cannot be masked by regenerating the
  // reference output, because reintroducing the basement/silt save() bug makes
  // the resumed run diverge by meters.
  ret = system ("cd output-fastscape_restart_basement_silt ; "
                "if diff -q topography.000061 topography.000062 > /dev/null ; then "
                "  echo 'restart reproduces straight run: topography identical' ; "
                "else "
                "  echo 'ERROR: restart diverged from straight run' ; "
                "fi");
  if (ret!=0)
    std::cout << "system() returned error " << ret << std::endl;

  // terminate current process:
  exit (0);
  return 42;
}


// run this function by initializing a global variable by it
int i = f();
