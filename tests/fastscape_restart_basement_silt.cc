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
                "mkdir -p topography;"
                "cp output1.tmp/topography/topography.00006 topography/topography.000061;"
                "cp output2.tmp/topography/topography.00006 topography/topography.000062;"
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
                "if diff -q topography/topography.000061 topography/topography.000062 > /dev/null ; then "
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
