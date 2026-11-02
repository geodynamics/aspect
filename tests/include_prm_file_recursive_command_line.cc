// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2025-2026 by the ASPECT authors.
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
 * Launch the following function when this plugin is created. Launch ASPECT
 * from the command line and then continue with the outer ASPECT run.
 */
int f()
{
  std::cout << "* Starting from command line ..." << std::endl;

  // call ASPECT with "--" and pipe an existing input file into it.
  int ret;
  std::string command;

  command = ("cd output-include_prm_file_recursive_command_line ; "
             "(cat " ASPECT_SOURCE_DIR "/tests/include_prm_file_recursive_command_line.prm "
             " ; "
             " echo 'set Output directory = output1.tmp' "
             " ; "
             " rm -rf output1.tmp ; mkdir output1.tmp "
             ") "
             "| ../../aspect -- > /dev/null");
  std::cout << "* Executing the following command:\n"
            << command
            << std::endl;
  ret = system (command.c_str());
  if (ret!=0)
    std::cout << "system() returned error " << ret << std::endl;

  std::cout << "* Copying output files ..." << std::endl;

  ret = system ("cd output-include_prm_file_recursive_command_line ; "
                "cp output1.tmp/log.txt log_command_line.txt;"
                "cp output1.tmp/statistics statistics_command_line;"
                "");
  if (ret!=0)
    std::cout << "system() returned error " << ret << std::endl;

  std::cout << "* Continuing run from parameter file ..." << std::endl;

  return 42;
}


// run this function by initializing a global variable by it
int i = f();
