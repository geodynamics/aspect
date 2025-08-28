/*
  Copyright (C) 2025 by the authors of the ASPECT code.

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
