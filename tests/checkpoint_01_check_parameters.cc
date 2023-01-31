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

/*
 * Launch the following function when this plugin is created. Launch ASPECT
 * twice to test checkpoint/resume and then terminate the outer ASPECT run.
 */
int f()
{
  std::cout << "* starting from beginning:" << std::endl;

  // call ASPECT with "--" and pipe an existing input file into it.
  int ret;
  std::string command;

  command = ("cd output-checkpoint_01_check_parameters ; "
             "(cat " ASPECT_SOURCE_DIR "/tests/checkpoint_01_check_parameters.prm "
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

  command = ("cd output-checkpoint_01_check_parameters ; "
             " rm -rf output2.tmp ; mkdir output2.tmp ; "
             " cp output1.tmp/restart* output2.tmp/");
  std::cout << "Executing the following command:\n"
            << command
            << std::endl;
  ret = system (command.c_str());
  if (ret!=0)
    std::cout << "system() returned error " << ret << std::endl;


  std::cout << "* now resuming:" << std::endl;
  command = ("cd output-checkpoint_01_check_parameters ; "
             "(cat " ASPECT_SOURCE_DIR "/tests/checkpoint_01_check_parameters.prm "
             " ; "
             " echo 'set Output directory = output2.tmp' "
             " ; "
             " echo 'set Resume computation = true' "
             " ; "
             //
             // now also prescribe an incompatible number of compositional fields
             //
             " echo 'subsection Compositional fields' "
             " ; "
             " echo '  set Number of fields = 3' "
             " ; "
             " echo 'end' "
             " ; "
             ") "
             "| ../../aspect -- > /dev/null");
  std::cout << "Executing the following command:\n"
            << command
            << std::endl;

  // now run ASPECT on this (invalid) input file and output
  // the expected error code
  ret = system (command.c_str());
  if (ret!=0)
    std::cout << "system() returned error " << ret << std::endl;


  std::cout << "* now comparing:" << std::endl;

  command = ("cd output-checkpoint_01_check_parameters ; "
             "cp output1.tmp/log.txt log.txt1;"
             "cp output1.tmp/statistics statistics1;");
  ret = system (command.c_str());
  std::cout << "Executing the following command:\n"
            << command
            << std::endl;

  if (ret!=0)
    std::cout << "system() returned error " << ret << std::endl;

  // terminate current process:
  exit (0);
  return 42;
}


// run this function by initializing a global variable by it
int i = f();
