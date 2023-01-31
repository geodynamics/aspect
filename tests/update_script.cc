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
 * Launch the following function when this plugin is created. Use the update
 * script to update this parameter file. Unfortunately at this point the file
 * has already been loaded, so the current ASPECT instance would get the old
 * file. Thus we start a new ASPECT instance that will load the new file.
 * To avoid an endless recursion we remove the shared library from the new
 * input file (otherwise the new instance would call this library again, and so
 * on). After finishing the new instance we exit to not continue the old one.
 */
int f()
{
  int ret;
  std::string command;

  command = ("cp update_script.x.prm output-update_script/updated.prm &&"
             "sed -i.bak 's:set Additional shared libraries = ./libupdate_script.so::' output-update_script/updated.prm &&"
             "bash " ASPECT_SOURCE_DIR "/contrib/utilities/update_prm_files.sh output-update_script/updated.prm &&"
             "rm output-update_script/updated.prm.bak");

  std::cout << "Executing the update script:\n"
            << command
            << std::endl;
  ret = system (command.c_str());
  if (ret!=0)
    std::cout << "system() returned error " << ret << std::endl;

  command = ("../aspect output-update_script/updated.prm");
  std::cout << "Running ASPECT with updated parameter file:\n"
            << command
            << std::endl;
  ret = system (command.c_str());
  if (ret!=0)
    std::cout << "system() returned error " << ret << std::endl;

  // abort current process:
  exit (0);
  return 42;
}


// run this function by initializing a global variable by it
int i = f();
