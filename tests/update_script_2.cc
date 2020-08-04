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
 *
 * This test in particular calls the update script multiple times on the same
 * file, and ensures that only the first time changes happen. Every subsequent
 * application should not change the file any more.
 */
int f()
{
  int ret;
  std::string command;

  command = ("cp update_script_2.x.prm output-update_script_2/updated2.prm;"
             "sed -i.bak 's:set Additional shared libraries = ./libupdate_script_2.so::' output-update_script_2/updated2.prm;"
             "bash " ASPECT_SOURCE_DIR "/contrib/utilities/update_prm_files.sh output-update_script_2/updated2.prm;"
             "bash " ASPECT_SOURCE_DIR "/contrib/utilities/update_prm_files.sh output-update_script_2/updated2.prm;"
             "rm output-update_script_2/updated2.prm.bak");

  std::cout << "Executing the update script:\n"
            << command
            << std::endl;
  ret = system (command.c_str());
  if (ret!=0)
    std::cout << "system() returned error " << ret << std::endl;

  command = ("../aspect output-update_script_2/updated2.prm");
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
