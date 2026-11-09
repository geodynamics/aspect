// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2022-2026 by the ASPECT authors.
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
 * twice to test checkpoint/resume and then terminate the outer ASPECT run.
 */
int f()
{
  std::cout << "* starting from beginning:" << std::endl;

  // call ASPECT with "--" and pipe an existing input file into it.
  int ret;
  std::string command;

  command = ("cd output-checkpoint_13_particle_generator_random ; "
             "rm -rf output1.tmp ; "
             "(cat " ASPECT_SOURCE_DIR "/tests/checkpoint_13_particle_generator_random.prm "
             " ; "
             " echo 'set Output directory = output1.tmp' "
             ") "
             "| env -i OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 PATH=\"$PATH\" mpirun --oversubscribe -np 2 ../../aspect -- > /dev/null");
  std::cout << "Executing the following command:\n"
            << command
            << std::endl;
  ret = system (command.c_str());
  if (ret!=0)
    std::cout << "system() returned error " << ret << std::endl;

  command = ("cd output-checkpoint_13_particle_generator_random ; "
             " rm -rf output2.tmp ; mkdir output2.tmp ; "
             " mv output1.tmp/restart output2.tmp/");
  std::cout << "Executing the following command:\n"
            << command
            << std::endl;
  ret = system (command.c_str());
  if (ret!=0)
    std::cout << "system() returned error " << ret << std::endl;


  std::cout << "* now resuming:" << std::endl;
  command = ("cd output-checkpoint_13_particle_generator_random ; "
             "(cat " ASPECT_SOURCE_DIR "/tests/checkpoint_13_particle_generator_random.prm "
             " ; "
             " echo 'set Output directory = output2.tmp' "
             " ; "
             " echo 'set Resume computation = true' "
             ") "
             "| env -i OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 PATH=\"$PATH\" mpirun --oversubscribe -np 2 ../../aspect -- > /dev/null");
  std::cout << "Executing the following command:\n"
            << command
            << std::endl;
  ret = system (command.c_str());
  if (ret!=0)
    std::cout << "system() returned error " << ret << std::endl;

  std::cout << "* now comparing:" << std::endl;

  ret = system ("cd output-checkpoint_13_particle_generator_random ; "
                "cp output1.tmp/particles/particles-00006.0000.gnuplot particles-00006.0000.gnuplot1;"
                "cp output2.tmp/particles/particles-00006.0000.gnuplot particles-00006.0000.gnuplot2;"
                "cp output1.tmp/log.txt log.txt1;"
                "cp output2.tmp/log.txt log.txt2;"
                "cp output1.tmp/statistics statistics1;"
                "cp output2.tmp/statistics statistics2;"
                "");
  if (ret!=0)
    std::cout << "system() returned error " << ret << std::endl;

  // terminate current process:
  exit (0);
  return 42;
}


// run this function by initializing a global variable by it
int i = f();
