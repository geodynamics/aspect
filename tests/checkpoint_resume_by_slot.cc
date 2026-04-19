/*
  Copyright (C) 2022 - 2024 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.
*/

#include <aspect/simulator.h>

#include <cstdlib>
#include <iostream>
#include <string>

namespace
{
  int run_command(const std::string &command)
  {
    std::cout << "Executing the following command:\n"
              << command
              << std::endl;

    const int ret = system(command.c_str());
    if (ret != 0)
      std::cout << "system() returned error " << ret << std::endl;

    return ret;
  }
}

int f()
{
  std::cout << "* starting from beginning:" << std::endl;

  run_command("cd output-checkpoint_resume_by_slot ; "
              "(cat " ASPECT_SOURCE_DIR "/tests/checkpoint_resume_by_slot.prm "
              " ; "
              " echo 'set Output directory = output1.tmp' "
              " ; "
              " rm -rf output1.tmp ; mkdir output1.tmp "
              ") "
              "| ../../aspect -- > /dev/null");

  run_command("cd output-checkpoint_resume_by_slot ; "
              " rm -rf output2.tmp ; cp -r output1.tmp output2.tmp ; "
              " rm -f output1.tmp/log.txt");

  std::cout << "* now resuming from slot 3:" << std::endl;
  run_command("cd output-checkpoint_resume_by_slot ; "
              "(cat " ASPECT_SOURCE_DIR "/tests/checkpoint_resume_by_slot.prm "
              " ; "
              " echo 'set Output directory = output1.tmp' "
              " ; "
              " echo 'set Resume computation = true' "
              " ; "
              " echo 'subsection Checkpointing' "
              " ; "
              " echo '  set Resume checkpoint = 3' "
              " ; "
              " echo 'end' "
              ") "
              "| ../../aspect -- > /dev/null");

  std::cout << "* now comparing:" << std::endl;
  run_command("cd output-checkpoint_resume_by_slot ; "
              "diff -c output1.tmp/statistics output2.tmp/statistics ; "
              "cp output1.tmp/log.txt log.txt1 ; "
              "cp output2.tmp/log.txt log.txt2 ; "
              "cp output1.tmp/statistics statistics1 ; "
              "cp output2.tmp/statistics statistics2 ; "
              "grep 'Resuming from snapshot' output1.tmp/log.txt > resume-summary.txt");

  exit(0);
  return 42;
}

int i = f();
