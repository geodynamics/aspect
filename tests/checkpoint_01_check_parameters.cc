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
