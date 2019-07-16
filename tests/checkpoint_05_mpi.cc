#include <aspect/simulator.h>
#include <iostream>

/*
 * Launch the following function when this plugin is created. Launch ASPECT
 * twice to test checkpoint/resume and then terminate the outer ASPECT run.
 */
int f()
{
  if (dealii::Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) != 0)
    return 0;

  std::cout << "* starting from beginning:" << std::endl;

  // call ASPECT with "--" and pipe an existing input file into it.
  int ret;
  std::string command;

  command = ("cd output-checkpoint_05_mpi ; "
             "(cat " ASPECT_SOURCE_DIR "/tests/checkpoint_05_mpi.prm "
             " ; "
             " echo 'set Output directory = output1.tmp' "
             " ; "
             " rm -rf output1.tmp ; mkdir output1.tmp "
             ") "
             "| mpirun -np 2 ../../aspect -- > /dev/null");
  std::cout << "Executing the following command:\n"
            << command
            << std::endl;
  ret = system (command.c_str());
  if (ret!=0)
    std::cout << "system() returned error " << ret << std::endl;

  command = ("cd output-checkpoint_05_mpi ; "
             " rm -rf output2.tmp ; mkdir output2.tmp ; "
             " cp output1.tmp/restart* output2.tmp/");
  std::cout << "Executing the following command:\n"
            << command
            << std::endl;
  ret = system (command.c_str());
  if (ret!=0)
    std::cout << "system() returned error " << ret << std::endl;


  std::cout << "* now resuming:" << std::endl;
  command = ("cd output-checkpoint_05_mpi ; "
             "(cat " ASPECT_SOURCE_DIR "/tests/checkpoint_05_mpi.prm "
             " ; "
             " echo 'set Output directory = output2.tmp' "
             " ; "
             " echo 'set Resume computation = true' "
             ") "
             "| mpirun -np 2 ../../aspect -- > /dev/null");
  std::cout << "Executing the following command:\n"
            << command
            << std::endl;
  ret = system (command.c_str());
  if (ret!=0)
    std::cout << "system() returned error " << ret << std::endl;

  std::cout << "* now comparing:" << std::endl;

  ret = system ("cd output-checkpoint_05_mpi ; "
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
