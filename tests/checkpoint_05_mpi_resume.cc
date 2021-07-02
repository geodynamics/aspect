#include <aspect/simulator.h>
#include <iostream>

/*
 * Launch the following function when this plugin is created. Copy checkpoint
 * files into the correct place to resume model.
 */
int f()
{
  if (dealii::Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) != 0)
    return 42;

  std::cout << "* Copying checkpoint files." << std::endl;

  std::string command;

  command = ("mkdir -p output-checkpoint_05_mpi_resume; "
             "cp -r output-checkpoint_05_mpi_create/restart* output-checkpoint_05_mpi_resume/");
  std::cout << "Executing the following command:\n"
            << command
            << std::endl;
  const int ret = system (command.c_str());
  if (ret!=0)
    {
      std::cout << "system() returned error " << ret << std::endl;
      exit(1);
    }

  std::cout << "* Finished copying files. Now resuming model." << std::endl;

  return 42;
}

// run this function by initializing a global variable by it
int i = f();
