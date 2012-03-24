#include <aspect/simulator.h>

#include <deal.II/base/utilities.h>


int main (int argc, char *argv[])
{
  using namespace dealii;
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);

  try
    {
      deallog.depth_console (0);

      // see which parameter file to use
      std::string parameter_filename;
      if (argc>=2)
        parameter_filename = argv[1];
      else
        parameter_filename = "box.prm";

      // declare parameters so that we can create a default file
      // if there is no parameter file
      ParameterHandler prm;
      aspect::Simulator<deal_II_dimension>::declare_parameters(prm);

      std::ifstream parameter_file (parameter_filename.c_str());
      if (!parameter_file)
        {
          parameter_file.close ();

          std::ostringstream message;
          message << "Input parameter file <"
                  << parameter_filename << "> not found. Creating a"
                  << std::endl
                  << "template file of the same name."
                  << std::endl;

          std::ofstream parameter_out (parameter_filename.c_str());
          prm.print_parameters (parameter_out,
                                ParameterHandler::Text);

          AssertThrow (false, ExcMessage (message.str().c_str()));
        }

      const bool success = prm.read_input (parameter_file);
      AssertThrow (success, ExcMessage ("Invalid input parameter file."));

      aspect::Simulator<deal_II_dimension> flow_problem (MPI_COMM_WORLD, prm);
      flow_problem.run ();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
