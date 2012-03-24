/*
  Copyright (C) 2011, 2012 by the authors of the ASPECT code.

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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/


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
