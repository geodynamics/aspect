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
/*  $Id$  */


#include <aspect/simulator.h>

#include <deal.II/base/utilities.h>
#include <deal.II/base/mpi.h>


int main (int argc, char *argv[])
{
  using namespace dealii;
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);

  try
    {
      deallog.depth_console(0);

      // see which parameter file to use
      std::string parameter_filename;
      if (argc >= 2)
        parameter_filename = argv[1];
      else
        parameter_filename = "box.prm";

      std::ifstream parameter_file(parameter_filename.c_str());
      if (!parameter_file)
        {
          parameter_file.close();

          const std::string message = (std::string("Input parameter file <")
              + parameter_filename + "> not found.");
          AssertThrow(false, ExcMessage (message));
        }

      // now that we know that the file can, at least in principle, be read
      // try to determine the dimension we want to work in. the default
      // is 2, but if we find a line of the kind "set Dimension = ..."
      // then the last such line wins
      unsigned int dim = 2;
        {
          std::ifstream x_file(parameter_filename.c_str());
          while (x_file)
            {
              // get one line and strip spaces at the front and back
              std::string line;
              std::getline(x_file, line);
              while ((line.size() > 0) && (line[0] == ' ' || line[0] == '\t'))
                line.erase(0, 1);
              while ((line.size() > 0)
                  && (line[line.size() - 1] == ' '
                      || line[line.size() - 1] == '\t'))
                line.erase(line.size() - 1, std::string::npos);

              // now see whether the line starts with 'set' followed by multiple spaces
              // if now, try next line
              if (line.size() < 4)
                continue;
              if ((line[0] != 's') || (line[1] != 'e') || (line[2] != 't')
                  || !(line[3] == ' ' || line[3] == '\t'))
                continue;

              line.erase(0, 4);
              while ((line.size() > 0) && (line[0] == ' ' || line[0] == '\t'))
                line.erase(0, 1);

              // now see whether the next word is "Dimension"
              if (line.size() < 4)
                continue;
              if (line.find("Dimension") != 0)
                continue;
              line.erase(0, 9);
              while ((line.size() > 0) && (line[0] == ' ' || line[0] == '\t'))
                line.erase(0, 1);

              // we'd expect an equals size here
              if ((line.size() < 1) || (line[0] != '='))
                continue;
              line.erase(0, 1);
              while ((line.size() > 0) && (line[0] == ' ' || line[0] == '\t'))
                line.erase(0, 1);

              // the rest should now be an integer
              dim = Utilities::string_to_int(line);
            }
        }

      // now switch between the templates that code for 2d or 3d. it
      // would be nicer if we didn't have to duplicate code, but the
      // following needs to be known at compile time whereas the dimensionality
      // is only read at run-time
      ParameterHandler prm;

      switch (dim)
        {
      case 2:
        {
          aspect::Simulator<2>::declare_parameters(prm);

          const bool success = prm.read_input(parameter_file);
          AssertThrow(success, ExcMessage ("Invalid input parameter file."));

          aspect::Simulator<2> flow_problem(MPI_COMM_WORLD, prm);
          flow_problem.run();

          break;
        }

      case 3:
        {
          aspect::Simulator<3>::declare_parameters(prm);

          const bool success = prm.read_input(parameter_file);
          AssertThrow(success, ExcMessage ("Invalid input parameter file."));

          aspect::Simulator<3> flow_problem(MPI_COMM_WORLD, prm);
          flow_problem.run();

          break;
        }

      default:
        AssertThrow((dim >= 2) && (dim <= 3),
                    ExcMessage ("ASPECT can only be run in 2d and 3d."));
        }
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
