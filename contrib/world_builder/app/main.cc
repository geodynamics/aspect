/*
  Copyright (C) 2018 - 2021 by the authors of the World Builder code.

  This file is part of the World Builder.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published
   by the Free Software Foundation, either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#ifdef WB_WITH_MPI
#include <mpi.h>
#endif

#include <algorithm>
#include <exception>
#include <iostream>
#include <fstream>

#include <world_builder/assert.h>
#include <world_builder/utilities.h>
#include <world_builder/world.h>

#include <app/main.h>

using namespace WorldBuilder::Utilities;

bool find_command_line_option(char **begin, char **end, const std::string &option)
{
  return std::find(begin, end, option) != end;
}

int main(int argc, char **argv)
{
  /**
   * First parse the command line options
   */
  std::string wb_file;
  std::string data_file;

  unsigned int dim = 3;
  unsigned int compositions = 0;
  unsigned int grain_compositions = 0;
  size_t number_of_grains = 0;

  if (find_command_line_option(argv, argv+argc, "-h") || find_command_line_option(argv, argv+argc, "--help"))
    {
      std::cout << "This program allows to use the world builder library directly with a world builder file and a data file. "
                "The data file will be filled with intitial conditions from the world as set by the world builder file." << std::endl
                << "Besides providing two files, where the first is the world builder file and the second is the data file, the available options are: " << std::endl
                << "-h or --help to get this help screen." << std::endl;
      return 0;
    }

  if (argc == 1)
    {
      std::cout << "Error: There where no files passed to the World Builder, use --help for more " << std::endl
                << "information on how  to use the World Builder app." << std::endl;
      return 0;
    }


  if (argc == 2)
    {
      std::cout << "Error:  The World Builder app requires at least two files, a World Builder file " << std::endl
                << "and a data file to convert." << std::endl;
      return 0;
    }

  if (argc != 3)
    {
      std::cout << "Only two command line arguments may be given, which should be the world builder file location and the data file location (in that order). " << std::endl;
      return 0;
    }

  int MPI_RANK = 0;
  int MPI_SIZE = 1;
#ifdef WB_WITH_MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &MPI_RANK);
  MPI_Comm_size(MPI_COMM_WORLD, &MPI_SIZE);
#endif

  if (MPI_RANK == 0)
    {
      wb_file = argv[1];
      data_file = argv[2];

      /**
       * Try to start the world builder
       */
      std::unique_ptr<WorldBuilder::World> world;
      //try
      {
        std::string output_dir = wb_file.substr(0,wb_file.find_last_of("/\\") + 1);
        world = std::unique_ptr<WorldBuilder::World>(new WorldBuilder::World(wb_file, true, output_dir));
      }
      /*catch (std::exception &e)
        {
          std::cerr << "Could not start the World builder, error: " << e.what() << "\n";
          return 1;
        }
      catch (...)
        {
          std::cerr << "Exception of unknown type!\n";
          return 1;
        }*/


      /**
       * Read the data from the data files
       */
      std::ifstream data_stream(data_file);

      // move the data into a vector of strings
      std::vector<std::vector<std::string> > data;
      std::string temp;

      while (std::getline(data_stream, temp))
        {
          std::istringstream buffer(temp);
          std::vector<std::string> line((std::istream_iterator<std::string>(buffer)),
                                        std::istream_iterator<std::string>());

          // remove the comma's in case it is a comma separated file.
          // TODO: make it split for comma's and/or spaces
          for (auto &line_i : line)
            line_i.erase(std::remove(line_i.begin(), line_i.end(), ','), line_i.end());

          data.push_back(line);
        }

      // Read config from data if pressent
      for (auto &line_i : data)
        {
          if (!line_i.empty() && line_i[0] == "#" && line_i[1] == "dim" && line_i[2] == "=")
            {
              dim = string_to_unsigned_int(line_i[3]);
            }

          if (!line_i.empty() && line_i[0] == "#" && line_i[1] == "compositions" && line_i[2] == "=")
            compositions = string_to_unsigned_int(line_i[3]);

          if (!line_i.empty() && line_i[0] == "#" && line_i[1] == "grain" && line_i[2] == "compositions" && line_i[3] == "=")
            grain_compositions = string_to_unsigned_int(line_i[4]);

          if (!line_i.empty() && line_i[0] == "#" && line_i[1] == "number" && line_i[2] == "of" && line_i[3] == "grains" && line_i[4] == "=")
            number_of_grains = string_to_unsigned_int(line_i[5]);

        }

      switch (dim)
        {
          case 2:
            // set the header
            std::cout << "# x z d g T ";

            for (unsigned int c = 0; c < compositions; ++c)
              std::cout << "c" << c << " ";

            for (unsigned int gc = 0; gc < grain_compositions; ++gc)
              for (size_t g = 0; g < number_of_grains; g++)
                std::cout << "gs" << gc << "-" << g << " " // gs = grain size, gm = grain rotation matrix
                          << "gm" << gc << "-" << g << "[0:0] " << "gm" << gc << "-" << g << "[0:1] " << "gm" << gc << "-" << g << "[0:2] "
                          << "gm" << gc << "-" << g << "[1:0] " << "gm" << gc << "-" << g << "[1:1] " << "gm" << gc << "-" << g << "[1:2] "
                          << "gm" << gc << "-" << g << "[2:0] " << "gm" << gc << "-" << g << "[2:1] " << "gm" << gc << "-" << g << "[2:2] ";

            std::cout <<std::endl;

            // set the values
            for (unsigned int i = 0; i < data.size(); ++i)
              if (data[i][0] != "#")
                {

                  WBAssertThrow(data[i].size() == dim + 2, "The file needs to contain dim + 2 entries, but contains " << data[i].size() << " entries "
                                " on line " << i+1 << " of the data file.  Dim is " << dim << ".");

                  std::array<double,2> coords = {{
                      string_to_double(data[i][0]),
                      string_to_double(data[i][1])
                    }
                  };
                  std::cout << data[i][0] << " " << data[i][1] << " " << data[i][2] << " " << data[i][3] << " ";
                  std::cout << world->temperature(coords, string_to_double(data[i][2]), string_to_double(data[i][3]))  << " ";

                  for (unsigned int c = 0; c < compositions; ++c)
                    {
                      std::cout << world->composition(coords, string_to_double(data[i][2]), c)  << " ";
                    }

                  for (unsigned int gc = 0; gc < grain_compositions; ++gc)
                    {
                      WorldBuilder::grains grains = world->grains(coords, string_to_double(data[i][2]), gc, number_of_grains);
                      for (unsigned int g = 0; g < number_of_grains; ++g)
                        {
                          std::cout << grains.sizes[g]  << " "
                                    << grains.rotation_matrices[g][0][0] << " " << grains.rotation_matrices[g][0][1] << " " << grains.rotation_matrices[g][0][2] << " "
                                    << grains.rotation_matrices[g][1][0] << " " << grains.rotation_matrices[g][1][1] << " " << grains.rotation_matrices[g][1][2] << " "
                                    << grains.rotation_matrices[g][2][0] << " " << grains.rotation_matrices[g][2][1] << " " << grains.rotation_matrices[g][2][2] << " ";
                        }
                    }
                  std::cout << std::endl;

                }
            break;
          case 3:
            // set the header
            std::cout << "# x y z d g T ";

            for (unsigned int c = 0; c < compositions; ++c)
              std::cout << "c" << c << " ";

            for (unsigned int gc = 0; gc < grain_compositions; ++gc)
              for (size_t g = 0; g < number_of_grains; g++)
                std::cout << "gs" << gc << "-" << g << " " // gs = grain size, gm = grain rotation matrix
                          << "gm" << gc << "-" << g << "[0:0] " << "gm" << gc << "-" << g << "[0:1] " << "gm" << gc << "-" << g << "[0:2] "
                          << "gm" << gc << "-" << g << "[1:0] " << "gm" << gc << "-" << g << "[1:1] " << "gm" << gc << "-" << g << "[1:2] "
                          << "gm" << gc << "-" << g << "[2:0] " << "gm" << gc << "-" << g << "[2:1] " << "gm" << gc << "-" << g << "[2:2] ";

            std::cout <<std::endl;

            // set the values
            for (unsigned int i = 0; i < data.size(); ++i)
              if (data[i][0] != "#")
                {
                  WBAssertThrow(data[i].size() == dim + 2, "The file needs to contain dim + 2 entries, but contains " << data[i].size() << " entries "
                                " on line " << i+1 << " of the data file. Dim is " << dim << ".");
                  std::array<double,3> coords = {{
                      string_to_double(data[i][0]),
                      string_to_double(data[i][1]),
                      string_to_double(data[i][2])
                    }
                  };


                  std::cout << data[i][0] << " " << data[i][1] << " " << data[i][2] << " " << data[i][3] << " " << data[i][4] << " ";
                  std::cout << world->temperature(coords, string_to_double(data[i][3]), string_to_double(data[i][4]))  << " ";

                  for (unsigned int c = 0; c < compositions; ++c)
                    {
                      std::cout << world->composition(coords, string_to_double(data[i][3]), c)  << " ";
                    }

                  for (unsigned int gc = 0; gc < grain_compositions; ++gc)
                    {
                      WorldBuilder::grains grains = world->grains(coords, string_to_double(data[i][3]), gc, number_of_grains);
                      for (unsigned int g = 0; g < number_of_grains; ++g)
                        {
                          std::cout << grains.sizes[g]  << " "
                                    << grains.rotation_matrices[g][0][0] << " " << grains.rotation_matrices[g][0][1] << " " << grains.rotation_matrices[g][0][2] << " "
                                    << grains.rotation_matrices[g][1][0] << " " << grains.rotation_matrices[g][1][1] << " " << grains.rotation_matrices[g][1][2] << " "
                                    << grains.rotation_matrices[g][2][0] << " " << grains.rotation_matrices[g][2][1] << " " << grains.rotation_matrices[g][2][2] << " ";
                        }
                    }
                  std::cout << std::endl;

                }
            break;
          default:
            std::cout << "The World Builder can only be run in 2d and 3d but a different space dimension " << std::endl
                      << "is given: dim = " << dim << ".";

#ifdef WB_WITH_MPI
            MPI_Finalize();
#endif
            return 0;
        }
    }
#ifdef WB_WITH_MPI
  MPI_Finalize();
#endif
  return 0;
}
