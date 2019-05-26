/*
  Copyright (C) 2018 by the authors of the World Builder code.

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

#include <algorithm>
#include <exception>
#include <iostream>
#include <fstream>

#include <world_builder/assert.h>
#include <world_builder/utilities.h>
#include <world_builder/world.h>

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

  std::string line;
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
      for (unsigned int i = 0; i < line.size(); ++i)
        line[i].erase(std::remove(line[i].begin(), line[i].end(), ','), line[i].end());

      data.push_back(line);
    }

  // Read config from data if pressent
  for (unsigned int i = 0; i < data.size(); ++i)
    {
      if (data[i].size() > 0 && data[i][0] == "#" && data[i][1] == "dim" && data[i][2] == "=")
        {
          dim = string_to_unsigned_int(data[i][3]);
        }

      if (data[i].size() > 0 && data[i][0] == "#" && data[i][1] == "compositions" && data[i][2] == "=")
        compositions = string_to_unsigned_int(data[i][3]);

    }

  switch (dim)
    {
      case 2:
        // set the header
        std::cout << "# x z d g T ";

        for (unsigned int c = 0; c < compositions; ++c)
          std::cout << "c" << c << " ";

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
              std::cout << std::endl;

            }
        break;
      case 3:
        // set the header
        std::cout << "# x y z d g T ";

        for (unsigned int c = 0; c < compositions; ++c)
          std::cout << "c" << c << " ";

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
              std::cout << std::endl;

            }
        break;
      default:
        std::cout << "The World Builder can only be run in 2d and 3d but a different space dimension " << std::endl
                  << "is given: dim = " << dim << ".";
        return 0;
    }

  return 0;
}
