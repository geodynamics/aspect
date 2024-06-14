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


#include "app/main.h"

#include "world_builder/assert.h"
#include "world_builder/consts.h"
#include "world_builder/point.h"
#include "world_builder/utilities.h"
#include "world_builder/world.h"

#ifdef WB_WITH_MPI
// we don't need the c++ MPI wrappers
#define OMPI_SKIP_MPICXX 1
#define MPICH_SKIP_MPICXX
#include <mpi.h>
#endif

#include <algorithm>
#include <fstream>
#include <iostream>
#include <memory>

#ifndef NDEBUG
#ifdef WB_USE_FP_EXCEPTIONS
#include <cfenv>
#endif
#endif

using namespace WorldBuilder::Utilities;

bool find_command_line_option(char **begin, char **end, const std::string &option)
{
  return std::find(begin, end, option) != end;
}

int main(int argc, char **argv)
{

#ifndef NDEBUG
#ifdef WB_USE_FP_EXCEPTIONS
  // Some implementations seem to not initialize the floating point exception
  // bits to zero. Make sure we start from a clean state.
  feclearexcept(FE_DIVBYZERO|FE_INVALID);

  // enable floating point exceptions
  feenableexcept(FE_DIVBYZERO|FE_INVALID);
#endif
#endif

  /**
   * First parse the command line options
   */
  std::string wb_file;
  std::string data_file;

  unsigned int dim = 3;
  unsigned int compositions = 0;
  unsigned int grain_compositions = 0;
  size_t n_grains = 0;
  bool convert_spherical = false;
  bool limit_debug_consistency_checks = false;
  bool output_json_files = false;

  if (find_command_line_option(argv, argv+argc, "-h") || find_command_line_option(argv, argv+argc, "--help"))
    {
      std::cout << "This program allows to use the world builder library directly with a world builder file and a data file. "
                "The data file will be filled with initial conditions from the world as set by the world builder file." << std::endl
                << "Besides providing two files, where the first is the world builder file and the second is the data file, the available options are: " << std::endl
                << "-h or --help to get this help screen." << std::endl;
      return 0;
    }

  if (find_command_line_option(argv, argv+argc, "-ldcc") || find_command_line_option(argv, argv+argc, "--limit-debug-consistency-checks"))
    limit_debug_consistency_checks = true;

  if (find_command_line_option(argv, argv+argc, "--output-json-files"))
    output_json_files = true;

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

  if ((argc == 3 && limit_debug_consistency_checks && output_json_files) || (argc == 4 && !(!limit_debug_consistency_checks != !output_json_files)) || (argc == 5 && (!limit_debug_consistency_checks && !output_json_files)) || argc > 5)
    {
      std::cout << "Only exactly two command line arguments may be given, which should be the world builder file location and the data file location (in that order) "
                << "or exactly three command line arguments, which should be the world builder file location, the data file location and --limit-debug-consistency-checks or --output-json-files (in that order),"
                "or exactly four command line arguments, which should be the world builder file location, the data file location and --limit-debug-consistency-checks and --output-json-files (in that order),"
                << ", argc = " << argc << ", limit_debug_consistency_checks = " << (limit_debug_consistency_checks ? "true" : "false") << ", output_json_files = " << (output_json_files ? "true" : "false") << std::endl;
      return 0;
    }

  int MPI_RANK = 0;
#ifdef WB_WITH_MPI
  int MPI_SIZE = 1;
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
        const std::string output_dir = wb_file.substr(0,wb_file.find_last_of("/\\") + 1);
        world = std::make_unique<WorldBuilder::World>(wb_file, output_json_files, output_dir,1,limit_debug_consistency_checks);
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

      // Read config from data if present
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
            n_grains = string_to_unsigned_int(line_i[5]);

          if (!line_i.empty() && line_i[0] == "#" && line_i[1] == "convert" && line_i[2] == "spherical" && line_i[3] == "=" && line_i[4] == "true")
            convert_spherical = true;
        }

      // set properties
      std::vector<std::array<unsigned ,3>> properties;
      properties.push_back({{1,0,0}}); // temperature

      for (size_t c = 0; c < compositions; ++c)
        properties.push_back({{2,static_cast<unsigned int>(c),0}}); // composition c

      for (size_t gc = 0; gc < grain_compositions; ++gc)
        properties.push_back({{3,static_cast<unsigned int>(gc),static_cast<unsigned int>(n_grains)}}); // grains gc

      properties.push_back({{4,0,0}}); // tag


      switch (dim)
        {
          case 2:
            WBAssertThrow(!convert_spherical, "Converting to spherical values is only available in 3D.");
            // set the header
            std::cout << "# x z d T ";

            for (unsigned int c = 0; c < compositions; ++c)
              std::cout << 'c' << c << ' ';

            for (unsigned int gc = 0; gc < grain_compositions; ++gc)
              for (size_t g = 0; g < n_grains; g++)
                std::cout << "gs" << gc << '-' << g << ' ' // gs = grain size, gm = grain rotation matrix
                          << "gm" << gc << '-' << g << "[0:0] " << "gm" << gc << '-' << g << "[0:1] " << "gm" << gc << '-' << g << "[0:2] "
                          << "gm" << gc << '-' << g << "[1:0] " << "gm" << gc << '-' << g << "[1:1] " << "gm" << gc << '-' << g << "[1:2] "
                          << "gm" << gc << '-' << g << "[2:0] " << "gm" << gc << '-' << g << "[2:1] " << "gm" << gc << '-' << g << "[2:2] ";

            std::cout << "tag ";
            std::cout <<std::endl;

            // set the values
            for (unsigned int i = 0; i < data.size(); ++i)
              if (data[i].size() > 0 && data[i][0] != "#")
                {

                  WBAssertThrow(data[i].size() == dim + 1, "The file needs to contain dim + 1 entries, but contains " << data[i].size() << " entries "
                                " on line " << i+1 << " of the data file (" << data_file << ").  Dim is " << dim << '.');
                  const std::array<double,2> coords = {{
                      string_to_double(data[i][0]),
                      string_to_double(data[i][1])
                    }
                  };
                  std::cout << data[i][0] << ' ' << data[i][1] << ' ' << data[i][2] << ' ';
                  std::vector<double> output = world->properties(coords, string_to_double(data[i][2]),properties);
                  std::cout << output[0]  << ' ';

                  for (unsigned int c = 0; c < compositions; ++c)
                    {
                      std::cout << output[1+c]  << ' ';
                    }

                  for (unsigned int gc = 0; gc < grain_compositions; ++gc)
                    {
                      const size_t start = 1+compositions+gc*n_grains*10;
                      for (unsigned int g = 0; g < n_grains; ++g)
                        {
                          std::cout << output[start+g]  << ' '
                                    << output[start+n_grains+g*9] << ' ' << output[start+n_grains+g*9+1] << ' ' << output[start+n_grains+g*9+2] << ' '
                                    << output[start+n_grains+g*9+3] << ' ' << output[start+n_grains+g*9+4] << ' ' << output[start+n_grains+g*9+5] << ' '
                                    << output[start+n_grains+g*9+6] << ' ' << output[start+n_grains+g*9+7] << ' ' << output[start+n_grains+g*9+8] << ' ';

                        }
                    }
                  std::cout << " " << output[output.size()-1] << std::endl;

                }
            break;
          case 3:
            // set the header
            std::cout << "# x y z d g T ";

            for (unsigned int c = 0; c < compositions; ++c)
              std::cout << 'c' << c << ' ';

            for (unsigned int gc = 0; gc < grain_compositions; ++gc)
              for (size_t g = 0; g < n_grains; g++)
                std::cout << "gs" << gc << '-' << g << ' ' // gs = grain size, gm = grain rotation matrix
                          << "gm" << gc << '-' << g << "[0:0] " << "gm" << gc << '-' << g << "[0:1] " << "gm" << gc << '-' << g << "[0:2] "
                          << "gm" << gc << '-' << g << "[1:0] " << "gm" << gc << '-' << g << "[1:1] " << "gm" << gc << '-' << g << "[1:2] "
                          << "gm" << gc << '-' << g << "[2:0] " << "gm" << gc << '-' << g << "[2:1] " << "gm" << gc << '-' << g << "[2:2] ";

            std::cout << "tag ";
            std::cout <<std::endl;

            // set the values
            for (unsigned int i = 0; i < data.size(); ++i)
              if (data[i].size() > 0 && data[i][0] != "#")
                {
                  WBAssertThrow(data[i].size() == dim + 1, "The file needs to contain dim + 1 entries, but contains " << data[i].size() << " entries "
                                " on line " << i+1 << " of the data file (" << data_file << "). Dim is " << dim << '.');
                  std::array<double,3> coords = {{
                      string_to_double(data[i][0]), // x or R
                      string_to_double(data[i][1]) *(convert_spherical ? (WorldBuilder::Consts::PI/180.): 1.), // y or long
                      string_to_double(data[i][2]) *(convert_spherical ? (WorldBuilder::Consts::PI/180.): 1.) // z or lat
                    }
                  };

                  if (convert_spherical)
                    {
                      coords = spherical_to_cartesian_coordinates(coords).get_array();
                    }

                  std::cout << data[i][0] << ' ' << data[i][1] << ' ' << data[i][2] << ' ' << data[i][3] << ' ';
                  std::vector<double> output = world->properties(coords, string_to_double(data[i][3]),properties);
                  std::cout << output[0]  << ' ';

                  for (unsigned int c = 0; c < compositions; ++c)
                    {
                      std::cout << output[1+c]  << ' ';
                    }

                  for (unsigned int gc = 0; gc < grain_compositions; ++gc)
                    {
                      const size_t start = 1+compositions+gc*n_grains*10;
                      for (unsigned int g = 0; g < n_grains; ++g)
                        {
                          std::cout << output[start+g]  << ' '
                                    << output[start+n_grains+g*9] << ' ' << output[start+n_grains+g*9+1] << ' ' << output[start+n_grains+g*9+2] << ' '
                                    << output[start+n_grains+g*9+3] << ' ' << output[start+n_grains+g*9+4] << ' ' << output[start+n_grains+g*9+5] << ' '
                                    << output[start+n_grains+g*9+6] << ' ' << output[start+n_grains+g*9+7] << ' ' << output[start+n_grains+g*9+8] << ' ';

                        }
                    }
                  std::cout << " " << output[output.size()-1] << std::endl;

                }
            break;
          default:
            std::cout << "The World Builder can only be run in 2d and 3d but a different space dimension " << std::endl
                      << "is given: dim = " << dim << '.';

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
