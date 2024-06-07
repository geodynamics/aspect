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

/**
 * Much of this code was based on GHOST (https://github.com/cedrict/GHOST),
 * and was contribute to the World Builder with the permission and help of
 * the author of GHOST.
 */

#include "visualization/main.h"

#include "world_builder/assert.h"
#include "world_builder/coordinate_system.h"
#include "world_builder/nan.h"
#include "world_builder/point.h"
#include "world_builder/utilities.h"
#include "world_builder/world.h"
#include "world_builder/config.h"

#include <algorithm>
#include <limits>


#include "vtu11/vtu11.hpp"
#undef max
#undef min

#ifdef WB_WITH_MPI
// we don't need the c++ MPI wrappers
#define OMPI_SKIP_MPICXX 1
#define MPICH_SKIP_MPICXX
#include <mpi.h>
#endif

#ifndef NDEBUG
#ifdef WB_USE_FP_EXCEPTIONS
#include <cfenv>
#endif
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <iterator>
#include <memory>
#include <string>
#include <thread>
#include <vector>

using namespace WorldBuilder;
using namespace WorldBuilder::Utilities;



/**
 * Filter the cells of a VTU mesh based on a given tag. All tags with smaller value than @p first_tag will be removed
 */
void filter_vtu_mesh(int dim,
                     const std::vector<bool> &include_tag,
                     vtu11::Vtu11UnstructuredMesh &input_mesh,
                     const std::vector<vtu11::DataSetData> &input_data,
                     vtu11::Vtu11UnstructuredMesh &output_mesh,
                     std::vector<vtu11::DataSetData> &output_data);

void filter_vtu_mesh(int dim,
                     const std::vector<bool> &include_tag,
                     vtu11::Vtu11UnstructuredMesh &input_mesh,
                     const std::vector<vtu11::DataSetData> &input_data,
                     vtu11::Vtu11UnstructuredMesh &output_mesh,
                     std::vector<vtu11::DataSetData> &output_data)
{
  output_data.resize(input_data.size());
  const std::int64_t invalid = -1;
  std::vector<std::int64_t> vertex_index_map(input_mesh.points().size(), invalid);

  const unsigned int n_vert_per_cell = (dim==3)?8:4;

  const std::size_t n_cells = input_mesh.types().size();
  std::uint64_t dst_cellid = 0;
  for (std::size_t cellidx = 0; cellidx <n_cells; ++cellidx)
    {
      int highest_tag = -1;
      for (size_t idx=cellidx*n_vert_per_cell; idx<(cellidx+1)*n_vert_per_cell; ++idx)
        {
          const std::size_t src_vid = static_cast<size_t>(input_mesh.connectivity()[idx]);
          highest_tag = std::max(highest_tag,static_cast<int>(input_data[2][src_vid]));
        }
      if (highest_tag < 0 || include_tag[static_cast<size_t>(highest_tag)]==false)
        continue;

      ++dst_cellid;

      for (size_t idx=cellidx*n_vert_per_cell; idx<(cellidx+1)*n_vert_per_cell; ++idx)
        {
          const size_t src_vid = static_cast<size_t>(input_mesh.connectivity()[idx]);

          std::int64_t dst_vid = vertex_index_map[src_vid];
          if (dst_vid == invalid)
            {
              dst_vid = static_cast<std::int64_t>(output_mesh.points().size()/3);
              vertex_index_map[src_vid] = dst_vid;

              for (unsigned int i=0; i<3; ++i)
                output_mesh.points().push_back(input_mesh.points()[src_vid*3+i]);

              for (unsigned int d=0; d<input_data.size(); ++d)
                output_data[d].push_back(input_data[d][src_vid]);
            }

          output_mesh.connectivity().push_back(dst_vid);

        }
      output_mesh.offsets().push_back(static_cast<long int>(dst_cellid*n_vert_per_cell));

      output_mesh.types().push_back(input_mesh.types()[cellidx]);
    }
}

/**
 * A very simple threadpool class. The threadpool currently only supports a
 * parallel for function, to easily parallelize the for function.
 *
 * Todo: might need something more advanced and safe in the future.
 */
class ThreadPool
{
  public:

    /**
     * Constructor
     */
    explicit ThreadPool(size_t number_of_threads)
    {
      pool.resize(number_of_threads);
    }

    /**
     * A function which allows to parallelize for loops.
     */
    template<typename Callable>
    void parallel_for(size_t start, size_t end, Callable func)
    {
      // Determine the size of the slice for the loop
      const size_t n = end - start + 1;
      size_t slice = static_cast<size_t>(std::round(n / static_cast<size_t>((pool.size()))));
      slice = std::max(slice, static_cast<size_t>(1));

      // Function which loops the passed function func
      auto loop_function = [&func] (size_t k1, size_t k2)
      {
        for (size_t k = k1; k < k2; k++)
          {
            func(k);
          }
      };

      // Launch jobs
      size_t i1 = start;
      size_t i2 = std::min(start + slice, end);
      for (size_t i = 0; i + 1 < pool.size() && i1 < end; ++i)
        {
          pool[i] = std::thread(loop_function, i1, i2);
          i1 = i2;
          i2 = std::min(i2 + slice, end);
        }
      if (i1 < end)
        {
          pool[pool.size()-1] = std::thread(loop_function, i1, end);
        }

      // Wait for jobs to finish
      for (std::thread &t : pool)
        {
          if (t.joinable())
            {
              t.join();
            }
        }
    }

  private:
    std::vector<std::thread> pool;

};


void project_on_sphere(double radius, double &x_, double &y_, double &z_)
{
  double x = x_;
  double y = y_;
  double z = z_;
  const WorldBuilder::Point<3> in_point(std::array<double,3> {{x,y,z}}, WorldBuilder::CoordinateSystem::cartesian);
  const WorldBuilder::Point<3> output_point(std::array<double,3> {{0,0,0}}, WorldBuilder::CoordinateSystem::cartesian);
  const double r = in_point.norm();
  const double theta = std::atan2(in_point[1],in_point[0]);
  const double phi = std::acos(in_point[2]/r);

  x_ = radius * std::cos(theta) * std::sin(phi);
  y_ = radius * std::sin(theta) * std::sin(phi);
  z_ = radius * std::cos(phi);

}

void lay_points(double x1, double y1, double z1,
                double x2, double y2, double z2,
                double x3, double y3, double z3,
                double x4, double y4, double z4,
                std::vector<double> &x, std::vector<double> &y, std::vector<double> &z,
                std::vector<bool> &hull, size_t level)
{
  // TODO: Assert that the vectors have the correct size;
  size_t counter = 0;
  for (size_t j = 0; j < level+1; ++j)
    {
      for (size_t i = 0; i < level+1; ++i)
        {
          // equidistant (is this irrelevant?)
          // double r = -1.0 + (2.0 / level) * i;
          // double s = -1.0 + (2.0 / level) * j;

          // equiangular
          const double pi4 = Consts::PI*0.25;
          const double x0 = -pi4 + static_cast<double>(i) * 2.0 * static_cast<double>(pi4)/static_cast<double>(level);
          const double y0 = -pi4 + static_cast<double>(j) * 2.0 * static_cast<double>(pi4)/static_cast<double>(level);
          const double r = std::tan(x0);
          const double s = std::tan(y0);


          const double N1 = 0.25 * (1.0 - r) * (1.0 - s);
          const double N2 = 0.25 * (1.0 + r) * (1.0 - s);
          const double N3 = 0.25 * (1.0 + r) * (1.0 + s);
          const double N4 = 0.25 * (1.0 - r) * (1.0 + s);

          x[counter] = x1 * N1 + x2 * N2 + x3 * N3 + x4 * N4;
          y[counter] = y1 * N1 + y2 * N2 + y3 * N3 + y4 * N4;
          z[counter] = z1 * N1 + z2 * N2 + z3 * N3 + z4 * N4;

          if (i == 0) hull[counter] = true;
          if (j == 0) hull[counter] = true;
          if (i == level) hull[counter] = true;
          if (j == level) hull[counter] = true;
          counter++;
        }
    }
}


std::vector<std::string> get_command_line_options_vector(int argc, char **argv)
{
  std::vector<std::string> vector;
  for (int i=1; i < argc; ++i)
    vector.emplace_back(argv[i]);

  return vector;
}

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

  // If set to true, we will output one visualization file per "tag"
  // with only the cells corresponding to that tag included.
  bool output_by_tag = false;
  // If set to true, we output a .filtered.vtu file without the background/mantle
  bool output_filtered = false;

  size_t dim = 3;
  size_t compositions = 0;

  // common
  std::string grid_type = "chunk";

  size_t n_cell_x = NaN::ISNAN; // x or long
  size_t n_cell_y = NaN::ISNAN; // y or lat
  size_t n_cell_z = NaN::ISNAN; // z or depth


  // spherical
  double x_min = NaN::DSNAN;  // x or long
  double x_max = NaN::DSNAN;  // x or long
  double y_min = NaN::DSNAN; // y or lat
  double y_max = NaN::DSNAN; // y or lat
  double z_min = NaN::DSNAN; // z or inner_radius
  double z_max = NaN::DSNAN; // z or outer_radius

  // Conservative choice for the number of threads to use:
  size_t number_of_threads = std::min(20u,1+std::thread::hardware_concurrency()/2);
  unsigned int max_resolution = std::numeric_limits<unsigned int>::max();

  try
    {
      if (find_command_line_option(argv, argv+argc, "-v") || find_command_line_option(argv, argv+argc, "--version"))
        {
          std::cout << "World Builder Grid Visualization tool.\n"
                    << "GWB Version: " << WorldBuilder::Version::MAJOR << "."
                    << WorldBuilder::Version::MINOR << "."
                    << WorldBuilder::Version::PATCH << "."
                    << WorldBuilder::Version::LABEL << "\n"
                    << "git hash: " << WorldBuilder::Version::GIT_SHA1 << " branch: " << WorldBuilder::Version::GIT_BRANCH
                    << std::endl;
          return 0;
        }

      if (find_command_line_option(argv, argv+argc, "-h") || find_command_line_option(argv, argv+argc, "--help"))
        {
          std::cout << "World Builder Grid Visualization tool.\n"
                    <<  "This program loads a world builder file and generates a visualization on a structured grid "
                    << "based on information specified in a separate .grid configuration file.\n\n"
                    << "Usage:\n"
                    << argv[0] << " [-j N] [--filtered] [--by-tag] example.wb example.grid\n\n"
                    << "Available options:\n"
                    << "  -j N                  Specify the number of threads the visualizer is allowed to use. Default: " << number_of_threads << ".\n"
                    << "  --filtered            Also produce a .filtered.vtu that removes cells only containing mantle or background.\n"
                    << "  --by-tag              Also produce a sequence of .N.vtu files that only contain cells of a specific tag.\n"
                    << "  --resolution-limit X  Specify a maximum resolution."
                    << "  -h or --help          To get this help screen.\n"
                    << "  -v or --version       To see version information.\n";
          return 0;
        }

      std::vector<std::string> options_vector = get_command_line_options_vector(argc, argv);

      for (size_t i = 0; i < options_vector.size(); ++i)
        {
          if (options_vector[i] == "-j")
            {
              number_of_threads = Utilities::string_to_unsigned_int(options_vector[i+1]);
              options_vector.erase(options_vector.begin()+static_cast<std::vector<std::string>::difference_type>(i));
              options_vector.erase(options_vector.begin()+static_cast<std::vector<std::string>::difference_type>(i));
              --i;
              continue;
            }
          if (options_vector[i] == "--filtered")
            {
              output_filtered = true;
              options_vector.erase(options_vector.begin()+static_cast<std::vector<std::string>::difference_type>(i));
              --i;
              continue;
            }
          if (options_vector[i] == "--by-tag")
            {
              output_by_tag = true;
              options_vector.erase(options_vector.begin()+static_cast<std::vector<std::string>::difference_type>(i));
              --i;
              continue;
            }
          if (options_vector[i] == "--resolution-limit")
            {
              max_resolution = Utilities::string_to_unsigned_int(options_vector[i+1]);
              options_vector.erase(options_vector.begin()+static_cast<std::vector<std::string>::difference_type>(i));
              options_vector.erase(options_vector.begin()+static_cast<std::vector<std::string>::difference_type>(i));
              --i;
              continue;
            }
        }

      if (options_vector.size() != 2)
        {
          std::cout << "World Builder Grid Visualization tool.\n"
                    << "Usage: " << argv[0] << " example.wb example.grid\n"
                    << "Try '" << argv[0] << " --help' for more information." << std::endl;
          return 0;
        }


      wb_file = options_vector[0];
      data_file = options_vector[1];

    }
  catch (std::exception &e)
    {
      std::cerr << "error: " << e.what() << "\n";
      return 1;
    }
  catch (...)
    {
      std::cerr << "Exception of unknown type!\n";
      return 1;
    }


  int MPI_RANK = 0;
#ifdef WB_WITH_MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &MPI_RANK);
#endif

  std::cout << "[1/6] Parsing file...                         \r";

  if (MPI_RANK ==  0)
    {
      /**
       * Try to start the world builder
       */
      std::cout << "[2/6] Starting the world builder with " << number_of_threads << " threads...                         \r";
      std::cout.flush();

      std::unique_ptr<WorldBuilder::World> world;
      try
        {
          world = std::make_unique<WorldBuilder::World>(wb_file);
        }
      catch (std::exception &e)
        {
          std::cerr << "Could not start the World builder from file '" << wb_file << "', error: " << e.what() << "\n";

#ifdef WB_WITH_MPI
          MPI_Finalize();
#endif
          return 1;
        }
      catch (...)
        {
          std::cerr << "Exception of unknown type!\n";

#ifdef WB_WITH_MPI
          MPI_Finalize();
#endif
          return 1;
        }

      /**
       * start the thread pool
       */
      ThreadPool pool(number_of_threads);

      /**
       * Read the data from the data files
       */
      std::cout << "[3/6] Reading grid file...                        \r";
      std::cout.flush();


      std::ifstream data_stream(data_file);

      // if config file is available, parse it
      WBAssertThrow(data_stream.good(),
                    "Could not find the provided config file at the specified location: " + data_file);


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

      std::string vtu_output_format = "RawBinaryCompressed";
      // Read config from data if present
      for (auto &line_i : data)
        {
          if (line_i.empty())
            continue;

          if (line_i[0] == "#")
            continue;

          if (line_i[0] == "grid_type" && line_i[1] == "=")
            {
              grid_type = line_i[2];
            }

          if (line_i[0] == "dim" && line_i[1] == "=")
            {
              dim = string_to_unsigned_int(line_i[2]);
            }

          if (line_i[0] == "vtu_output_format" && line_i[1] == "=")
            {
              vtu_output_format = line_i[2];
            }

          if (line_i[0] == "compositions" && line_i[1] == "=")
            compositions = string_to_unsigned_int(line_i[2]);

          if (line_i[0] == "x_min" && line_i[1] == "=")
            x_min = string_to_double(line_i[2]);
          if (line_i[0] == "x_max" && line_i[1] == "=")
            x_max = string_to_double(line_i[2]);
          if (line_i[0] == "y_min" && line_i[1] == "=")
            y_min = string_to_double(line_i[2]);
          if (line_i[0] == "y_max" && line_i[1] == "=")
            y_max = string_to_double(line_i[2]);
          if (line_i[0] == "z_min" && line_i[1] == "=")
            z_min = string_to_double(line_i[2]);
          if (line_i[0] == "z_max" && line_i[1] == "=")
            z_max = string_to_double(line_i[2]);

          if (line_i[0] == "n_cell_x" && line_i[1] == "=")
            n_cell_x = std::min(string_to_unsigned_int(line_i[2]),max_resolution);
          if (line_i[0] == "n_cell_y" && line_i[1] == "=")
            n_cell_y = std::min(string_to_unsigned_int(line_i[2]),max_resolution);
          if (line_i[0] == "n_cell_z" && line_i[1] == "=")
            n_cell_z = std::min(string_to_unsigned_int(line_i[2]),max_resolution);

        }

      WBAssertThrow(dim == 2 || dim == 3, "dim should be set in the grid file and can only be 2 or 3.");

      WBAssertThrow(!std::isnan(x_min), "x_min is not a number:" << x_min << ". This value has probably not been provided in the grid file.");
      WBAssertThrow(!std::isnan(x_max), "x_max is not a number:" << x_max << ". This value has probably not been provided in the grid file.");
      WBAssertThrow(dim == 2 || !std::isnan(y_min), "y_min is not a number:" << y_min << ". This value has probably not been provided in the grid file.");
      WBAssertThrow(dim == 2 || !std::isnan(y_max), "y_max is not a number:" << y_max << ". This value has probably not been provided in the grid file.");
      WBAssertThrow(!std::isnan(z_min), "z_min is not a number:" << z_min << ". This value has probably not been provided in the grid file.");
      WBAssertThrow(!std::isnan(z_max), "z_max is not a number:" << z_max << ". This value has probably not been provided in the grid file.");


      WBAssertThrow(n_cell_x != 0, "n_cell_z may not be equal to zero: " << n_cell_x << '.');
      // int's cannot generally be nan's (see https://stackoverflow.com/questions/3949457/can-an-integer-be-nan-in-c),
      // but visual studio is giving problems over this, so it is taken out for now.
      //WBAssertThrow(!std::isnan(n_cell_x), "n_cell_z is not a number:" << n_cell_x << '.');

      WBAssertThrow(dim == 3 || n_cell_z != 0, "In 3d n_cell_z may not be equal to zero: " << n_cell_y << '.');
      // int's cannot generally be nan's (see https://stackoverflow.com/questions/3949457/can-an-integer-be-nan-in-c),
      // but visual studio is giving problems over this, so it is taken out for now.
      //WBAssertThrow(!std::isnan(n_cell_z), "n_cell_z is not a number:" << n_cell_y << '.');

      WBAssertThrow(n_cell_z != 0, "n_cell_z may not be equal to zero: " << n_cell_z << '.');
      // int's cannot generally be nan's (see https://stackoverflow.com/questions/3949457/can-an-integer-be-nan-in-c),
      // but visual studio is giving problems over this, so it is taken out for now.
      //WBAssertThrow(!std::isnan(n_cell_z), "n_cell_z is not a number:" << n_cell_z << '.');





      /**
       * All variables set by the user
       */

      if (grid_type == "sphere")
        WBAssert(n_cell_x == n_cell_y, "For the sphere grid the amount of cells in the x (long) and y (lat) direction have to be the same.");

      if (grid_type == "spherical" ||
          grid_type == "chunk" ||
          grid_type == "annulus")
        {
          x_min *= (Consts::PI/180);
          x_max *= (Consts::PI/180);
          y_min *= (Consts::PI/180);
          y_max *= (Consts::PI/180);
        }



      /**
       * All variables needed for the visualization
       */
      size_t n_cell = NaN::ISNAN;
      size_t n_p = NaN::ISNAN;

      std::vector<double> grid_x(0);
      std::vector<double> grid_y(0);
      std::vector<double> grid_z(0);
      std::vector<double> grid_depth(0);

      std::vector<std::vector<size_t> > grid_connectivity(0);


      const bool compress_size = true;




      /**
       * Begin making the grid
       */
      std::cout << "[4/6] Building the grid...                        \r";
      std::cout.flush();
      WBAssertThrow(dim == 2 || dim == 3, "Dimension should be 2d or 3d.");
      if (grid_type == "cartesian")
        {
          n_cell = n_cell_x * n_cell_z * (dim == 3 ? n_cell_y : 1);
          if (!compress_size && dim == 3)
            n_p = n_cell * 8 ; // it shouldn't matter for 2d in the output, so just do 3d.
          else
            n_p = (n_cell_x + 1) * (n_cell_z + 1) * (dim == 3 ? (n_cell_y + 1) : 1);


          const double dx = (x_max - x_min) / static_cast<double>(n_cell_x);
          const double dy = dim == 2 ? 0 : (y_max - y_min) / static_cast<double>(n_cell_y);
          const double dz = (z_max - z_min) / static_cast<double>(n_cell_z);


          WBAssertThrow(!std::isnan(dx), "dz is not a number:" << dz << '.');
          WBAssertThrow(dim == 2 || !std::isnan(dy), "dz is not a number:" << dz << '.');
          WBAssertThrow(!std::isnan(dz), "dz is not a number:" << dz << '.');

          // todo: determine whether a input variable is desirable for this.
          const double surface = z_max;

          grid_x.resize(n_p);
          grid_z.resize(n_p);

          if (dim == 3)
            grid_y.resize(n_p);

          grid_depth.resize(n_p);

          // compute positions
          size_t counter = 0;
          if (dim == 2)
            {
              for (size_t j = 0; j <= n_cell_z; ++j)
                {
                  for (size_t i = 0; i <= n_cell_x; ++i)
                    {
                      grid_x[counter] = x_min + static_cast<double>(i) * dx;
                      grid_z[counter] = z_min + static_cast<double>(j) * dz;
                      grid_depth[counter] = (surface - z_min) - static_cast<double>(j) * dz;
                      counter++;
                    }
                }
            }
          else
            {
              if (compress_size)
                {
                  for (size_t i = 0; i <= n_cell_x; ++i)
                    {
                      for (size_t j = 0; j <= n_cell_y; ++j)
                        {
                          for (size_t k = 0; k <= n_cell_z; ++k)
                            {
                              grid_x[counter] = x_min + static_cast<double>(i) * dx;
                              grid_y[counter] = y_min + static_cast<double>(j) * dy;
                              grid_z[counter] = z_min + static_cast<double>(k) * dz;
                              grid_depth[counter] = (surface - z_min) - static_cast<double>(k) * dz;
                              counter++;
                            }
                        }
                    }
                }
              else
                {
                  for (size_t i = 0; i < n_cell_x; ++i)
                    {
                      for (size_t j = 0; j < n_cell_y; ++j)
                        {
                          for (size_t k = 0; k < n_cell_z; ++k)
                            {
                              // position is defined by the vtk file format
                              // position 0 of this cell
                              grid_x[counter] = x_min + static_cast<double>(i) * dx;
                              grid_y[counter] = y_min + static_cast<double>(j) * dy;
                              grid_z[counter] = z_min + static_cast<double>(k) * dz;
                              grid_depth[counter] = (surface - z_min) - static_cast<double>(k) * dz;
                              counter++;
                              // position 1 of this cell
                              grid_x[counter] = x_min + (static_cast<double>(i) + 1.0) * dx;
                              grid_y[counter] = y_min + static_cast<double>(j) * dy;
                              grid_z[counter] = z_min + static_cast<double>(k) * dz;
                              grid_depth[counter] = (surface - z_min) - static_cast<double>(k) * dz;
                              counter++;
                              // position 2 of this cell
                              grid_x[counter] = x_min + (static_cast<double>(i) + 1.0) * dx;
                              grid_y[counter] = y_min + (static_cast<double>(j) + 1.0) * dy;
                              grid_z[counter] = z_min + static_cast<double>(k) * dz;
                              grid_depth[counter] = (surface - z_min) - static_cast<double>(k) * dz;
                              counter++;
                              // position 3 of this cell
                              grid_x[counter] = x_min + static_cast<double>(i) * dx;
                              grid_y[counter] = y_min + (static_cast<double>(j) + 1.0) * dy;
                              grid_z[counter] = z_min + static_cast<double>(k) * dz;
                              grid_depth[counter] = (surface - z_min) - static_cast<double>(k) * dz;
                              counter++;
                              // position 0 of this cell
                              grid_x[counter] = x_min + static_cast<double>(i) * dx;
                              grid_y[counter] = y_min + static_cast<double>(j) * dy;
                              grid_z[counter] = z_min + (static_cast<double>(k) + 1.0) * dz;
                              grid_depth[counter] = (surface - z_min) - (static_cast<double>(k) + 1.0) * dz;
                              counter++;
                              // position 1 of this cell
                              grid_x[counter] = x_min + (static_cast<double>(i) + 1.0) * dx;
                              grid_y[counter] = y_min + static_cast<double>(j) * dy;
                              grid_z[counter] = z_min + (static_cast<double>(k) + 1.0) * dz;
                              grid_depth[counter] = (surface - z_min) - (static_cast<double>(k) + 1.0) * dz;
                              counter++;
                              // position 2 of this cell
                              grid_x[counter] = x_min + (static_cast<double>(i) + 1.0) * dx;
                              grid_y[counter] = y_min + (static_cast<double>(j) + 1.0) * dy;
                              grid_z[counter] = z_min + (static_cast<double>(k) + 1.0) * dz;
                              grid_depth[counter] = (surface - z_min) - (static_cast<double>(k) + 1.0) * dz;
                              counter++;
                              // position 3 of this cell
                              grid_x[counter] = x_min + static_cast<double>(i) * dx;
                              grid_y[counter] = y_min + (static_cast<double>(j) + 1.0) * dy;
                              grid_z[counter] = z_min + (static_cast<double>(k) + 1.0) * dz;
                              grid_depth[counter] = (surface - z_min) - (static_cast<double>(k) + 1.0) * dz;
                              WBAssert(counter < n_p, "Assert counter smaller then n_P: counter = " << counter << ", n_p = " << n_p);
                              counter++;
                            }
                        }
                    }
                }
            }

          // compute connectivity. Local to global mapping.
          grid_connectivity.resize(n_cell,std::vector<size_t>((dim-1)*4));

          counter = 0;
          if (dim == 2)
            {
              for (size_t j = 1; j <= n_cell_z; ++j)
                {
                  for (size_t i = 1; i <= n_cell_x; ++i)
                    {
                      grid_connectivity[counter][0] = i + (j - 1) * (n_cell_x + 1) - 1;
                      grid_connectivity[counter][1] = i + 1 + (j - 1) * (n_cell_x + 1) - 1;
                      grid_connectivity[counter][2] = i + 1  + j * (n_cell_x + 1) - 1;
                      grid_connectivity[counter][3] = i + j * (n_cell_x + 1) - 1;
                      counter++;
                    }
                }
            }
          else
            {
              if (compress_size)
                {
                  for (size_t i = 1; i <= n_cell_x; ++i)
                    {
                      for (size_t j = 1; j <= n_cell_y; ++j)
                        {
                          for (size_t k = 1; k <= n_cell_z; ++k)
                            {
                              grid_connectivity[counter][0] = (n_cell_y + 1) * (n_cell_z + 1) * (i - 1) + (n_cell_z + 1) * (j - 1) + k - 1;
                              grid_connectivity[counter][1] = (n_cell_y + 1) * (n_cell_z + 1) * (i    ) + (n_cell_z + 1) * (j - 1) + k - 1;
                              grid_connectivity[counter][2] = (n_cell_y + 1) * (n_cell_z + 1) * (i    ) + (n_cell_z + 1) * (j    ) + k - 1;
                              grid_connectivity[counter][3] = (n_cell_y + 1) * (n_cell_z + 1) * (i - 1) + (n_cell_z + 1) * (j    ) + k - 1;
                              grid_connectivity[counter][4] = (n_cell_y + 1) * (n_cell_z + 1) * (i - 1) + (n_cell_z + 1) * (j - 1) + k;
                              grid_connectivity[counter][5] = (n_cell_y + 1) * (n_cell_z + 1) * (i    ) + (n_cell_z + 1) * (j - 1) + k;
                              grid_connectivity[counter][6] = (n_cell_y + 1) * (n_cell_z + 1) * (i    ) + (n_cell_z + 1) * (j    ) + k;
                              grid_connectivity[counter][7] = (n_cell_y + 1) * (n_cell_z + 1) * (i - 1) + (n_cell_z + 1) * (j    ) + k;
                              counter++;
                            }
                        }
                    }
                }
              else
                {
                  for (size_t i = 0; i < n_cell; ++i)
                    {
                      grid_connectivity[i][0] = counter;
                      grid_connectivity[i][1] = counter + 1;
                      grid_connectivity[i][2] = counter + 2;
                      grid_connectivity[i][3] = counter + 3;
                      grid_connectivity[i][4] = counter + 4;
                      grid_connectivity[i][5] = counter + 5;
                      grid_connectivity[i][6] = counter + 6;
                      grid_connectivity[i][7] = counter + 7;
                      counter = counter + 8;
                    }
                }
            }
        }
      else if (grid_type == "annulus")
        {
          /**
           * An annulus which is a 2d hollow sphere.
           * TODO: make it so you can determine your own cross section.
           */
          WBAssertThrow(dim == 2, "The annulus only works in 2d.");


          const double inner_radius = z_min;
          const double outer_radius = z_max;

          const double l_outer = 2.0 * Consts::PI * outer_radius;

          const double lr = outer_radius - inner_radius;
          const double dr = lr / static_cast<double>(n_cell_z);

          const size_t n_cell_t = static_cast<size_t>((2.0 * Consts::PI * outer_radius)/dr);

          // compute the amount of cells
          n_cell = n_cell_t *n_cell_z;
          n_p = n_cell_t *(n_cell_z + 1);  // one less then cartesian because two cells overlap.

          const double sx = l_outer / static_cast<double>(n_cell_t);
          const double sz = dr;

          grid_x.resize(n_p);
          grid_z.resize(n_p);
          grid_depth.resize(n_p);

          size_t counter = 0;
          for (size_t j = 0; j <= n_cell_z; ++j)
            {
              for (size_t i = 1; i <= n_cell_t; ++i)
                {
                  grid_x[counter] = (static_cast<double>(i) - 1.0) * sx;
                  grid_z[counter] = static_cast<double>(j) * sz;
                  counter++;
                }
            }

          counter = 0;
          for (size_t j = 1; j <= n_cell_z+1; ++j)
            {
              for (size_t i = 1; i <= n_cell_t; ++i)
                {
                  const double xi = grid_x[counter];
                  const double zi = grid_z[counter];
                  const double theta = xi / l_outer * 2.0 * Consts::PI;
                  grid_x[counter] = std::cos(theta) * (inner_radius + zi);
                  grid_z[counter] = std::sin(theta) * (inner_radius + zi);
                  grid_depth[counter] = outer_radius - std::sqrt(grid_x[counter] * grid_x[counter] + grid_z[counter] * grid_z [counter]);
                  grid_depth[counter] = (std::fabs(grid_depth[counter]) < 1e-8 ? 0 : grid_depth[counter]);
                  counter++;
                }
            }

          grid_connectivity.resize(n_cell,std::vector<size_t>(4));
          counter = 0;
          for (size_t j = 1; j <= n_cell_z; ++j)
            {
              for (size_t i = 1; i <= n_cell_t; ++i)
                {
                  std::vector<size_t> cell_connectivity(4);
                  cell_connectivity[0] = counter + 1;
                  cell_connectivity[1] = counter + 1 + 1;
                  cell_connectivity[2] = i + j * n_cell_t + 1;
                  cell_connectivity[3] = i + j * n_cell_t;
                  if (i == n_cell_t)
                    {
                      cell_connectivity[1] = cell_connectivity[1] - n_cell_t;
                      cell_connectivity[2] = cell_connectivity[2] - n_cell_t;
                    }
                  grid_connectivity[counter][0] = cell_connectivity[1] - 1;
                  grid_connectivity[counter][1] = cell_connectivity[0] - 1;
                  grid_connectivity[counter][2] = cell_connectivity[3] - 1;
                  grid_connectivity[counter][3] = cell_connectivity[2] - 1;
                  counter++;
                }
            }
        }
      else if (grid_type == "chunk")
        {
          const double inner_radius = z_min;
          const double outer_radius = z_max;

          WBAssertThrow(x_min <= x_max, "The minimum longitude must be less than the maximum longitude.");
          WBAssertThrow(y_min <= y_max, "The minimum latitude must be less than the maximum latitude.");
          WBAssertThrow(inner_radius < outer_radius, "The inner radius must be less than the outer radius.");

          WBAssertThrow(x_min - x_max <= 2.0 * Consts::PI, "The difference between the minimum and maximum longitude "
                        " must be less than or equal to 360 degree.");

          WBAssertThrow(y_min >= - 0.5 * Consts::PI, "The minimum latitude must be larger then or equal to -90 degree.");
          WBAssertThrow(y_min <= 0.5 * Consts::PI, "The maximum latitude must be smaller then or equal to 90 degree.");

          const double opening_angle_long_rad = (x_max - x_min);
          const double opening_angle_lat_rad =  (y_max - y_min);

          n_cell = n_cell_x * n_cell_z * (dim == 3 ? n_cell_y : 1);
          if (!compress_size && dim == 3)
            n_p = n_cell * 8 ; // it shouldn't matter for 2d in the output, so just do 3d.
          else
            n_p = (n_cell_x + 1) * (n_cell_z + 1) * (dim == 3 ? (n_cell_y + 1) : 1);

          const double dlong = opening_angle_long_rad / static_cast<double>(n_cell_x);
          const double dlat = dim == 3 ? opening_angle_lat_rad / static_cast<double>(n_cell_y) : 0.;
          const double lr = outer_radius - inner_radius;
          const double dr = lr / static_cast<double>(n_cell_z);

          grid_x.resize(n_p);
          grid_y.resize(dim == 3 ? n_p : 0);
          grid_z.resize(n_p);
          grid_depth.resize(n_p);

          std::cout << "[4/6] Building the grid: stage 1 of 3                        \r";
          std::cout.flush();
          size_t counter = 0;
          if (dim == 2)
            {
              for (size_t i = 1; i <= n_cell_x + 1; ++i)
                for (size_t j = 1; j <= n_cell_z + 1; ++j)
                  {
                    grid_x[counter] = x_min + (static_cast<double>(i) - 1.0) * dlong;
                    grid_z[counter] = inner_radius + (static_cast<double>(j) - 1.0) * dr;
                    grid_depth[counter] = lr - (static_cast<double>(j) - 1.0) * dr;
                    counter++;
                  }
            }
          else
            {
              if (compress_size)
                {
                  for (size_t i = 1; i <= n_cell_x + 1; ++i)
                    for (size_t j = 1; j <= n_cell_y + 1; ++j)
                      for (size_t k = 1; k <= n_cell_z + 1; ++k)
                        {
                          grid_x[counter] = x_min + (static_cast<double>(i) - 1.0) * dlong;
                          grid_y[counter] = y_min + (static_cast<double>(j) - 1.0) * dlat;
                          grid_z[counter] = inner_radius + (static_cast<double>(k) - 1.0) * dr;
                          grid_depth[counter] = lr - (static_cast<double>(k) - 1.0) * dr;
                          counter++;
                        }
                }
              else
                {
                  for (size_t i = 0; i < n_cell_x; ++i)
                    {
                      for (size_t j = 0; j < n_cell_y; ++j)
                        {
                          for (size_t k = 0; k < n_cell_z; ++k)
                            {
                              // position is defined by the vtk file format
                              // position 0 of this cell
                              grid_x[counter] = x_min + static_cast<double>(i) * dlong;
                              grid_y[counter] = y_min + static_cast<double>(j) * dlat;
                              grid_z[counter] = inner_radius + static_cast<double>(k) * dr;
                              grid_depth[counter] = lr - static_cast<double>(k) * dr;
                              counter++;
                              // position 1 of this cell
                              grid_x[counter] = x_min + (static_cast<double>(i) + 1.0) * dlong;
                              grid_y[counter] = y_min + static_cast<double>(j) * dlat;
                              grid_z[counter] = inner_radius + static_cast<double>(k) * dr;
                              grid_depth[counter] = lr - static_cast<double>(k) * dr;
                              counter++;
                              // position 2 of this cell
                              grid_x[counter] = x_min + (static_cast<double>(i) + 1.0) * dlong;
                              grid_y[counter] = y_min + (static_cast<double>(j) + 1.0) * dlat;
                              grid_z[counter] = inner_radius + static_cast<double>(k) * dr;
                              grid_depth[counter] = lr - static_cast<double>(k) * dr;
                              counter++;
                              // position 3 of this cell
                              grid_x[counter] = x_min + static_cast<double>(i) * dlong;
                              grid_y[counter] = y_min + (static_cast<double>(j) + 1.0) * dlat;
                              grid_z[counter] = inner_radius + static_cast<double>(k) * dr;
                              grid_depth[counter] = lr - static_cast<double>(k) * dr;
                              counter++;
                              // position 0 of this cell
                              grid_x[counter] = x_min + static_cast<double>(i) * dlong;
                              grid_y[counter] = y_min + static_cast<double>(j) * dlat;
                              grid_z[counter] = inner_radius + (static_cast<double>(k) + 1.0) * dr;
                              grid_depth[counter] = lr - (static_cast<double>(k) + 1.0) * dr;
                              counter++;
                              // position 1 of this cell
                              grid_x[counter] = x_min + (static_cast<double>(i) + 1.0) * dlong;
                              grid_y[counter] = y_min + static_cast<double>(j) * dlat;
                              grid_z[counter] = inner_radius + (static_cast<double>(k) + 1.0) * dr;
                              grid_depth[counter] = lr - (static_cast<double>(k) + 1.0) * dr;
                              counter++;
                              // position 2 of this cell
                              grid_x[counter] = x_min + (static_cast<double>(i) + 1.0) * dlong;
                              grid_y[counter] = y_min + (static_cast<double>(j) + 1.0) * dlat;
                              grid_z[counter] = inner_radius + (static_cast<double>(k) + 1.0) * dr;
                              grid_depth[counter] = lr - (static_cast<double>(k) + 1.0) * dr;
                              counter++;
                              // position 3 of this cell
                              grid_x[counter] = x_min + static_cast<double>(i) * dlong;
                              grid_y[counter] = y_min + (static_cast<double>(j) + 1.0) * dlat;
                              grid_z[counter] = inner_radius + (static_cast<double>(k) + 1.0) * dr;
                              grid_depth[counter] = lr - (static_cast<double>(k) + 1.0) * dr;
                              WBAssert(counter < n_p, "Assert counter smaller then n_P: counter = " << counter << ", n_p = " << n_p);
                              counter++;
                            }
                        }
                    }
                }
            }

          std::cout << "[4/6] Building the grid: stage 2 of 3                        \r";
          std::cout.flush();
          if (dim == 2)
            {
              for (size_t i = 0; i < n_p; ++i)
                {

                  const double longitude = grid_x[i];
                  const double radius = grid_z[i];

                  grid_x[i] = radius * std::cos(longitude);
                  grid_z[i] = radius * std::sin(longitude);
                }
            }
          else
            {
              for (size_t i = 0; i < n_p; ++i)
                {

                  const double longitude = grid_x[i];
                  const double latitutde = grid_y[i];
                  const double radius = grid_z[i];

                  grid_x[i] = radius * std::cos(latitutde) * std::cos(longitude);
                  grid_y[i] = radius * std::cos(latitutde) * std::sin(longitude);
                  grid_z[i] = radius * std::sin(latitutde);
                }
            }
          std::cout << "[4/6] Building the grid: stage 3 of 3                        \r";
          std::cout.flush();
          // compute connectivity. Local to global mapping.
          grid_connectivity.resize(n_cell,std::vector<size_t>((dim-1)*4));

          counter = 0;
          if (dim == 2)
            {
              for (size_t i = 1; i <= n_cell_x; ++i)
                {
                  for (size_t j = 1; j <= n_cell_z; ++j)
                    {
                      grid_connectivity[counter][0] = (n_cell_z + 1) * (i - 1) + j - 1;
                      grid_connectivity[counter][1] = (n_cell_z + 1) * (i - 1) + j;
                      grid_connectivity[counter][2] = (n_cell_z + 1) * (i    ) + j;
                      grid_connectivity[counter][3] = (n_cell_z + 1) * (i    ) + j - 1;

                      counter = counter+1;
                      std::cout << "[4/6] Building the grid: stage 3 of 3 [" << (static_cast<double>(i)/static_cast<double>(n_cell))*100.0 << "%]                       \r";
                      std::cout.flush();
                    }
                }
            }
          else
            {
              if (compress_size)
                {
                  for (size_t i = 1; i <= n_cell_x; ++i)
                    {
                      for (size_t j = 1; j <= n_cell_y; ++j)
                        {
                          for (size_t k = 1; k <= n_cell_z; ++k)
                            {
                              grid_connectivity[counter][0] = (n_cell_y + 1) * (n_cell_z + 1) * (i - 1) + (n_cell_z + 1) * (j - 1) + k - 1;
                              grid_connectivity[counter][1] = (n_cell_y + 1) * (n_cell_z + 1) * (i    ) + (n_cell_z + 1) * (j - 1) + k - 1;
                              grid_connectivity[counter][2] = (n_cell_y + 1) * (n_cell_z + 1) * (i    ) + (n_cell_z + 1) * (j    ) + k - 1;
                              grid_connectivity[counter][3] = (n_cell_y + 1) * (n_cell_z + 1) * (i - 1) + (n_cell_z + 1) * (j    ) + k - 1;
                              grid_connectivity[counter][4] = (n_cell_y + 1) * (n_cell_z + 1) * (i - 1) + (n_cell_z + 1) * (j - 1) + k;
                              grid_connectivity[counter][5] = (n_cell_y + 1) * (n_cell_z + 1) * (i    ) + (n_cell_z + 1) * (j - 1) + k;
                              grid_connectivity[counter][6] = (n_cell_y + 1) * (n_cell_z + 1) * (i    ) + (n_cell_z + 1) * (j    ) + k;
                              grid_connectivity[counter][7] = (n_cell_y + 1) * (n_cell_z + 1) * (i - 1) + (n_cell_z + 1) * (j    ) + k;
                              counter++;
                            }
                        }
                    }
                }
              else
                {
                  for (size_t i = 0; i < n_cell; ++i)
                    {
                      grid_connectivity[i][0] = counter;
                      grid_connectivity[i][1] = counter + 1;
                      grid_connectivity[i][2] = counter + 2;
                      grid_connectivity[i][3] = counter + 3;
                      grid_connectivity[i][4] = counter + 4;
                      grid_connectivity[i][5] = counter + 5;
                      grid_connectivity[i][6] = counter + 6;
                      grid_connectivity[i][7] = counter + 7;
                      counter = counter + 8;
                      std::cout << "[4/6] Building the grid: stage 3 of 3 [" << (static_cast<double>(i)/static_cast<double>(n_cell))*100.0 << "%]                       \r";
                      std::cout.flush();
                    }
                }
            }
        }
      else if (grid_type == "sphere")
        {

          WBAssertThrow(dim == 3, "The sphere only works in 3d.");


          const double inner_radius = z_min;
          const double outer_radius = z_max;

          const size_t n_block = 12;

          const size_t block_n_cell = n_cell_x*n_cell_x;
          const size_t block_n_p = (n_cell_x + 1) * (n_cell_x + 1);
          const size_t block_n_v = 4;


          std::vector<std::vector<double> > block_grid_x(n_block,std::vector<double>(block_n_p));
          std::vector<std::vector<double> > block_grid_y(n_block,std::vector<double>(block_n_p));
          std::vector<std::vector<double> > block_grid_z(n_block,std::vector<double>(block_n_p));
          std::vector<std::vector<std::vector<size_t> > > block_grid_connectivity(n_block,std::vector<std::vector<size_t> >(block_n_cell,std::vector<size_t>(block_n_v)));
          std::vector<std::vector<bool> > block_grid_hull(n_block,std::vector<bool>(block_n_p));

          /**
           * block node layout
           */
          for (size_t i_block = 0; i_block < n_block; ++i_block)
            {
              const size_t block_n_cell_x = n_cell_x;
              const size_t block_n_cell_y = n_cell_x;
              const double Lx = 1.0;
              const double Ly = 1.0;

              size_t counter = 0;
              for (size_t j = 0; j <= block_n_cell_y; ++j)
                {
                  for (size_t i = 0; i <= block_n_cell_y; ++i)
                    {
                      block_grid_x[i_block][counter] = static_cast<double>(i) * Lx / static_cast<double>(block_n_cell_x);
                      block_grid_y[i_block][counter] = static_cast<double>(j) * Ly / static_cast<double>(block_n_cell_y);
                      block_grid_z[i_block][counter] = 0.0;
                      counter++;
                    }
                }

              counter = 0;
              // using i=1 and j=1 here because i an j are not used in lookup and storage
              // so the code can remain very similar to ghost and the cartesian code.
              for (size_t j = 1; j <= block_n_cell_y; ++j)
                {
                  for (size_t i = 1; i <= block_n_cell_x; ++i)
                    {
                      block_grid_connectivity[i_block][counter][0] = i + (j - 1) * (block_n_cell_x + 1) - 1;
                      block_grid_connectivity[i_block][counter][1] = i + 1 + (j - 1) * (block_n_cell_x + 1) - 1;
                      block_grid_connectivity[i_block][counter][2] = i + 1  + j * (block_n_cell_x + 1) - 1;
                      block_grid_connectivity[i_block][counter][3] = i + j * (block_n_cell_x + 1) - 1;
                      counter++;
                    }
                }
            }

          /**
           * map blocks
           */
          double radius = 1;

          // four corners
          double xA = -1.0;
          double yA = 0.0;
          double zA = -1.0 / std::sqrt(2.0);

          double xB = 1.0;
          double yB = 0.0;
          double zB = -1.0 / std::sqrt(2.0);

          double xC = 0.0;
          double yC = -1.0;
          double zC = 1.0 / std::sqrt(2.0);

          double xD = 0.0;
          double yD = 1.0;
          double zD = 1.0 / std::sqrt(2.0);

          // middles of faces
          double xM = (xA+xB+xC)/3.0;
          double yM = (yA+yB+yC)/3.0;
          double zM = (zA+zB+zC)/3.0;

          double xN = (xA+xD+xC)/3.0;
          double yN = (yA+yD+yC)/3.0;
          double zN = (zA+zD+zC)/3.0;

          double xP = (xA+xD+xB)/3.0;
          double yP = (yA+yD+yB)/3.0;
          double zP = (zA+zD+zB)/3.0;

          double xQ = (xC+xD+xB)/3.0;
          double yQ = (yC+yD+yB)/3.0;
          double zQ = (zC+zD+zB)/3.0;

          // middle of edges
          double xF = (xB+xC)/2.0;
          double yF = (yB+yC)/2.0;
          double zF = (zB+zC)/2.0;

          double xG = (xA+xC)/2.0;
          double yG = (yA+yC)/2.0;
          double zG = (zA+zC)/2.0;

          double xE = (xB+xA)/2.0;
          double yE = (yB+yA)/2.0;
          double zE = (zB+zA)/2.0;

          double xH = (xD+xC)/2.0;
          double yH = (yD+yC)/2.0;
          double zH = (zD+zC)/2.0;

          double xJ = (xD+xA)/2.0;
          double yJ = (yD+yA)/2.0;
          double zJ = (zD+zA)/2.0;

          double xK = (xD+xB)/2.0;
          double yK = (yD+yB)/2.0;
          double zK = (zD+zB)/2.0;

          // Making sure points A..Q are on a sphere
          project_on_sphere(radius,xA,yA,zA);
          project_on_sphere(radius,xB,yB,zB);
          project_on_sphere(radius,xC,yC,zC);
          project_on_sphere(radius,xD,yD,zD);
          project_on_sphere(radius,xE,yE,zE);
          project_on_sphere(radius,xF,yF,zF);
          project_on_sphere(radius,xG,yG,zG);
          project_on_sphere(radius,xH,yH,zH);
          project_on_sphere(radius,xJ,yJ,zJ);
          project_on_sphere(radius,xK,yK,zK);
          project_on_sphere(radius,xM,yM,zM);
          project_on_sphere(radius,xN,yN,zN);
          project_on_sphere(radius,xP,yP,zP);
          project_on_sphere(radius,xQ,yQ,zQ);

          lay_points(xM,yM,zM,xG,yG,zG,xA,yA,zA,xE,yE,zE,block_grid_x[0], block_grid_y[0], block_grid_z[0],block_grid_hull[0], n_cell_x);
          lay_points(xF,yF,zF,xM,yM,zM,xE,yE,zE,xB,yB,zB,block_grid_x[1], block_grid_y[1], block_grid_z[1],block_grid_hull[1], n_cell_x);
          lay_points(xC,yC,zC,xG,yG,zG,xM,yM,zM,xF,yF,zF,block_grid_x[2], block_grid_y[2], block_grid_z[2],block_grid_hull[2], n_cell_x);
          lay_points(xG,yG,zG,xN,yN,zN,xJ,yJ,zJ,xA,yA,zA,block_grid_x[3], block_grid_y[3], block_grid_z[3],block_grid_hull[3], n_cell_x);
          lay_points(xC,yC,zC,xH,yH,zH,xN,yN,zN,xG,yG,zG,block_grid_x[4], block_grid_y[4], block_grid_z[4],block_grid_hull[4], n_cell_x);
          lay_points(xH,yH,zH,xD,yD,zD,xJ,yJ,zJ,xN,yN,zN,block_grid_x[5], block_grid_y[5], block_grid_z[5],block_grid_hull[5], n_cell_x);
          lay_points(xA,yA,zA,xJ,yJ,zJ,xP,yP,zP,xE,yE,zE,block_grid_x[6], block_grid_y[6], block_grid_z[6],block_grid_hull[6], n_cell_x);
          lay_points(xJ,yJ,zJ,xD,yD,zD,xK,yK,zK,xP,yP,zP,block_grid_x[7], block_grid_y[7], block_grid_z[7],block_grid_hull[7], n_cell_x);
          lay_points(xP,yP,zP,xK,yK,zK,xB,yB,zB,xE,yE,zE,block_grid_x[8], block_grid_y[8], block_grid_z[8],block_grid_hull[8], n_cell_x);
          lay_points(xQ,yQ,zQ,xK,yK,zK,xD,yD,zD,xH,yH,zH,block_grid_x[9], block_grid_y[9], block_grid_z[9],block_grid_hull[9], n_cell_x);
          lay_points(xQ,yQ,zQ,xH,yH,zH,xC,yC,zC,xF,yF,zF,block_grid_x[10], block_grid_y[10], block_grid_z[10],block_grid_hull[10], n_cell_x);
          lay_points(xQ,yQ,zQ,xF,yF,zF,xB,yB,zB,xK,yK,zK,block_grid_x[11], block_grid_y[11], block_grid_z[11],block_grid_hull[11], n_cell_x);

          // make sure all points end up on a sphere
          for (size_t i_block = 0; i_block < n_block; ++i_block)
            {
              for (size_t i_point = 0; i_point < block_n_p; ++i_point)
                {
                  project_on_sphere(radius,block_grid_x[i_block][i_point],block_grid_y[i_block][i_point],block_grid_z[i_block][i_point]);
                }
            }

          /**
           * merge blocks
           */
          std::vector<double> temp_x(n_block * block_n_p);
          std::vector<double> temp_y(n_block * block_n_p);
          std::vector<double> temp_z(n_block * block_n_p);
          std::vector<bool> sides(n_block * block_n_p);

          for (size_t i = 0; i < n_block; ++i)
            {
              size_t counter = 0;
              for (size_t j = i * block_n_p; j < i * block_n_p + block_n_p; ++j)
                {
                  WBAssert(j < temp_x.size(), "j should be smaller then the size of the array temp_x.");
                  WBAssert(j < temp_y.size(), "j should be smaller then the size of the array temp_y.");
                  WBAssert(j < temp_z.size(), "j should be smaller then the size of the array temp_z.");
                  temp_x[j] = block_grid_x[i][counter];
                  temp_y[j] = block_grid_y[i][counter];
                  temp_z[j] = block_grid_z[i][counter];
                  sides[j] = block_grid_hull[i][counter];
                  counter++;
                }
            }


          std::vector<bool> double_points(n_block * block_n_p,false);
          std::vector<size_t> point_to(n_block * block_n_p);

          for (size_t i = 0; i < n_block * block_n_p; ++i)
            point_to[i] = i;

          // TODO: This becomes problematic with too large values of outer radius. Find a better way, maybe through an epsilon.
          const double distance = 1e-12*outer_radius;

          size_t counter = 0;
          size_t amount_of_double_points = 0;
          for (size_t i = 1; i < n_block * block_n_p; ++i)
            {
              if (sides[i])
                {
                  const double gxip = temp_x[i];
                  const double gyip = temp_y[i];
                  const double gzip = temp_z[i];
                  for (size_t j = 0; j < i-1; ++j)
                    {
                      if (sides[j])
                        {
                          if (std::fabs(gxip-temp_x[j]) < distance &&
                              std::fabs(gyip-temp_y[j]) < distance &&
                              std::fabs(gzip-temp_z[j]) < distance)
                            {
                              double_points[i] = true;
                              point_to[i] = j;
                              amount_of_double_points++;
                              break;
                            }
                        }
                    }
                }
            }


          const size_t shell_n_p = n_block * block_n_p - amount_of_double_points;
          const size_t shell_n_cell = n_block * block_n_cell;
          const size_t shell_n_v = block_n_v;

          std::vector<double> shell_grid_x(shell_n_p);
          std::vector<double> shell_grid_y(shell_n_p);
          std::vector<double> shell_grid_z(shell_n_p);
          std::vector<std::vector<size_t> > shell_grid_connectivity(shell_n_cell,std::vector<size_t>(shell_n_v));

          counter = 0;
          for (size_t i = 0; i < n_block * block_n_p; ++i)
            {
              if (!double_points[i])
                {
                  shell_grid_x[counter] = fabs(temp_x[i]) < 1e-8 ? 0. : temp_x[i];
                  shell_grid_y[counter] = fabs(temp_y[i]) < 1e-8 ? 0. : temp_y[i];
                  shell_grid_z[counter] = fabs(temp_z[i]) < 1e-8 ? 0. : temp_z[i];

                  counter++;
                }
            }

          for (size_t i = 0; i < n_block; ++i)
            {
              counter = 0;
              for (size_t j = i * block_n_cell; j < i * block_n_cell + block_n_cell; ++j)
                {
                  for (size_t k = 0; k < shell_n_v; ++k)
                    {
                      shell_grid_connectivity[j][k] = block_grid_connectivity[i][counter][k] + i * block_n_p;
                    }
                  counter++;
                }
            }

          for (size_t i = 0; i < shell_n_cell; ++i)
            {
              for (size_t j = 0; j < shell_n_v; ++j)
                {
                  shell_grid_connectivity[i][j] = point_to[shell_grid_connectivity[i][j]];
                }
            }

          std::vector<size_t> compact(n_block * block_n_p);

          counter = 0;
          for (size_t i = 0; i < n_block * block_n_p; ++i)
            {
              if (!double_points[i])
                {
                  compact[i] = counter;
                  counter++;
                }
            }


          for (size_t i = 0; i < shell_n_cell; ++i)
            {
              for (size_t j = 0; j < shell_n_v; ++j)
                {
                  shell_grid_connectivity[i][j] = compact[shell_grid_connectivity[i][j]];
                }
            }


          /**
           * build hollow sphere
           */

          std::vector<double> temp_shell_grid_x(shell_n_p);
          std::vector<double> temp_shell_grid_y(shell_n_p);
          std::vector<double> temp_shell_grid_z(shell_n_p);

          const size_t n_v = shell_n_v * 2;
          n_p = (n_cell_z + 1) * shell_n_p;
          n_cell = (n_cell_z) * shell_n_cell;

          grid_x.resize(n_p);
          grid_y.resize(n_p);
          grid_z.resize(n_p);
          grid_depth.resize(n_p);
          grid_connectivity.resize(n_cell,std::vector<size_t>(n_v));


          for (size_t i = 0; i < n_cell_z + 1; ++i)
            {
              temp_shell_grid_x = shell_grid_x;
              temp_shell_grid_y = shell_grid_y;
              temp_shell_grid_z = shell_grid_z;

              // We do not need to copy the shell_grid_connectivity into a
              // temperorary variable, because we do not change it. We can
              // directly use it.

              radius = inner_radius + ((outer_radius - inner_radius) / static_cast<double>(n_cell_z)) * static_cast<double>(i);
              for (size_t j = 0; j < shell_n_p; ++j)
                {
                  WBAssert(j < temp_shell_grid_x.size(), "ERROR: j = " << j << ", temp_shell_grid_x.size() = " << temp_shell_grid_x.size());
                  WBAssert(j < temp_shell_grid_y.size(), "ERROR: j = " << j << ", temp_shell_grid_y.size() = " << temp_shell_grid_y.size());
                  WBAssert(j < temp_shell_grid_z.size(), "ERROR: j = " << j << ", temp_shell_grid_z.size() = " << temp_shell_grid_z.size());
                  project_on_sphere(radius, temp_shell_grid_x[j], temp_shell_grid_y[j], temp_shell_grid_z[j]);
                }

              const size_t i_beg =  i * shell_n_p;
              const size_t i_end = (i+1) * shell_n_p;
              counter = 0;
              for (size_t j = i_beg; j < i_end; ++j)
                {
                  WBAssert(j < grid_x.size(), "ERROR: j = " << j << ", grid_x.size() = " << grid_x.size());
                  WBAssert(j < grid_y.size(), "ERROR: j = " << j << ", grid_y.size() = " << grid_y.size());
                  WBAssert(j < grid_z.size(), "ERROR: j = " << j << ", grid_z.size() = " << grid_z.size());
                  grid_x[j] = temp_shell_grid_x[counter];
                  grid_y[j] = temp_shell_grid_y[counter];
                  grid_z[j] = temp_shell_grid_z[counter];
                  grid_depth[j] = outer_radius - std::sqrt(grid_x[j] * grid_x[j] + grid_y[j] * grid_y[j] + grid_z[j] * grid_z[j]);
                  grid_depth[j] = (std::fabs(grid_depth[j]) < 1e-8 ? 0 : grid_depth[j]);

                  counter++;
                }
            }

          for (size_t i = 0; i < n_cell_z; ++i)
            {
              const size_t i_beg = i * shell_n_cell;
              const size_t i_end = (i+1) * shell_n_cell;
              counter = 0;
              for (size_t j = i_beg; j < i_end; ++j)
                {
                  for (size_t k = 0; k < shell_n_v; ++k)
                    {
                      grid_connectivity[j][k] = shell_grid_connectivity[counter][k] + i * shell_n_p;
                    }
                  counter++;
                }


              counter = 0;
              for (size_t j = i_beg; j < i_end; ++j)
                {
                  for (size_t k = shell_n_v ; k < 2 * shell_n_v; ++k)
                    {
                      WBAssert(k-shell_n_v < shell_grid_connectivity[counter].size(), "k - shell_n_v is larger then shell_grid_connectivity[counter]: k= " << k << ", shell_grid_connectivity[counter].size() = " << shell_grid_connectivity[counter].size());
                      grid_connectivity[j][k] = shell_grid_connectivity[counter][k-shell_n_v] + (i+1) * shell_n_p;
                    }
                  counter++;
                }
            }
        }
      else
        {
          WBAssertThrow(false, "Geometry type '" << grid_type << "' is not a valid geometry type. Valid geometry types are: "
                        << "'cartesian', 'annulus', 'chunk' and 'sphere'. "
                        << "Please note that the annulus can only be used in 2d and the sphere can only be used in 3d.");
        }

      // create paraview file.
      std::cout << "[5/6] Preparing to write the paraview file...                                                   \r";
      std::cout.flush();

      const std::string base_filename = wb_file.substr(wb_file.find_last_of("/\\") + 1);
      std::string::size_type  const p(base_filename.find_last_of('.'));
      const std::string file_without_extension = base_filename.substr(0, p);

      const std::stringstream buffer;
      const std::ofstream myfile;

      std::cout << "[5/6] Preparing to write the paraview file: stage 1 of 6, converting the points                              \r";
      std::cout.flush();
      std::vector<double> points(grid_x.size()*3, 0.0);
      if (dim == 2)
        for (size_t i = 0; i < n_p; ++i)
          {
            points[i*3] = grid_x[i];
            points[i*3+1] = grid_z[i];
            // third one is zero
          }
      else
        {
          for (size_t i = 0; i < n_p; ++i)
            {
              points[i*3] = grid_x[i];
              points[i*3+1] = grid_y[i];
              points[i*3+2] = grid_z[i];
            }
        }
      std::cout << "[5/6] Preparing to write the paraview file: stage 2 of 6, converting the connectivity                              \r";
      std::cout.flush();
      const size_t pow_2_dim = dim == 2 ? 4 : 8;
      std::vector<vtu11::VtkIndexType> connectivity(n_cell*pow_2_dim);
      if (dim == 2)
        for (size_t i = 0; i < n_cell; ++i)
          {
            connectivity[i*pow_2_dim] = static_cast<vtu11::VtkIndexType>(grid_connectivity[i][0]);
            connectivity[i*pow_2_dim+1] = static_cast<vtu11::VtkIndexType>(grid_connectivity[i][1]);
            connectivity[i*pow_2_dim+2] = static_cast<vtu11::VtkIndexType>(grid_connectivity[i][2]);
            connectivity[i*pow_2_dim+3] = static_cast<vtu11::VtkIndexType>(grid_connectivity[i][3]);
          }
      else
        for (size_t i = 0; i < n_cell; ++i)
          {

            connectivity[i*pow_2_dim] = static_cast<vtu11::VtkIndexType>(grid_connectivity[i][0]);
            connectivity[i*pow_2_dim+1] = static_cast<vtu11::VtkIndexType>(grid_connectivity[i][1]);
            connectivity[i*pow_2_dim+2] = static_cast<vtu11::VtkIndexType>(grid_connectivity[i][2]);
            connectivity[i*pow_2_dim+3] = static_cast<vtu11::VtkIndexType>(grid_connectivity[i][3]);
            connectivity[i*pow_2_dim+4] = static_cast<vtu11::VtkIndexType>(grid_connectivity[i][4]);
            connectivity[i*pow_2_dim+5] = static_cast<vtu11::VtkIndexType>(grid_connectivity[i][5]);
            connectivity[i*pow_2_dim+6] = static_cast<vtu11::VtkIndexType>(grid_connectivity[i][6]);
            connectivity[i*pow_2_dim+7] = static_cast<vtu11::VtkIndexType>(grid_connectivity[i][7]);
          }
      std::cout << "[5/6] Preparing to write the paraview file: stage 3 of 6, creating the offsets                              \r";
      std::cout.flush();
      std::vector<vtu11::VtkIndexType> offsets(n_cell);
      if (dim == 2)
        for (size_t i = 0; i < n_cell; ++i)
          offsets[i] = static_cast<vtu11::VtkIndexType>((i+1) * 4);
      else
        for (size_t i = 0; i < n_cell; ++i)
          offsets[i] = static_cast<vtu11::VtkIndexType>((i+1) * 8);

      std::cout << "[5/6] Preparing to write the paraview file: stage 4 of 6, creating the Data set info                              \r";
      std::cout.flush();
      std::vector<vtu11::VtkCellType> types(n_cell, dim == 2 ? 9 : 12);

      // Create tuples with (name, association, number of components) for each data set
      std::vector<vtu11::DataSetInfo> dataSetInfo
      {
        { "Depth", vtu11::DataSetType::PointData, 1 },
        { "Temperature", vtu11::DataSetType::PointData, 1 },
        { "Tag", vtu11::DataSetType::PointData, 1 },
      };
      for (size_t c = 0; c < compositions; ++c)
        {
          dataSetInfo.emplace_back( "Composition "+std::to_string(c), vtu11::DataSetType::PointData, 1 );
        }

      std::cout << "[5/6] Preparing to write the paraview file: stage 5 of 5, computing the properties                              \r";
      std::cout.flush();

      std::vector<std::array<unsigned ,3>> properties;
      properties.push_back({{1,0,0}}); // temperature

      properties.push_back({{4,0,0}}); // tag

      for (unsigned int c = 0; c < compositions; ++c)
        properties.push_back({{2,c,0}}); // composition c


      // compute temperature
      std::vector<vtu11::DataSetData> data_set(3+compositions);
      data_set[0] = grid_depth;
      data_set[1].resize(n_p);
      data_set[2].resize(n_p);
      for (size_t c = 0; c < compositions; ++c)
        data_set[3+c].resize(n_p);

      if (dim == 2)
        {
          pool.parallel_for(0, n_p, [&] (size_t i)
          {
            const std::array<double,2> coords = {{grid_x[i], grid_z[i]}};
            std::vector<double> output = world->properties(coords, grid_depth[i],properties);
            data_set[1][i] = output[0];
            data_set[2][i] = output[1];
            for (size_t c = 0; c < compositions; ++c)
              {
                data_set[3+c][i] = output[2+c];
              }
          });
        }
      else
        {
          pool.parallel_for(0, n_p, [&] (size_t i)
          {
            const std::array<double,3> coords = {{grid_x[i], grid_y[i], grid_z[i]}};
            std::vector<double> output = world->properties(coords, grid_depth[i],properties);
            data_set[1][i] = output[0];
            data_set[2][i] = output[1];
            for (size_t c = 0; c < compositions; ++c)
              {
                data_set[3+c][i] = output[2+c];
              }
          });
        }
      std::cout << "[6/6] Writing the paraview file                                                                                \r";
      std::cout.flush();

      {
        vtu11::Vtu11UnstructuredMesh mesh { points, connectivity, offsets, types };
        vtu11::writeVtu( file_without_extension + ".vtu", mesh, dataSetInfo, data_set, vtu_output_format );

        if (output_filtered)
          {
            std::vector<bool> include_tag(world->feature_tags.size(), true);
            for (unsigned int idx = 0; idx<include_tag.size(); ++idx)
              {
                if (world->feature_tags[idx]=="mantle layer")
                  include_tag[idx] = false;
              }
            std::vector<double> filtered_points;
            std::vector<vtu11::VtkIndexType> filtered_connectivity;
            std::vector<vtu11::VtkIndexType> filtered_offsets;
            std::vector<vtu11::VtkCellType> filtered_types;

            vtu11::Vtu11UnstructuredMesh filtered_mesh {filtered_points, filtered_connectivity, filtered_offsets, filtered_types};
            std::vector<vtu11::DataSetData> filtered_data_set;

            filter_vtu_mesh(static_cast<int>(dim), include_tag, mesh, data_set, filtered_mesh, filtered_data_set);
            vtu11::writeVtu( file_without_extension + ".filtered.vtu", filtered_mesh, dataSetInfo, filtered_data_set, vtu_output_format );
          }

        if (output_by_tag)
          {
            for (unsigned int idx = 0; idx<world->feature_tags.size(); ++idx)
              {
                if (world->feature_tags[idx]=="mantle layer")
                  continue;

                std::vector<double> filtered_points;
                std::vector<vtu11::VtkIndexType> filtered_connectivity;
                std::vector<vtu11::VtkIndexType> filtered_offsets;
                std::vector<vtu11::VtkCellType> filtered_types;

                vtu11::Vtu11UnstructuredMesh filtered_mesh {filtered_points, filtered_connectivity, filtered_offsets, filtered_types};
                std::vector<vtu11::DataSetData> filtered_data_set;

                std::vector<bool> include_tag(world->feature_tags.size(), false);
                include_tag[idx]=true;
                filter_vtu_mesh(static_cast<int>(dim), include_tag, mesh, data_set, filtered_mesh, filtered_data_set);
                const std::string filename = file_without_extension + "."+ std::to_string(idx)+".vtu";
                vtu11::writeVtu( filename, filtered_mesh, dataSetInfo, filtered_data_set, vtu_output_format );
              }
          }
      }
      std::cout << "                                                                                                               \r";
      std::cout.flush();
    }

#ifdef WB_WITH_MPI
  MPI_Finalize();
#endif
  return 0;
}
