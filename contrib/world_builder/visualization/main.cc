/*
  Copyright (C) 2018 - 2020 by the authors of the World Builder code.

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
#include <cmath>

#include <algorithm>
#include <exception>
#include <iostream>
#include <array>
#include <fstream>
#include <thread>

#include <world_builder/assert.h>
#include <world_builder/nan.h>
#include <world_builder/utilities.h>
#include <world_builder/world.h>
#include <world_builder/coordinate_system.h>

#include <visualization/main.h>

using namespace WorldBuilder;
using namespace WorldBuilder::Utilities;


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
      // Deterimine the size of the slice for the loop
      size_t n = end - start + 1;
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
  WorldBuilder::Point<3> output_point(std::array<double,3> {{0,0,0}}, WorldBuilder::CoordinateSystem::cartesian);
  double r = in_point.norm();
  double theta = std::atan2(in_point[1],in_point[0]);
  double phi = std::acos(in_point[2]/r);

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
          const double pi4 = const_pi*0.25;
          const double x0 = -pi4 + static_cast<double>(i) * 2.0 * static_cast<double>(pi4)/static_cast<double>(level);
          const double y0 = -pi4 + static_cast<double>(j) * 2.0 * static_cast<double>(pi4)/static_cast<double>(level);
          const double r = std::tan(x0);
          const double s = std::tan(y0);


          double N1 = 0.25 * (1.0 - r) * (1.0 - s);
          double N2 = 0.25 * (1.0 + r) * (1.0 - s);
          double N3 = 0.25 * (1.0 + r) * (1.0 + s);
          double N4 = 0.25 * (1.0 - r) * (1.0 + s);

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
    vector.push_back(std::string(argv[i]));

  return vector;
}

bool find_command_line_option(char **begin, char **end, const std::string &option)
{
  return std::find(begin, end, option) != end;
}

int main(int argc, char **argv)
{
  /**
   * First parse the command line options
   */
  std::cout << "[1/5] Parsing file...                         \r";
  std::string wb_file;
  std::string data_file;

  size_t dim = 3;
  size_t compositions = 0;
  double gravity = 10;

  //commmon
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

  size_t number_of_threads = 1;

  try
    {

      if (find_command_line_option(argv, argv+argc, "-h") || find_command_line_option(argv, argv+argc, "--help"))
        {
          std::cout << "This program allows to use the world builder library directly with a world builder file and a grid file. "
                    "The data file will be filled with intitial conditions from the world as set by the world builder file." << std::endl
                    << "Besides providing two files, where the first is the world builder file and the second is the grid file, the available options are: " << std::endl
                    << "-h or --help to get this help screen," << std::endl
                    << "-j the number of threads the visualizer is allowed to use." << std::endl;
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
            }
        }


      if (options_vector.size() == 0)
        {
          std::cout << "Error: There where no files passed to the World Builder, use --help for more " << std::endl
                    << "information on how  to use the World Builder app." << std::endl;
          return 0;
        }


      if (options_vector.size() == 1)
        {
          std::cout << "Error:  The World Builder app requires at least two files, a World Builder file " << std::endl
                    << "and a data file to convert." << std::endl;
          return 0;
        }

      if (options_vector.size() != 2)
        {
          std::cout << "Only two command line arguments may be given, which should be the world builder file location and the grid file location (in that order). "
                    << "command line options where given." << std::endl;
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

  /**
   * Try to start the world builder
   */
  std::cout << "[2/5] Starting the world builder with " << number_of_threads << " threads...                         \r";
  std::cout.flush();

  std::unique_ptr<WorldBuilder::World> world;
  try
    {
      world = std::unique_ptr<WorldBuilder::World>(new WorldBuilder::World(wb_file));
    }
  catch (std::exception &e)
    {
      std::cerr << "Could not start the World builder from file '" << wb_file << "', error: " << e.what() << "\n";
      return 1;
    }
  catch (...)
    {
      std::cerr << "Exception of unknown type!\n";
      return 1;
    }

  /**
   * start the thread pool
   */
  ThreadPool pool(number_of_threads);

  /**
   * Read the data from the data files
   */
  std::cout << "[3/5] Reading grid file...                        \r";
  std::cout.flush();


  std::ifstream data_stream(data_file);

  // if config file is available, parse it
  WBAssertThrow(data_stream.good(),
                "Could not find the provided convig file at the specified location: " + data_file);


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
      for (size_t i = 0; i < line.size(); ++i)
        line[i].erase(std::remove(line[i].begin(), line[i].end(), ','), line[i].end());

      data.push_back(line);
    }

  // Read config from data if pressent
  for (size_t i = 0; i < data.size(); ++i)
    {
      if (data[i].size() == 0)
        continue;

      if (data[i][0] == "#")
        continue;

      if (data[i][0] == "grid_type" && data[i][1] == "=")
        {
          grid_type = data[i][2];
        }

      if (data[i][0] == "dim" && data[i][1] == "=")
        {
          dim = string_to_unsigned_int(data[i][2]);
        }

      if (data[i][0] == "compositions" && data[i][1] == "=")
        compositions = string_to_unsigned_int(data[i][2]);

      if (data[i][0] == "x_min" && data[i][1] == "=")
        x_min = string_to_double(data[i][2]);
      if (data[i][0] == "x_max" && data[i][1] == "=")
        x_max = string_to_double(data[i][2]);
      if (data[i][0] == "y_min" && data[i][1] == "=")
        y_min = string_to_double(data[i][2]);
      if (data[i][0] == "y_max" && data[i][1] == "=")
        y_max = string_to_double(data[i][2]);
      if (data[i][0] == "z_min" && data[i][1] == "=")
        z_min = string_to_double(data[i][2]);
      if (data[i][0] == "z_max" && data[i][1] == "=")
        z_max = string_to_double(data[i][2]);

      if (data[i][0] == "n_cell_x" && data[i][1] == "=")
        n_cell_x = string_to_unsigned_int(data[i][2]);
      if (data[i][0] == "n_cell_y" && data[i][1] == "=")
        n_cell_y = string_to_unsigned_int(data[i][2]);
      if (data[i][0] == "n_cell_z" && data[i][1] == "=")
        n_cell_z = string_to_unsigned_int(data[i][2]);

    }

  WBAssertThrow(dim == 2 || dim == 3, "dim should be set in the grid file and can only be 2 or 3.");

  WBAssertThrow(!std::isnan(x_min), "x_min is not a number:" << x_min << ". This value has probably not been provided in the grid file.");
  WBAssertThrow(!std::isnan(x_max), "x_max is not a number:" << x_max << ". This value has probably not been provided in the grid file.");
  WBAssertThrow(dim == 2 || !std::isnan(y_min), "y_min is not a number:" << y_min << ". This value has probably not been provided in the grid file.");
  WBAssertThrow(dim == 2 || !std::isnan(y_max), "y_max is not a number:" << y_max << ". This value has probably not been provided in the grid file.");
  WBAssertThrow(!std::isnan(z_min), "z_min is not a number:" << z_min << ". This value has probably not been provided in the grid file.");
  WBAssertThrow(!std::isnan(z_max), "z_max is not a number:" << z_max << ". This value has probably not been provided in the grid file.");


  WBAssertThrow(n_cell_x != 0, "n_cell_z may not be equal to zero: " << n_cell_x << ".");
  // int's cannot generally be nan's (see https://stackoverflow.com/questions/3949457/can-an-integer-be-nan-in-c),
  // but visual studio is giving problems over this, so it is taken out for now.
  //WBAssertThrow(!std::isnan(n_cell_x), "n_cell_z is not a number:" << n_cell_x << ".");

  WBAssertThrow(dim == 3 || n_cell_z != 0, "In 3d n_cell_z may not be equal to zero: " << n_cell_y << ".");
  // int's cannot generally be nan's (see https://stackoverflow.com/questions/3949457/can-an-integer-be-nan-in-c),
  // but visual studio is giving problems over this, so it is taken out for now.
  //WBAssertThrow(!std::isnan(n_cell_z), "n_cell_z is not a number:" << n_cell_y << ".");

  WBAssertThrow(n_cell_z != 0, "n_cell_z may not be equal to zero: " << n_cell_z << ".");
  // int's cannot generally be nan's (see https://stackoverflow.com/questions/3949457/can-an-integer-be-nan-in-c),
  // but visual studio is giving problems over this, so it is taken out for now.
  //WBAssertThrow(!std::isnan(n_cell_z), "n_cell_z is not a number:" << n_cell_z << ".");





  /**
   * All variables set by the user
   */

  if (grid_type == "sphere")
    WBAssert(n_cell_x == n_cell_y, "For the sphere grid the amount of cells in the x (long) and y (lat) direction have to be the same.");

  if (grid_type == "spherical" ||
      grid_type == "chunk" ||
      grid_type == "anullus")
    {
      x_min *= (const_pi/180);
      x_max *= (const_pi/180);
      y_min *= (const_pi/180);
      y_max *= (const_pi/180);
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


  bool compress_size = false;




  /**
   * Begin making the grid
   */
  std::cout << "[4/5] Building the grid...                        \r";
  std::cout.flush();
  WBAssertThrow(dim == 2 || dim == 3, "Dimension should be 2d or 3d.");
  if (grid_type == "cartesian")
    {
      n_cell = n_cell_x * n_cell_z * (dim == 3 ? n_cell_y : 1);
      if (compress_size == false && dim == 3)
        n_p = n_cell * 8 ; // it shouldn't matter for 2d in the output, so just do 3d.
      else
        n_p = (n_cell_x + 1) * (n_cell_z + 1) * (dim == 3 ? (n_cell_y + 1) : 1);


      double dx = (x_max - x_min) / static_cast<double>(n_cell_x);
      double dy = (y_max - y_min) / static_cast<double>(n_cell_y);
      double dz = (z_max - z_min) / static_cast<double>(n_cell_z);


      WBAssertThrow(!std::isnan(dx), "dz is not a number:" << dz << ".");
      WBAssertThrow(dim == 2 || !std::isnan(dy), "dz is not a number:" << dz << ".");
      WBAssertThrow(!std::isnan(dz), "dz is not a number:" << dz << ".");

      // todo: determine wheter a input variable is desirable for this.
      double surface = z_max;

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
          if (compress_size == true)
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
          if (compress_size == true)
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


      double inner_radius = z_min;
      double outer_radius = z_max;

      double l_outer = 2.0 * const_pi * outer_radius;

      double lr = outer_radius - inner_radius;
      double dr = lr / static_cast<double>(n_cell_z);

      size_t n_cell_t = static_cast<size_t>((2.0 * const_pi * outer_radius)/dr);

      // compute the ammount of cells
      n_cell = n_cell_t *n_cell_z;
      n_p = n_cell_t *(n_cell_z + 1);  // one less then cartesian because two cells overlap.

      double sx = l_outer / static_cast<double>(n_cell_t);
      double sz = dr;

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
              double xi = grid_x[counter];
              double zi = grid_z[counter];
              double theta = xi / l_outer * 2.0 * const_pi;
              grid_x[counter] = std::cos(theta) * (inner_radius + zi);
              grid_z[counter] = std::sin(theta) * (inner_radius + zi);
              grid_depth[counter] = outer_radius - std::sqrt(grid_x[counter] * grid_x[counter] + grid_z[counter] * grid_z [counter]);
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
      double inner_radius = z_min;
      double outer_radius = z_max;

      WBAssertThrow(x_min <= x_max, "The minimum longitude must be less than the maximum longitude.");
      WBAssertThrow(y_min <= y_max, "The minimum latitude must be less than the maximum latitude.");
      WBAssertThrow(inner_radius < outer_radius, "The inner radius must be less than the outer radius.");

      WBAssertThrow(x_min - x_max <= 2.0 * const_pi, "The difference between the minimum and maximum longitude "
                    " must be less than or equal to 360 degree.");

      WBAssertThrow(y_min >= - 0.5 * const_pi, "The minimum latitude must be larger then or equal to -90 degree.");
      WBAssertThrow(y_min <= 0.5 * const_pi, "The maximum latitude must be smaller then or equal to 90 degree.");

      double opening_angle_long_rad = (x_max - x_min);
      double opening_angle_lat_rad =  (y_max - y_min);

      n_cell = n_cell_x * n_cell_z * (dim == 3 ? n_cell_y : 1);
      if (compress_size == false && dim == 3)
        n_p = n_cell * 8 ; // it shouldn't matter for 2d in the output, so just do 3d.
      else
        n_p = (n_cell_x + 1) * (n_cell_z + 1) * (dim == 3 ? (n_cell_y + 1) : 1);

      double dlong = opening_angle_long_rad / static_cast<double>(n_cell_x);
      double dlat = opening_angle_lat_rad / static_cast<double>(n_cell_y);
      double lr = outer_radius - inner_radius;
      double dr = lr / static_cast<double>(n_cell_z);

      grid_x.resize(n_p);
      grid_y.resize(dim == 3 ? n_p : 0);
      grid_z.resize(n_p);
      grid_depth.resize(n_p);

      std::cout << "[4/5] Building the grid: stage 1 of 3                        \r";
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
          if (compress_size == true)
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

      std::cout << "[4/5] Building the grid: stage 2 of 3                        \r";
      std::cout.flush();
      if (dim == 2)
        {
          for (size_t i = 0; i < n_p; ++i)
            {

              double longitude = grid_x[i];
              double radius = grid_z[i];

              grid_x[i] = radius * std::cos(longitude);
              grid_z[i] = radius * std::sin(longitude);
            }
        }
      else
        {
          for (size_t i = 0; i < n_p; ++i)
            {

              double longitude = grid_x[i];
              double latitutde = grid_y[i];
              double radius = grid_z[i];

              grid_x[i] = radius * std::cos(latitutde) * std::cos(longitude);
              grid_y[i] = radius * std::cos(latitutde) * std::sin(longitude);
              grid_z[i] = radius * std::sin(latitutde);
            }
        }
      std::cout << "[4/5] Building the grid: stage 3 of 3                        \r";
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
                  std::cout << "[4/5] Building the grid: stage 3 of 3 [" << (static_cast<double>(i)/static_cast<double>(n_cell))*100.0 << "%]                       \r";
                  std::cout.flush();
                }
            }
        }
      else
        {
          if (compress_size == true)
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
                  std::cout << "[4/5] Building the grid: stage 3 of 3 [" << (static_cast<double>(i)/static_cast<double>(n_cell))*100.0 << "%]                       \r";
                  std::cout.flush();
                }
            }
        }
    }
  else if (grid_type == "sphere")
    {

      WBAssertThrow(dim == 3, "The sphere only works in 3d.");


      double inner_radius = z_min;
      double outer_radius = z_max;

      size_t n_block = 12;

      size_t block_n_cell = n_cell_x*n_cell_x;
      size_t block_n_p = (n_cell_x + 1) * (n_cell_x + 1);
      size_t block_n_v = 4;


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
          size_t block_n_cell_x = n_cell_x;
          size_t block_n_cell_y = n_cell_x;
          double Lx = 1.0;
          double Ly = 1.0;

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
      double distance = 1e-12*outer_radius;

      size_t counter = 0;
      size_t amount_of_double_points = 0;
      for (size_t i = 1; i < n_block * block_n_p; ++i)
        {
          if (sides[i])
            {
              double gxip = temp_x[i];
              double gyip = temp_y[i];
              double gzip = temp_z[i];
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


      size_t shell_n_p = n_block * block_n_p - amount_of_double_points;
      size_t shell_n_cell = n_block * block_n_cell;
      size_t shell_n_v = block_n_v;

      std::vector<double> shell_grid_x(shell_n_p);
      std::vector<double> shell_grid_y(shell_n_p);
      std::vector<double> shell_grid_z(shell_n_p);
      std::vector<std::vector<size_t> > shell_grid_connectivity(shell_n_cell,std::vector<size_t>(shell_n_v));

      counter = 0;
      for (size_t i = 0; i < n_block * block_n_p; ++i)
        {
          if (!double_points[i])
            {
              shell_grid_x[counter] = temp_x[i];
              shell_grid_y[counter] = temp_y[i];
              shell_grid_z[counter] = temp_z[i];

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

      size_t n_v = shell_n_v * 2;
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

          size_t i_beg =  i * shell_n_p;
          size_t i_end = (i+1) * shell_n_p;
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

              counter++;
            }
        }

      for (size_t i = 0; i < n_cell_z; ++i)
        {
          size_t i_beg = i * shell_n_cell;
          size_t i_end = (i+1) * shell_n_cell;
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

  // create paraview file.
  std::cout << "[5/5] Writing the paraview file...                                               \r";
  std::cout.flush();


  std::cout << "[5/5] Writing the paraview file: stage 1 of 3, writing header part 1                              \r";
  std::cout.flush();

  std::string base_filename = wb_file.substr(wb_file.find_last_of("/\\") + 1);
  std::string::size_type const p(base_filename.find_last_of('.'));
  std::string file_without_extension = base_filename.substr(0, p);

  std::ofstream myfile;
  myfile.open (file_without_extension + ".vtu");
  myfile << "<?xml version=\"1.0\" ?> " << std::endl;
  myfile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
  myfile << "<UnstructuredGrid>" << std::endl;
  myfile << "<FieldData>" << std::endl;
  myfile << "<DataArray type=\"Float32\" Name=\"TIME\" NumberOfTuples=\"1\" format=\"ascii\">0</DataArray>" << std::endl;
  myfile << "</FieldData>" << std::endl;
  myfile << "<Piece NumberOfPoints=\""<< n_p << "\" NumberOfCells=\"" << n_cell << "\">" << std::endl;
  myfile << "  <Points>" << std::endl;
  myfile << "    <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
  if (dim == 2)
    for (size_t i = 0; i < n_p; ++i)
      myfile << grid_x[i] << " " << grid_z[i] << " " << "0.0" << std::endl;
  else
    for (size_t i = 0; i < n_p; ++i)
      {
        myfile << grid_x[i] << " " << grid_y[i] << " " << grid_z[i] << std::endl;
      }
  std::cout << "[5/5] Writing the paraview file: stage 1 of 3, writing header part 2                              \r";
  std::cout.flush();
  myfile << "    </DataArray>" << std::endl;
  myfile << "  </Points>" << std::endl;
  myfile << std::endl;
  myfile << "  <Cells>" << std::endl;
  myfile << "    <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << std::endl;
  if (dim == 2)
    for (size_t i = 0; i < n_cell; ++i)
      myfile << grid_connectivity[i][0] << " " <<grid_connectivity[i][1] << " " << grid_connectivity[i][2] << " " << grid_connectivity[i][3] << std::endl;
  else
    for (size_t i = 0; i < n_cell; ++i)
      myfile << grid_connectivity[i][0] << " " <<grid_connectivity[i][1] << " " << grid_connectivity[i][2] << " " << grid_connectivity[i][3]  << " "
             << grid_connectivity[i][4] << " " <<grid_connectivity[i][5] << " " << grid_connectivity[i][6] << " " << grid_connectivity[i][7]<< std::endl;
  myfile << "    </DataArray>" << std::endl;
  myfile << "    <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << std::endl;
  if (dim == 2)
    for (size_t i = 1; i <= n_cell; ++i)
      myfile << i * 4 << " ";
  else
    for (size_t i = 1; i <= n_cell; ++i)
      myfile << i * 8 << " ";
  myfile << std::endl << "    </DataArray>" << std::endl;
  myfile << "    <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << std::endl;
  if (dim == 2)
    for (size_t i = 0; i < n_cell; ++i)
      myfile << "9" << " ";
  else
    for (size_t i = 0; i < n_cell; ++i)
      myfile << "12" << " ";
  myfile <<  std::endl <<"    </DataArray>" << std::endl;
  myfile << "  </Cells>" << std::endl;

  myfile << "  <PointData Scalars=\"scalars\">" << std::endl;

  myfile << "<DataArray type=\"Float32\" Name=\"Depth\" format=\"ascii\">" << std::endl;

  for (size_t i = 0; i < n_p; ++i)
    {
      myfile <<  grid_depth[i] << std::endl;
    }
  myfile << "</DataArray>" << std::endl;


  std::cout << "[5/5] Writing the paraview file: stage 2 of 3, computing temperatures                    \r";
  std::cout.flush();

  myfile << "    <DataArray type=\"Float32\" Name=\"Temperature\" format=\"ascii\">" << std::endl;
  std::vector<double> temp_vector(n_p);
  if (dim == 2)
    {
      pool.parallel_for(0, n_p, [&] (size_t i)
      {
        std::array<double,2> coords = {{grid_x[i], grid_z[i]}};
        temp_vector[i] = world->temperature(coords, grid_depth[i], gravity);
      });
    }
  else
    {
      pool.parallel_for(0, n_p, [&] (size_t i)
      {
        std::array<double,3> coords = {{grid_x[i], grid_y[i], grid_z[i]}};
        temp_vector[i] = world->temperature(coords, grid_depth[i], gravity);
      });
    }

  std::cout << "[5/5] Writing the paraview file: stage 2 of 3, writing temperatures                    \r";
  std::cout.flush();

  for (size_t i = 0; i < n_p; ++i)
    myfile << temp_vector[i]  << std::endl;
  myfile << "    </DataArray>" << std::endl;


  std::cout << "[5/5] Writing the paraview file: stage 3 of 3, writing compositions                     \r";
  std::cout.flush();

  for (size_t c = 0; c < compositions; ++c)
    {
      std::cout << "[5/5] Writing the paraview file: stage 2 of 3, computing composition "
                << c << " of " << compositions-1 << "            \r";
      std::cout.flush();

      myfile << "<DataArray type=\"Float32\" Name=\"Composition " << c << "\" Format=\"ascii\">" << std::endl;

      if (dim == 2)
        {
          pool.parallel_for(0, n_p, [&] (size_t i)
          {
            std::array<double,2> coords = {{grid_x[i], grid_z[i]}};
            temp_vector[i] =  world->composition(coords, grid_depth[i], static_cast<unsigned int>(c));
          });
        }
      else
        {
          pool.parallel_for(0, n_p, [&] (size_t i)
          {
            std::array<double,3> coords = {{grid_x[i], grid_y[i], grid_z[i]}};
            temp_vector[i] =  world->composition(coords, grid_depth[i], static_cast<unsigned int>(c));
          });
        }


      std::cout << "[5/5] Writing the paraview file: stage 2 of 3, writing composition "
                << c << " of " << compositions-1 << "            \r";
      std::cout.flush();

      for (size_t i = 0; i < n_p; ++i)
        myfile << temp_vector[i]  << std::endl;

      myfile << "</DataArray>" << std::endl;
    }

  myfile << "  </PointData>" << std::endl;


  myfile << " </Piece>" << std::endl;
  myfile << " </UnstructuredGrid>" << std::endl;
  myfile << "</VTKFile>" << std::endl;

  std::cout << "                                                                                \r";
  std::cout.flush();

  return 0;
}
