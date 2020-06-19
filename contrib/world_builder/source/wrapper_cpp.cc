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

#include "world_builder/wrapper_cpp.h"
#include "world_builder/world.h"
#include "iostream"

using namespace WorldBuilder;
namespace wrapper_cpp
{
  WorldBuilderWrapper::WorldBuilderWrapper(std::string filename, bool has_output_dir, std::string output_dir, const unsigned long random_number_seed)
    : ptr_ptr_world(NULL)
  {
    WorldBuilder::World *a = new WorldBuilder::World(filename, has_output_dir, output_dir, random_number_seed);
    ptr_ptr_world = reinterpret_cast<void *>(a);
  }


  WorldBuilderWrapper::~WorldBuilderWrapper()
  {
    WorldBuilder::World *a = reinterpret_cast<WorldBuilder::World *>(ptr_ptr_world);
    delete a;
  }


  double
  WorldBuilderWrapper::temperature_2d(double x, double z, double depth, double gravity)
  {
    std::array<double,2> position = {{x,z}};
    return reinterpret_cast<WorldBuilder::World *>(ptr_ptr_world)->temperature(position,depth,gravity);
  }

  double WorldBuilderWrapper::temperature_3d(double x, double y, double z, double depth, double gravity)
  {
    std::array<double,3> position = {{x,y,z}};
    return reinterpret_cast<WorldBuilder::World *>(ptr_ptr_world)->temperature(position,depth,gravity);
  }

  double WorldBuilderWrapper::composition_2d(double x, double z, double depth, unsigned int composition_number)
  {
    std::array<double,2> position = {{x,z}};
    return reinterpret_cast<WorldBuilder::World *>(ptr_ptr_world)->composition(position,depth,composition_number);
  }

  double WorldBuilderWrapper::composition_3d(double x, double y, double z, double depth, unsigned int composition_number)
  {
    std::array<double,3> position = {{x,y,z}};
    return reinterpret_cast<WorldBuilder::World *>(ptr_ptr_world)->composition(position,depth,composition_number);
  }
}