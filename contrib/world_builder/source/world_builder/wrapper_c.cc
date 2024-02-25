/*
  Copyright (C) 2018-2024 by the authors of the World Builder code.

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

#include "world_builder/wrapper_c.h"

#include "world_builder/world.h"

extern "C" {
  /**
   * This function creates an object of the world builder and returns a pointer
   * to it. This pointer can then be used to call the temperature and composition
   * functions. When done call the release world function to destroy the object.
   */
  void create_world(void **ptr_ptr_world, const char *world_builder_file, const bool *has_output_dir_, const char *output_dir_, const unsigned long random_number_seed)
  {
    bool has_output_dir = false;

    if (has_output_dir_ != nullptr)
      {
        has_output_dir = *has_output_dir_;
      }

    std::string output_dir;
    if (output_dir_ != nullptr)
      {
        output_dir = *output_dir_;
      }

    WorldBuilder::World *a = new WorldBuilder::World(std::string(world_builder_file), has_output_dir, output_dir,random_number_seed);

    *ptr_ptr_world = reinterpret_cast<void *>(a);
  }


  /**
   * This function return the temperature at a specific location given x, z, depth and
   * gravity.
   */
  void temperature_2d(void *ptr_ptr_world, double x, double z, double depth, double *temperature)
  {
    WorldBuilder::World *a = reinterpret_cast<WorldBuilder::World *>(ptr_ptr_world);
    const std::array<double,2> position = {{x,z}};
    *temperature = a->temperature(position,depth);
  }


  /**
   * This function return the temperature at a specific location given x, y, z, depth and
   * gravity.
   */
  void temperature_3d(void *ptr_ptr_world, double x, double y, double z, double depth, double *temperature)
  {
    WorldBuilder::World *a = reinterpret_cast<WorldBuilder::World *>(ptr_ptr_world);
    const std::array<double,3> position = {{x,y,z}};
    *temperature = a->temperature(position,depth);
  }


  /**
   * This function return the composition at a specific location given x, z, depth and
   * composition number.
   */
  void composition_2d(void *ptr_ptr_world, double x, double z, double depth, unsigned int composition_number, double *composition)
  {
    WorldBuilder::World *a = reinterpret_cast<WorldBuilder::World *>(ptr_ptr_world);
    const std::array<double,2> position = {{x,z}};
    *composition = a->composition(position,depth,composition_number);
  }


  /**
   * This function return the composition at a specific location given x, y, z, depth and
   * composition number.
   */
  void composition_3d(void *ptr_ptr_world, double x, double y, double z, double depth, unsigned int composition_number, double *composition)
  {
    WorldBuilder::World *a = reinterpret_cast<WorldBuilder::World *>(ptr_ptr_world);
    const std::array<double,3> position = {{x,y,z}};
    *composition = a->composition(position,depth,composition_number);
  }


  /**
   * The destructor for the world builder class. Call this function when done with the
   * world builder.
   */
  void release_world(void *ptr_ptr_world)
  {
    WorldBuilder::World *a = reinterpret_cast<WorldBuilder::World *>(ptr_ptr_world);
    delete a;
  }
}
