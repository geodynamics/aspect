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

#ifndef _world_builder_wrapper_cpp_h
#define _world_builder_wrapper_cpp_h

#include <string>

namespace wrapper_cpp
{
  /**
   * This class is to be used by SWIG. To make it easy for SWIG we do not use any
   * other world builder header file in this header file. This means that the class
   * stores a void pointer internally. The cpp implementation can use the world builder
   * header files.
   */
  class WorldBuilderWrapper
  {
    public:
      /**
       * constructor
       */
      WorldBuilderWrapper(std::string filename, bool has_output_dir = false, std::string output_dir = "", const unsigned long random_number_seed = 1.0);

      /**
       * destructor
       */
      ~WorldBuilderWrapper();

      /**
       * This function return the temperature at a specific location given x, z, depth and
       * gravity.
       */
      double temperature_2d(double x, double z, double depth, double gravity);

      /**
       * This function return the temperature at a specific location given x, y, z, depth and
       * gravity.
       */
      double temperature_3d(double x, double y, double z, double depth, double gravity);

      /**
       * This function return the composition at a specific location given x, z, depth and
       * composition number.
       */
      double composition_2d(double x, double z, double depth, unsigned int composition_number);

      /**
       * This function return the composition at a specific location given x, y, z, depth and
       * composition number.
       */
      double composition_3d(double x, double y, double z, double depth, unsigned int composition_number);


    private:
      void *ptr_ptr_world;


  };
}
#endif