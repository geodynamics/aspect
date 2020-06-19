/*
  Copyright (C) 2020 by the authors of the World Builder code.

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

#ifndef _world_builder_grains_h
#define _world_builder_grains_h
#include <vector>
#include <array>

namespace WorldBuilder
{
  /**
   * This is a simple structure to store information about grains.
   * The advantage of storing all grains in seperate vectors, compared
   * to having a vector of individual grains, is that the vectors can
   * be empty if the information is not needed.
   */
  struct grains
  {
    // The sizes of the grains
    std::vector<double> sizes;

    // the rotation matrices of the latices of the grains.
    // todo: convention.
    std::vector<std::array<std::array<double,3>,3> > rotation_matrices;
  };
}

#endif