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

#include <limits>
#include <iostream>

#include <world_builder/point.h>

namespace WorldBuilder
{
  /**
   * The 2d version of the point class.
   */
  template class Point<2>;


  /**
   * The 3d version of the point class.
   */
  template class Point<3>;

  /**
   * Multiplies a 2d point with a scalar.
   */
  template Point<2> operator*(const double scalar, const Point<2> &point);

  /**
   * Multiplies a 3d point with a scalar.
   */
  template Point<3> operator*(const double scalar, const Point<3> &point);
}
