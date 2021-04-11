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
#include <world_builder/types/point.h>

namespace WorldBuilder
{
  namespace Types
  {
    template class Point<2>;
    template class Point<3>;

    /**
     * Multiplies a Types::Point<2> with a scalr and returns a
     * WorldBuilder::Point<2>.
     */
    template WorldBuilder::Point<2> operator*(const double scalar, const Point<2> &point);


    /**
     * Multiplies a Types::Point<3> with a scalr and returns a
     * WorldBuilder::Point<3>.
     */
    template WorldBuilder::Point<3> operator*(const double scalar, const Point<3> &point);
  }
}

