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

#include "world_builder/point.h"

#include <limits>

namespace WorldBuilder
{

  template<>
  Point<2>::Point(const double x, const double y, const CoordinateSystem coordinate_system_)
    :
    point({{x,y}}),
  coordinate_system(coordinate_system_)
  {}

  template<>
  Point<3>::Point(const double /*x*/, const double /*y*/, CoordinateSystem coordinate_system_)
    :
    point({{std::numeric_limits<double>::signaling_NaN(),std::numeric_limits<double>::signaling_NaN(), std::numeric_limits<double>::signaling_NaN()}}),
  coordinate_system(coordinate_system_)
  {
    WBAssertThrow(false,"Can't use the 2d constructor in 3d.");
  }


  template<>
  Point<2>::Point(const double /*x*/, const double /*y*/, const double /*z*/, CoordinateSystem coordinate_system_)
    :
    point({{std::numeric_limits<double>::signaling_NaN(),std::numeric_limits<double>::signaling_NaN()}}),
  coordinate_system(coordinate_system_)
  {
    WBAssertThrow(false,"Can't use the 3d constructor in 2d.");
  }


  template<>
  Point<3>::Point(const double x, const double y, const double z, CoordinateSystem coordinate_system_)
    :
    point({{x,y,z}}),
  coordinate_system(coordinate_system_)
  {}







  template<unsigned int dim>
  double
  Point<dim>::distance(const Point<2> &two) const
  {
    if (this->coordinate_system == spherical)
      {
        // spherical
        const double d_longitude = two[0] - this->point[0];
        const double d_latitude = two[1] - this->point[1];
        const double sin_d_lat = std::sin(d_latitude * 0.5);
        const double sin_d_long = std::sin(d_longitude * 0.5);
        return 2.0 * asin(sqrt((sin_d_lat * sin_d_lat) + (sin_d_long*sin_d_long) * std::cos(this->point[1]) * std::cos(two[1])));
      }

    // cartesian
    const double x_distance_to_reference_point = point[0]-two[0];
    const double y_distance_to_reference_point = point[1]-two[1];
    return sqrt((x_distance_to_reference_point*x_distance_to_reference_point) + (y_distance_to_reference_point*y_distance_to_reference_point));

  }

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
} // namespace WorldBuilder
