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

#ifndef WORLD_BUILDER_OBJECTS_BEZIER_CURVE_H
#define WORLD_BUILDER_OBJECTS_BEZIER_CURVE_H

#include "world_builder/objects/closest_point_on_curve.h"
#include "world_builder/point.h"
#include <array>
#include <vector>

namespace WorldBuilder
{
  namespace Objects
  {

    /**
     * @brief Class for circle line/spline, including interpolation on it
     *
     */
    class BezierCurve
    {
      public:
        /**
         * @brief Construct a new Bezier Curve object
         *
         * @param p
         * @param angle_constrains
         */
        BezierCurve() = default;

        /**
         * @brief Construct a new Bezier Curve object
         *
         * @param p
         * @param angle_constrains
         */
        BezierCurve(const std::vector<Point<2> > &p, const std::vector<double> &angle_constrains = {});

        /**
         * @brief Finds the closest point on the curve. If the the closest point
         *        doesn't fall on the segment, return a point with x and y being nan.
         *
         * @param p
         * @param verbose Whether this function should be outputting its Newton iteration
         * to std::cout while running. This is very expensive, but useful for debugging
         * purposes.
         * @return ClosestPointOnCurve
         */
        ClosestPointOnCurve closest_point_on_curve_segment(const Point<2> &p, const bool verbose = false) const;

        /**
         * @brief
         *
         * @param i
         * @param x
         * @return Point<2>
         */
        Point<2> operator()(const size_t i, const double x) const;

      private:
        std::vector<Point<2> > points;
        std::vector<std::array<Point<2>,2 > > control_points;
        std::vector<double> lengths;
        std::vector<double> angles;

    };
  }

}


#endif
