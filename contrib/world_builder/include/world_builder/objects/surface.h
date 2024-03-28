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

#ifndef WORLD_BUILDER_OBJECTS_SURFACE_H
#define WORLD_BUILDER_OBJECTS_SURFACE_H

#include "world_builder/utilities.h"
#include "world_builder/kd_tree.h"

namespace WorldBuilder
{
  namespace Objects
  {
    struct SurfaceValueInfo
    {
      size_t triangle_index;
      double interpolated_value;
      double interpolator_s;
      double interpolator_t;

      SurfaceValueInfo(
        size_t triangle_index_,
        double interpolated_value_,
        double interpolator_s_,
        double interpolator_t_)
        :
        triangle_index(triangle_index_),
        interpolated_value(interpolated_value_),
        interpolator_s(interpolator_s_),
        interpolator_t(interpolator_t_) {};

      SurfaceValueInfo(double interpolated_value_)
        :
        triangle_index(NaN::IQNAN),
        interpolated_value(interpolated_value_),
        interpolator_s(NaN::DQNAN),
        interpolator_t(NaN::DQNAN) {};
    };

    class Surface
    {
      public:
        /**
         * Constructor to create an empty surface.
         */
        Surface();

        /**
         * Constructor to create a surface from value at points object output.
         */
        Surface(std::pair<std::vector<double>,std::vector<double>> values_at_points);

        /**
         * Returns the value of the surface at the check point.
         */
        SurfaceValueInfo local_value(const Point<2> &check_point) const;

        /**
         * Whether the surface is a constant value or not. This is used for optimalization.
         */
        bool constant_value;

        /**
         * The minimum value of all provided points.
         */
        double minimum;

        /**
         * The maximum value of all provided points.
         */
        double maximum;

        /**
         * The KD tree which stores the centroids of all triangles and an index to the triangle points
         * and values stored in the triangles member variable.
         */
        KDTree::KDTree tree;

        /**
         * Stores the triangles as a list of three points.
         */
        std::vector<std::array<std::array<double,3>,3> > triangles;

        /**
         * @brief Stores precomputed values
         */
        std::vector<std::array<double,8> > in_triangle_precomputed;

      private:

    };
  }

}

#endif