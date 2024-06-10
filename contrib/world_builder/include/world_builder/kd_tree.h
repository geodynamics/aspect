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

#ifndef KD_TREE_H
#define KD_TREE_H

#include <vector>
#include <cstddef>
#include <iostream>
#include <limits>

#include "world_builder/assert.h"
#include "world_builder/point.h"

namespace WorldBuilder
{
  namespace KDTree
  {
    /**
     * A struct to store an index for a point and the
     * distance to that point.
     */
    struct IndexDistance
    {
      size_t index;
      double distance;
    };


    /**
     * A struct to store an index for a point and the
     * distance to that point.
     */
    struct IndexDistances
    {
      size_t min_index;
      double min_distance;
      std::vector<IndexDistance> vector;
    };

    /**
     * Each item contains an index and an x and y coordinate
     */
    struct Node
    {
      size_t index;
      double x;
      double y;

      Node(size_t index_, double x_, double y_)
        :
        index(index_),
        x(x_),
        y(y_)
      {}


      /**
       * access index (const)
       */
      inline
      const double &operator[](const bool y_axis) const
      {
        WBAssert(std::fabs((y_axis ? y : x) - *(&x+y_axis)) < std::numeric_limits<double>::epsilon(),
                 "Internal error: y_axis=" << y_axis << ", x=" << x << ", y=" << y <<", *(&x+y_axis)=" << *(&x+y_axis)
                 << ", ((bool)y_axis ? x : y) - *(&x+y_axis)=" << fabs((y_axis ? x : y) - *(&x+y_axis)));
        return *(&x+y_axis);
      }
    };

    class KDTree
    {
      public:
        /**
         * Constructor.
         */
        KDTree() = default;
        /**
         * Constructor. Requires a list of Nodes.
         */
        KDTree(const std::vector<Node> point_list);

        /**
         * Create a tree based on the current information in the
         * node vector.
         */
        void create_tree(const size_t left,
                         const size_t right,
                         const bool x_axis);

        /**
         * Return a reference to the vector containing the nodes.
         */
        const std::vector<Node> &get_nodes() const;

        /**
         * Returns the index of the closest point and the distance
         * of that point to the check point.
         */
        IndexDistance find_closest_point(const Point<2> &check_point) const;



        /**
         * Returns the index of the closest point and the distance
         * of that point to the check point. Stores the points searched
         * through in a unsorted vector.
         * Note: I can only guarantee that the point with the least distance
         * is the closest point. The point in the list with the second/third/etc.
         * smallest distance may or may not actually be the global point
         * with the second/third/etc. smallest distance. It will most likely be a
         * very good guess though.
         */
        IndexDistances find_closest_points(const Point<2> &check_point) const;


      private:
        /**
         * Returns the index of the closest point and the distance
         * of that point to the check point. This function is used
         * by find_closest_point to find the correct point.
         */
        void find_closest_point_recursive(const Point<2> &check_point,
                                          const size_t left,
                                          const size_t right,
                                          const bool y_axis,
                                          IndexDistance &index_distance) const;


        /**
         * Returns the index of the closest point and the distance
         * of that point to the check point. This function is used
         * by find_closest_point to find the correct point.
         */
        void find_closest_points_recursive(const Point<2> &check_point,
                                           const size_t left,
                                           const size_t right,
                                           const bool y_axis,
                                           IndexDistances &index_distances) const;

        std::vector<Node> nodes;
    };
  }
}
#endif