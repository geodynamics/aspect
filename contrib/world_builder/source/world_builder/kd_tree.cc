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

#include "world_builder/kd_tree.h"
#include <algorithm>
#include <utility>

namespace WorldBuilder
{
  namespace KDTree
  {
    KDTree::KDTree(std::vector<Node> point_list)
      : nodes(std::move(std::move(point_list)))
    {}


    void
    KDTree::create_tree(const size_t left,
                        const size_t right,
                        const bool y_axis)
    {
      const size_t mid=(left+right)>>1;
      std::nth_element(nodes.begin()+static_cast<long int>(left),nodes.begin()+static_cast<long int>(mid),nodes.begin()+static_cast<long int>(right)+1,
                       [y_axis](Node& i, Node& j) -> bool
      {
        return i[y_axis] < j[y_axis];
      });

      if (left<mid)
        create_tree(left,mid-1,!y_axis);
      if (right>mid)
        create_tree(mid+1,right,!y_axis);
    }


    const std::vector<Node> &
    KDTree::get_nodes() const
    {
      return nodes;
    }


    IndexDistance
    KDTree::find_closest_point(const Point<2> &check_point) const
    {
      const size_t start_node_index = 0;
      IndexDistance index_distance = {start_node_index,
                                      std::numeric_limits<double>::max()
                                     };

      find_closest_point_recursive(check_point,0,nodes.size()-1,false,index_distance);

      return index_distance;
    }


    void
    KDTree::find_closest_point_recursive(const Point<2> &check_point,
                                         const size_t left,
                                         const size_t right,
                                         const bool y_axis,
                                         IndexDistance &index_distance) const
    {
      // Calculate the index of this node
      const size_t mid=(left+right)>>1;
      const Node &node = nodes[mid];
      if (check_point[static_cast<size_t>(y_axis)] < node[y_axis])
        {
          // Traverse left child
          if (left<mid)
            {
              find_closest_point_recursive(check_point,left,mid-1,!y_axis,index_distance);
            }

          // Compare node's point to current closest point
          const double distance = sqrt((node[false]-check_point[0])*(node[false]-check_point[0])
                                       +(node[true]-check_point[1])*(node[true]-check_point[1]));
          if (index_distance.distance > distance)
            {
              index_distance.index = mid;
              index_distance.distance = distance;
            }

          // Traverse right child
          if (right>mid)
            {
              if ((node[y_axis]-check_point[static_cast<size_t>(y_axis)]) < index_distance.distance)
                {
                  find_closest_point_recursive(check_point,mid+1,right,!y_axis,index_distance);
                }
            }
        }
      else
        {
          // Traverse right child
          if (right>mid)
            find_closest_point_recursive(check_point,mid+1,right,!y_axis,index_distance);

          // Compare node's point to current closest point
          const double distance = sqrt((node[false]-check_point[0])*(node[false]-check_point[0])
                                       +(node[true]-check_point[1])*(node[true]-check_point[1]));
          if (index_distance.distance > distance)
            {
              index_distance.index = mid;
              index_distance.distance = distance;
            }

          // Traverse left child
          if (left<mid)
            if ((node[y_axis]-check_point[static_cast<size_t>(y_axis)]) < index_distance.distance)
              find_closest_point_recursive(check_point,left,mid-1,!y_axis,index_distance);
        }
    }



    IndexDistances
    KDTree::find_closest_points(const Point<2> &check_point) const
    {
      const size_t start_node_index = 0;
      IndexDistances index_distances = {start_node_index,
                                        std::numeric_limits<double>::max(),
                                        {}
                                       };
      index_distances.vector.reserve(nodes.size());

      find_closest_points_recursive(check_point,0,nodes.size()-1,false,index_distances);

      return index_distances;
    }


    void
    KDTree::find_closest_points_recursive(const Point<2> &check_point,
                                          const size_t left,
                                          const size_t right,
                                          const bool y_axis,
                                          IndexDistances &index_distances) const
    {
      // Calculate the index of this node
      const size_t mid=(left+right)>>1;
      const Node &node = nodes[mid];
      if (check_point[static_cast<size_t>(y_axis)] < node[y_axis])
        {
          // Traverse left child
          if (left<mid)
            {
              find_closest_points_recursive(check_point,left,mid-1,!y_axis,index_distances);
            }

          // Compare node's point to current closest point
          const double distance = sqrt((node[false]-check_point[0])*(node[false]-check_point[0])
                                       +(node[true]-check_point[1])*(node[true]-check_point[1]));
          if (index_distances.min_distance > distance)
            {
              index_distances.min_index = mid;
              index_distances.min_distance = distance;
            }

          index_distances.vector.emplace_back(IndexDistance {mid, distance});

          // Traverse right child
          if (right>mid)
            {
              if ((node[y_axis]-check_point[static_cast<size_t>(y_axis)]) < index_distances.min_distance)
                {
                  find_closest_points_recursive(check_point,mid+1,right,!y_axis,index_distances);
                }
            }
        }
      else
        {
          // Traverse right child
          if (right>mid)
            find_closest_points_recursive(check_point,mid+1,right,!y_axis,index_distances);

          // Compare node's point to current closest point
          const double distance = sqrt((node[false]-check_point[0])*(node[false]-check_point[0])
                                       +(node[true]-check_point[1])*(node[true]-check_point[1]));
          if (index_distances.min_distance > distance)
            {
              index_distances.min_index = mid;
              index_distances.min_distance = distance;
            }

          index_distances.vector.emplace_back(IndexDistance {mid, distance});

          // Traverse left child
          if (left<mid)
            if ((node[y_axis]-check_point[static_cast<size_t>(y_axis)]) < index_distances.min_distance)
              find_closest_points_recursive(check_point,left,mid-1,!y_axis,index_distances);
        }
    }
  } // namespace KDTree
} // namespace WorldBuilder
