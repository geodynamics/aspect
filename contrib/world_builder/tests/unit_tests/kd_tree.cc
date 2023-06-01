/*
  Copyright (C) 2022 by the authors of the World Builder code.

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

#define DOCTEST_CONFIG_SUPER_FAST_ASSERTS

#include "doctest/doctest.h"
#include "world_builder/coordinate_system.h"
#include "world_builder/kd_tree.h"

#include <iostream>

using namespace WorldBuilder::KDTree;
using doctest::Approx;

/**
 * Compare the given two std::array<double,3> entries with an epsilon (using Catch::Approx)
 */
inline void compare_node_trees(
  const std::vector<Node> &computed,
  const std::vector<Node> &expected)
{
  for (size_t i = 0; i < computed.size(); ++i)
    {
      CHECK(computed[i].index == expected[i].index);
      CHECK(computed[i].x == Approx(expected[i].x));
      CHECK(computed[i].y == Approx(expected[i].y));
    }
}

TEST_CASE("create kd-tree: 15 nodes")
{
  const std::vector<Node> reference = {Node(8,2,3),Node(4,4,4),Node(9,5,2),
                                       Node(2,1,5),Node(10,2,9),Node(5,5,8),
                                       Node(11,6,7),Node(1,7,3),Node(12,9,2.5),
                                       Node(6,10,2),Node(13,11.5,1),Node(3,11,4),
                                       Node(14,8,6),Node(7,9,7),Node(15,11,8)
                                      };
  {
    const std::vector<Node> nodes = {Node(1,7,3),Node(2,1,5),Node(3,11,4),
                                     Node(4,4,4),Node(5,5,8),Node(6,10,2),
                                     Node(7,9,7),Node(8,2,3),Node(9,5,2),
                                     Node(10,2,9),Node(11,6,7),Node(12,9,2.5),
                                     Node(13,11.5,1),Node(14,8,6),Node(15,11,8)
                                    };

    KDTree tree = KDTree(nodes);

    tree.create_tree(0, nodes.size()-1, false);

    compare_node_trees(tree.get_nodes(), reference);
  }

  {
    // shuffle nodes
    const std::vector<Node> nodes = {Node(2,1,5),Node(3,11,4),Node(9,5,2),
                                     Node(4,4,4),Node(1,7,3),Node(5,5,8),
                                     Node(7,9,7),Node(8,2,3),Node(11,6,7),
                                     Node(10,2,9),Node(12,9,2.5),Node(14,8,6),
                                     Node(13,11.5,1),Node(6,10,2),Node(15,11,8)
                                    };

    KDTree tree = KDTree(nodes);

    tree.create_tree(0, nodes.size()-1, false);

    compare_node_trees(tree.get_nodes(), reference);
  }
}


TEST_CASE("create kd-tree: 6 nodes")
{
  const std::vector<Node> reference = {Node(4,4,2),Node(5,2,6),Node(2,7,6),
                                       Node(1,9,2),Node(6,12,3),Node(3,11,6)
                                      };
  {
    const std::vector<Node> nodes = {Node(1,9,2),Node(2,7,6),Node(3,11,6),
                                     Node(4,4,2),Node(5,2,6),Node(6,12,3)
                                    };

    KDTree tree = KDTree(nodes);

    tree.create_tree(0, nodes.size()-1, false);

    compare_node_trees(tree.get_nodes(), reference);
  }

  {
    // shuffle nodes
    const std::vector<Node> nodes = {Node(5,2,6),Node(2,7,6),Node(1,9,2),
                                     Node(6,12,3),Node(4,4,2),Node(3,11,6)
                                    };

    KDTree tree = KDTree(nodes);

    tree.create_tree(0, nodes.size()-1, false);

    compare_node_trees(tree.get_nodes(), reference);
  }
}


TEST_CASE("create kd-tree: 10 nodes")
{
  const std::vector<Node> reference = {Node(3,2,1),Node(1,4,2),Node(4,1,4),
                                       Node(7,3,6),Node(0,7,3),Node(5,9,2),
                                       Node(8,11,1),Node(2,12,3),Node(9,10,5),
                                       Node(6,11,6)
                                      };
  {
    const std::vector<Node> nodes = {Node(0,7,3),Node(1,4,2),Node(2,12,3),
                                     Node(3,2,1),Node(4,1,4),Node(5,9,2),
                                     Node(7,3,6),Node(6,11,6),Node(9,10,5),
                                     Node(8,11,1)
                                    };

    KDTree tree = KDTree(nodes);

    tree.create_tree(0, nodes.size()-1, false);

    compare_node_trees(tree.get_nodes(), reference);

    // check some values
    using namespace WorldBuilder;
    std::vector<Point<2>> check_points = {Point<2>(7,2,CoordinateSystem::cartesian),
                                          Point<2>(8,2,CoordinateSystem::cartesian),
                                          Point<2>(9,2,CoordinateSystem::cartesian),
                                          Point<2>(8,2.5001,CoordinateSystem::cartesian),
                                          Point<2>(8,2.4999,CoordinateSystem::cartesian),
                                          Point<2>(7,3,CoordinateSystem::cartesian),
                                          Point<2>(11,1.5,CoordinateSystem::cartesian),
                                          Point<2>(10,2,CoordinateSystem::cartesian),
                                          Point<2>(12,2,CoordinateSystem::cartesian),
                                          Point<2>(11,0.5,CoordinateSystem::cartesian),
                                          Point<2>(9.25,1.75,CoordinateSystem::cartesian),
                                          Point<2>(11,3.5,CoordinateSystem::cartesian),
                                          Point<2>(10.5,4,CoordinateSystem::cartesian),
                                          Point<2>(9.5,4,CoordinateSystem::cartesian),
                                          Point<2>(9.5,5.5,CoordinateSystem::cartesian),
                                          Point<2>(10.75,5.5,CoordinateSystem::cartesian),
                                          Point<2>(11.75,5.5,CoordinateSystem::cartesian),
                                          Point<2>(5,1,CoordinateSystem::cartesian),
                                          Point<2>(3,1,CoordinateSystem::cartesian),
                                          Point<2>(1,1,CoordinateSystem::cartesian),
                                          Point<2>(2,3,CoordinateSystem::cartesian),
                                          Point<2>(0.5,3,CoordinateSystem::cartesian),
                                          Point<2>(0.5,5,CoordinateSystem::cartesian),
                                          Point<2>(2,4.5,CoordinateSystem::cartesian),
                                          Point<2>(3,5,CoordinateSystem::cartesian),
                                          Point<2>(3,7,CoordinateSystem::cartesian)
                                         };

    std::vector<size_t> index_results = {4,5,5,4,5,4,6,5,7,6,5,7,8,8,8,9,9,1,0,0,2,2,2,2,3,3};

    for (size_t i = 0; i < check_points.size(); ++i)
      {
        IndexDistance const index_distance = tree.find_closest_point(check_points[i]);

        CHECK(index_distance.index == index_results[i]);
      }

    std::vector<std::vector<size_t> > index_distances_results =
    {
      {5,7,8,4,3,2,1,0},{5,7,4,3,2,1,0},{6,5,7,4,3,2,1},{5,6,7,8,4,3,2,1,0},
      {5,6,7,8,4,3,2,1,0},{8,9,7,5,6,4,},{6,5,7,4,0,1},{6,5,7,4,3,2,1,0},
      {6,5,7,4,3,2,1,0},{6,5,7,4,0,1},{6,5,7,4,0,1,3,2},{9,8,7,6,5,4,3,2,1,0},
      {9,8,7,6,5,4,3,2,1,0},{8,9,7,6,5,4,3,2,1,0},{8,9,7,6,5,4,3,2,1,0},
      {9,8,7,6,5,4,3,2,1,0},{9,8,7,6,5,4,3,2,1,0},{0,1,3,2,4},{0,1,4},{0,1,4},
      {3,2,1,0,4},{2,3,1,0,4},{2,3,1,0,4},{3,2,1,0,4},{3,2,1,0,4},{3,2,1,0,4}
    };

    for (size_t i = 0; i < check_points.size(); ++i)
      {
        IndexDistances index_distances = tree.find_closest_points(check_points[i]);
        CHECK(index_distances.min_index == index_results[i]);

        for (size_t j = 0; j < index_distances.vector.size(); ++j)
          {
            CHECK(index_distances.vector[j].index == index_distances_results[i][j]);
          }
      }
  }

  {
    // shuffle nodes
    const std::vector<Node> nodes = {Node(3,2,1),Node(4,1,4),Node(5,9,2),
                                     Node(0,7,3),Node(6,11,6),Node(2,12,3),
                                     Node(7,3,6),Node(9,10,5),Node(1,4,2),
                                     Node(8,11,1)
                                    };

    KDTree tree = KDTree(nodes);

    tree.create_tree(0, nodes.size()-1, false);

    compare_node_trees(tree.get_nodes(), reference);
  }
}