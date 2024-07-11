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
#include "world_builder/assert.h"
#include "world_builder/objects/surface.h"

#include <iostream>

#include "delaunator-cpp/delaunator.hpp"


namespace WorldBuilder
{
  namespace Objects
  {

    namespace
    {
      /**
       * Test whether a point is in a triangle. If that is the case is stores the interpolated
       * value of the triangle into `interpolated_value` and returns true.
       */
      bool in_triangle(const std::array<std::array<double,3>,3> &points,
                       const std::array<double,8> &precomputed,
                       const Point<2> &check_point,
                       double &interpolate_value,
                       double &interpolator_s,
                       double &interpolator_t)
      {
        double factor = 1e4;
        // based on https://stackoverflow.com/questions/2049582/how-to-determine-if-a-point-is-in-a-2d-triangle
        // compute s, t and area
        const double s_no_area = -(precomputed[0] + precomputed[1]*check_point[0] + precomputed[2]*check_point[1]);
        const double t_no_area = -(precomputed[3] + precomputed[4]*check_point[0] + precomputed[5]*check_point[1]);

        if (s_no_area >= -factor*std::numeric_limits<double>::epsilon() && t_no_area >= -factor*std::numeric_limits<double>::epsilon() && s_no_area+t_no_area-precomputed[6]<=precomputed[6]*factor*std::numeric_limits<double>::epsilon())
          {
            // point is in this triangle
            interpolator_s = precomputed[7]*s_no_area;
            interpolator_t = precomputed[7]*t_no_area;
            interpolate_value = points[0][2]*(1-interpolator_s-interpolator_t)+points[1][2]*interpolator_s+points[2][2]*interpolator_t;
            return true;
          }
        return false;
      }
    }

    Surface::Surface()
      = default;


    Surface::Surface(std::pair<std::vector<double>,std::vector<double>> values_at_points)
    {
      // first find the min and max. This need always to be done;
      WBAssertThrow(!values_at_points.first.empty(), "internal error: no values in values_at_points.first.");
      minimum = values_at_points.first[0];
      maximum = values_at_points.first[0];

      for (const auto &value: values_at_points.first)
        {
          if (value < minimum)
            {
              minimum = value;
            }
          if (value > maximum)
            {
              maximum = value;
            }
        }

      // if there are points defined there is not a constant value,
      // if there are no points defined there is a constant value.
      if (!values_at_points.second.empty())
        {
          constant_value = false;

          // Next make the triangulation and the kd-tree.
          // Feed the delaunator
          delaunator::Delaunator triangulation(values_at_points.second);

          // now loop over all the triangles and add them to the correct vectors
          triangles.resize(triangulation.triangles.size()/3);
          std::vector<KDTree::Node> node_list;
          for (std::size_t i = 0; i < triangulation.triangles.size(); i+=3)
            {
              delaunator::Delaunator &t = triangulation;
              // 1. We need to create a triangle list where each entry contains the coordinates
              //    and value of the points of the triangle:
              //    [[[x0,y0,v0],[x1,y1,v1],[x2,y2,v2]],[[x1,y1,v1],[x0,y0,v0],[x3,y3,v3]],...]
              // 2. we need to compute the centeroid locations and add those to the KD-tree
              //    with the index stored on the triangle list
              triangles[i/3][0][0] = t.coords[2 * t.triangles[i]]; // tx0
              triangles[i/3][0][1] = t.coords[2 * t.triangles[i]+1]; // ty0
              triangles[i/3][0][2] = values_at_points.first[t.triangles[i]]; // v0
              triangles[i/3][1][0] = t.coords[2 * t.triangles[i+1]]; // tx1
              triangles[i/3][1][1] = t.coords[2 * t.triangles[i+1]+1]; // ty1
              triangles[i/3][1][2] = values_at_points.first[t.triangles[i+1]]; // v1
              triangles[i/3][2][0] = t.coords[2 * t.triangles[i+2]]; // tx2
              triangles[i/3][2][1] = t.coords[2 * t.triangles[i+2]+1]; // ty2
              triangles[i/3][2][2] = values_at_points.first[t.triangles[i+2]]; // v2

              node_list.emplace_back(i/3,
                                     (t.coords[2*t.triangles[i]]+t.coords[2*t.triangles[i+1]]+t.coords[2*t.triangles[i+2]])/3.,
                                     (t.coords[2*t.triangles[i]+1]+t.coords[2*t.triangles[i+1]+1]+t.coords[2*t.triangles[i+2]+1])/3.);

            }
          // now feed and create the KD-tree
          tree = KDTree::KDTree(node_list);
          tree.create_tree(0, node_list.size()-1, false);


          // precompute values for in_triangle
          in_triangle_precomputed.resize(tree.get_nodes().size());
          for (size_t iii = 0; iii < tree.get_nodes().size(); iii++)
            {
              in_triangle_precomputed[iii][0] = triangles[iii][0][1]*triangles[iii][2][0] - triangles[iii][0][0]*triangles[iii][2][1];
              in_triangle_precomputed[iii][1] = triangles[iii][2][1] - triangles[iii][0][1];
              in_triangle_precomputed[iii][2] = triangles[iii][0][0] - triangles[iii][2][0];
              in_triangle_precomputed[iii][3] = triangles[iii][0][0]*triangles[iii][1][1] - triangles[iii][0][1]*triangles[iii][1][0];
              in_triangle_precomputed[iii][4] = triangles[iii][0][1] - triangles[iii][1][1];
              in_triangle_precomputed[iii][5] = triangles[iii][1][0] - triangles[iii][0][0];
              in_triangle_precomputed[iii][6] = -(-triangles[iii][1][1]*triangles[iii][2][0] + triangles[iii][0][1]*(-triangles[iii][1][0] + triangles[iii][2][0]) + triangles[iii][0][0]*(triangles[iii][1][1] - triangles[iii][2][1]) + triangles[iii][1][0]*triangles[iii][2][1]);
              in_triangle_precomputed[iii][7] = 1./in_triangle_precomputed[iii][6];
            }
        }
      else
        {
          constant_value = true;
        }

    }

    SurfaceValueInfo
    Surface::local_value(const Point<2> &check_point) const
    {
      if (constant_value)
        {
          // just min and max are the same since it is constant. Just return min.
          return minimum;
        }
      // first find the closest centeroids
      const KDTree::IndexDistances index_distances = tree.find_closest_points(check_point);

      // try triangle of the closest centroid
      double interpolated_value = 0;
      double interpolator_s = NaN::DQNAN;
      double interpolator_t = NaN::DQNAN;

      if (in_triangle(triangles[tree.get_nodes()[index_distances.min_index].index],in_triangle_precomputed[tree.get_nodes()[index_distances.min_index].index],check_point,interpolated_value,interpolator_s,interpolator_t))
        {
          return SurfaceValueInfo {tree.get_nodes()[index_distances.min_index].index,interpolated_value,interpolator_s,interpolator_t};
        }

      Point<2> other_point = check_point;
      KDTree::IndexDistances index_distances_other;
      const bool spherical = check_point.get_coordinate_system() == CoordinateSystem::spherical;
      if (spherical)
        {
          other_point[0] += check_point[0] < 0 ? 2.0 * WorldBuilder::Consts::PI : -2.0 * WorldBuilder::Consts::PI;
          index_distances_other = tree.find_closest_points(other_point);

          if (in_triangle(triangles[tree.get_nodes()[index_distances_other.min_index].index],in_triangle_precomputed[tree.get_nodes()[index_distances_other.min_index].index],other_point,interpolated_value,interpolator_s,interpolator_t))
            {
              return SurfaceValueInfo {tree.get_nodes()[index_distances_other.min_index].index,interpolated_value,interpolator_s,interpolator_t};
            }
        }


      {
        // if not found go to closets nodes
        // Todo: could remove the cosest node, because it was already tested. Could also sort based no distance.
        for (const auto &index_distance: index_distances.vector)
          {
            if (in_triangle(triangles[tree.get_nodes()[index_distance.index].index],in_triangle_precomputed[tree.get_nodes()[index_distance.index].index],check_point,interpolated_value,interpolator_s,interpolator_t))
              {
                return SurfaceValueInfo {tree.get_nodes()[index_distance.index].index,interpolated_value,interpolator_s,interpolator_t};
              }
          }
        if (spherical)
          {
            for (auto &index_distance: index_distances_other.vector)
              {
                if (in_triangle(triangles[tree.get_nodes()[index_distance.index].index],in_triangle_precomputed[tree.get_nodes()[index_distance.index].index],other_point,interpolated_value,interpolator_s,interpolator_t))
                  {
                    return SurfaceValueInfo {tree.get_nodes()[index_distance.index].index,interpolated_value,interpolator_s,interpolator_t};
                  }
              }
          }

        // if still not found, go through all nodes
        // Todo: Although this shouldonly very rearly happen, could remove already tested nodes.
        for (const auto &nodes: tree.get_nodes())
          {
            if (in_triangle(triangles[nodes.index],in_triangle_precomputed[nodes.index],check_point,interpolated_value,interpolator_s,interpolator_t))
              {
                return SurfaceValueInfo {nodes.index,interpolated_value,interpolator_s,interpolator_t};
              }
            if (spherical && in_triangle(triangles[nodes.index],in_triangle_precomputed[nodes.index],other_point,interpolated_value,interpolator_s,interpolator_t))
              {
                return SurfaceValueInfo {nodes.index,interpolated_value,interpolator_s,interpolator_t};
              }
          }
        WBAssertThrow(false, "Internal error: The requested point was not in any triangle. "
                      << "This could be due to rounding errors if the difference between the check point and triangle points are small, "
                      << "or you are requesting a point outside the boundaries defined by the additional points. The check point was "
                      << check_point[0] <<  ":" << check_point[1] << ".");
      }
      return 0;
    }
  } // namespace Objects
} // namespace WorldBuilder

