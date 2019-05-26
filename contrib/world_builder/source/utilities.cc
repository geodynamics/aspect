/*
  Copyright (C) 2018 by the authors of the World Builder code.

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

#include <iostream>
#include <iomanip>
#include <algorithm>

#include <world_builder/assert.h>
#include <world_builder/coordinate_systems/interface.h>
#include <world_builder/nan.h>
#include <world_builder/utilities.h>


namespace WorldBuilder
{
  namespace Utilities
  {
    bool
    polygon_contains_point(const std::vector<Point<2> > &point_list,
                           const Point<2> &point)
    {
      if (point.get_coordinate_system() == CoordinateSystem::spherical)
        {
          Point<2> other_point = point;
          other_point[0] += point[0] < 0 ? 2.0 * const_pi : -2.0 * const_pi;

          return (polygon_contains_point_implementation(point_list, point) ||
                  polygon_contains_point_implementation(point_list, other_point));
        }
      else
        {
          return polygon_contains_point_implementation(point_list, point);
        }
    }

    bool
    polygon_contains_point_implementation(const std::vector<Point<2> > &point_list,
                                          const Point<2> &point)
    {
      /**
       * This code has been based on http://geomalgorithms.com/a03-_inclusion.html,
       * and therefore requires the following copyright notice:
       *
       * Copyright 2000 softSurfer, 2012 Dan Sunday
       * This code may be freely used and modified for any purpose
       * providing that this copyright notice is included with it.
       * SoftSurfer makes no warranty for this code, and cannot be held
       * liable for any real or imagined damage resulting from its use.
       * Users of this code must verify correctness for their application.
       *
       * The main functional difference between the original code and this
       * code is that all the boundaries are considered to be inside the
       * polygon. One should of course realize that with floating point
       * arithmetic no guarantees can be made for the borders, but for
       * exact arithmetic this algorithm would work (also see polygon
       * in point test).
       */
      int pointNo = point_list.size();
      int    wn = 0;    // the  winding number counter
      int   j=pointNo-1;

      // loop through all edges of the polygon
      for (int i=0; i<pointNo; i++)
        {
          // edge from V[i] to  V[i+1]
          if (point_list[j][1] <= point[1])
            {
              // start y <= P.y
              if (point_list[i][1] >= point[1])      // an upward crossing
                {
                  const double is_left = (point_list[i][0] - point_list[j][0]) * (point[1] - point_list[j][1])
                                         - (point[0] -  point_list[j][0]) * (point_list[i][1] - point_list[j][1]);

                  if ( is_left > 0 && point_list[i][1] > point[1])
                    {
                      // P left of  edge
                      ++wn;            // have  a valid up intersect
                    }
                  else if ( is_left == 0)
                    {
                      // The point is exactly on the infinite line.
                      // determine if it is on the segment
                      const double dot_product = (point - point_list[j])*(point_list[i] - point_list[j]);

                      if (dot_product >= 0)
                        {
                          const double squaredlength = (point_list[i] - point_list[j]).norm_square();

                          if (dot_product <= squaredlength)
                            {
                              return true;
                            }
                        }
                    }
                }
            }
          else
            {
              // start y > P.y (no test needed)
              if (point_list[i][1]  <= point[1])     // a downward crossing
                {
                  const double is_left = (point_list[i][0] - point_list[j][0]) * (point[1] - point_list[j][1])
                                         - (point[0] -  point_list[j][0]) * (point_list[i][1] - point_list[j][1]);

                  if ( is_left < 0)
                    {
                      // P right of  edge
                      --wn;            // have  a valid down intersect
                    }
                  else if ( is_left == 0)
                    {
                      // This code is to make sure that the boundaries are included in the polygon.
                      // The point is exactly on the infinite line.
                      // determine if it is on the segment
                      const double dot_product = (point - point_list[j])*(point_list[i] - point_list[j]);

                      if (dot_product >= 0)
                        {
                          const double squaredlength = (point_list[i] - point_list[j]).norm_square();

                          if (dot_product <= squaredlength)
                            {
                              return true;
                            }
                        }
                    }
                }
            }
          j=i;
        }

      return (wn != 0);
    }

    double
    signed_distance_to_polygon(const std::vector<Point<2> > &point_list,
                               const Point<2> &point)
    {
      // If the point lies outside polygon, we give it a negative sign,
      // inside a positive sign.
      const double sign = polygon_contains_point(point_list, point) ? 1.0 : -1.0;

      /**
       * This code is based on http://geomalgorithms.com/a02-_lines.html#Distance-to-Infinite-Line,
       * and therefore requires the following copyright notice:
       *
       * Copyright 2000 softSurfer, 2012 Dan Sunday
       * This code may be freely used and modified for any purpose
       * providing that this copyright notice is included with it.
       * SoftSurfer makes no warranty for this code, and cannot be held
       * liable for any real or imagined damage resulting from its use.
       * Users of this code must verify correctness for their application.
       *
       */

      const unsigned int n_poly_points = point_list.size();
      WBAssertThrow(n_poly_points >= 3, "Not enough polygon points were specified.");

      // Initialize a vector of distances for each point of the polygon with a very large distance
      std::vector<double> distances(n_poly_points, 1e23);

      // Create another polygon but with all points shifted 1 position to the right
      std::vector<Point<2> > shifted_point_list(n_poly_points, Point<2>(point.get_coordinate_system()));
      shifted_point_list[0] = point_list[n_poly_points-1];

      for (unsigned int i = 0; i < n_poly_points-1; ++i)
        shifted_point_list[i+1] = point_list[i];

      for (unsigned int i = 0; i < n_poly_points; ++i)
        {
          // Create vector along the polygon line segment
          Point<2> vector_segment = shifted_point_list[i] - point_list[i];
          // Create vector from point to the second segment point
          Point<2> vector_point_segment = point - point_list[i];

          // Compute dot products to get angles
          const double c1 = vector_point_segment * vector_segment;
          const double c2 = vector_segment * vector_segment;

          // point lies closer to not-shifted polygon point, but perpendicular base line lies outside segment
          if (c1 <= 0.0)
            distances[i] = (Point<2> (point_list[i] - point)).norm();
          // point lies closer to shifted polygon point, but perpendicular base line lies outside segment
          else if (c2 <= c1)
            distances[i] = (Point<2> (shifted_point_list[i] - point)).norm();
          // perpendicular base line lies on segment
          else
            {
              const Point<2> point_on_segment = point_list[i] + (c1/c2) * vector_segment;
              distances[i] = (Point<2> (point - point_on_segment)).norm();
            }
        }

      // Return the minimum of the distances of the point to all polygon segments
      return *std::min_element(distances.begin(),distances.end()) * sign;
    }


    NaturalCoordinate::NaturalCoordinate(const std::array<double,3> &position,
                                         const CoordinateSystems::Interface &coordinate_system_)
    {
      coordinate_system = coordinate_system_.natural_coordinate_system();
      coordinates = coordinate_system_.cartesian_to_natural_coordinates(position);
    }

    // todo, should be possible to make this without the interface, since the Point knows the coord system.
    NaturalCoordinate::NaturalCoordinate(const Point<3> &position,
                                         const CoordinateSystems::Interface &coordinate_system_)
    {
      coordinate_system = coordinate_system_.natural_coordinate_system();
      coordinates = coordinate_system_.cartesian_to_natural_coordinates(position.get_array());
    }

    const std::array<double,3> &NaturalCoordinate::get_coordinates()
    {
      return coordinates;
    }


    /*std::array<double,1> NaturalCoordinate::get_surface_coordinates() const
    {
      std::array<double,1> coordinate;

      switch (coordinate_system)
        {
          case Coordinates::CoordinateSystem::cartesian:
            coordinate[0] = coordinates[0];
            break;

          case Coordinates::CoordinateSystem::spherical:
            coordinate[0] = coordinates[1];
            break;

          default:
            coordinate[0] = 0;
            WBAssert (false, ExcNotImplemented());
            break;
        }

      return coordinate;
    }*/


    const std::array<double,2> NaturalCoordinate::get_surface_coordinates() const
    {
      std::array<double,2> coordinate;

      switch (coordinate_system)
        {
          case CoordinateSystem::cartesian:
            coordinate[0] = coordinates[0];
            coordinate[1] = coordinates[1];
            break;

          case CoordinateSystem::spherical:
            coordinate[0] = coordinates[1];
            coordinate[1] = coordinates[2];
            break;

          default:
            WBAssert (false, "Coordinate system not implemented.");
        }

      return coordinate;
    }


    CoordinateSystem
    NaturalCoordinate::get_coordinate_system() const
    {
      return coordinate_system;
    }


    double NaturalCoordinate::get_depth_coordinate() const
    {
      switch (coordinate_system)
        {
          case CoordinateSystem::cartesian:
            return coordinates[2];

          case CoordinateSystem::spherical:
            return coordinates[0];

          default:
            WBAssert (false, "Coordinate system not implemented.");
        }

      return 0;
    }


    std::array<double,3>
    cartesian_to_spherical_coordinates(const Point<3> &position)
    {
      std::array<double,3> scoord;

      scoord[0] = position.norm(); // R
      scoord[1] = std::atan2(position[1],position[0]); // Phi
      //if (scoord[1] < 0.0)
      //scoord[1] += 2.0*const_pi; // correct phi to [0,2*pi]

      if (scoord[0] > std::numeric_limits<double>::min())
        scoord[2] = 0.5 * const_pi - std::acos(position[2]/scoord[0]);
      else
        scoord[2] = 0.0;

      return scoord;
    }

    Point<3>
    spherical_to_cartesian_coordinates(const std::array<double,3> &scoord)
    {
      Point<3> ccoord(cartesian);

      ccoord[0] = scoord[0] * std::sin(0.5 * const_pi - scoord[2]) * std::cos(scoord[1]); // X
      ccoord[1] = scoord[0] * std::sin(0.5 * const_pi - scoord[2]) * std::sin(scoord[1]); // Y
      ccoord[2] = scoord[0] * std::cos(0.5 * const_pi - scoord[2]); // Z


      return ccoord;
    }



    CoordinateSystem
    string_to_coordinate_system(const std::string &coordinate_system)
    {
      if (coordinate_system == "cartesian")
        return CoordinateSystem::cartesian;
      else if (coordinate_system == "spherical")
        return CoordinateSystem::spherical;
      else
        WBAssertThrow(false, "Coordinate system not implemented.");

      return invalid;
    }


    template<int dim>
    const std::array<double,dim>
    convert_point_to_array(const Point<dim> &point_)
    {
      std::array<double,dim> array;
      for (unsigned int i = 0; i < dim; ++i)
        array[i] = point_[i];
      return array;
    }

    double
    string_to_double(const std::string &string)
    {
      // trim whitespace on either side of the text if necessary
      std::string s = string;
      while ((s.size() > 0) && (s[0] == ' '))
        s.erase(s.begin());
      while ((s.size() > 0) && (s[s.size() - 1] == ' '))
        s.erase(s.end() - 1);

      std::istringstream i(s);
      double d;
      char c;
      if (!(i >> d) || i.get(c))
        WBAssertThrow(false, "Could not convert \"" + s + "\" to a double.");

      return d;
    }

    int
    string_to_int(const std::string &string)
    {
      // trim whitespace on either side of the text if necessary
      std::string s = string;
      while ((s.size() > 0) && (s[0] == ' '))
        s.erase(s.begin());
      while ((s.size() > 0) && (s[s.size() - 1] == ' '))
        s.erase(s.end() - 1);

      std::istringstream i(s);
      int d;
      char c;
      if (!(i >> d) || i.get(c))
        WBAssertThrow(false, "Could not convert \"" + s + "\" to an int.");

      return d;
    }


    unsigned int
    string_to_unsigned_int(const std::string &string)
    {
      // trim whitespace on either side of the text if necessary
      std::string s = string;
      while ((s.size() > 0) && (s[0] == ' '))
        s.erase(s.begin());
      while ((s.size() > 0) && (s[s.size() - 1] == ' '))
        s.erase(s.end() - 1);


      std::istringstream i(s);
      unsigned int d;
      char c;
      if (!(i >> d) || i.get(c))
        WBAssertThrow(false, "Could not convert \"" + s + "\" to an unsigned int.");

      return d;
    }


    Point<3>
    cross_product(const Point<3> &a, const Point<3> &b)
    {
      WBAssert(a.get_coordinate_system() == b.get_coordinate_system(), "Trying to do a cross product of points of a different coordinate system.");
      const double x = a[1] * b[2] - b[1] * a[2];
      const double y = a[2] * b[0] - b[2] * a[0];
      const double z = a[0] * b[1] - b[0] * a[1];
      return Point<3>(x,y,z,a.get_coordinate_system());
    }

    std::map<std::string,double>
    distance_point_from_curved_planes(const Point<3> &check_point, // cartesian point in spherical system
                                      const Point<2> &reference_point, // in (rad) spherical coordinates in spherical system
                                      const std::vector<Point<2> > &point_list, // in  (rad) spherical coordinates in spherical system
                                      const std::vector<std::vector<double> > &plane_segment_lengths,
                                      const std::vector<std::vector<Point<2> > > &plane_segment_angles,
                                      const double start_radius,
                                      const std::unique_ptr<CoordinateSystems::Interface> &coordinate_system,
                                      const bool only_positive,
                                      std::vector<double> global_x_list)
    {
      // TODO: Assert that point_list, plane_segment_angles and plane_segment_lenghts have the same size.
      /*WBAssert(point_list.size() == plane_segment_lengths.size(),
               "Internal error: The size of point_list (" << point_list.size()
               << ") and plane_segment_lengths (" << plane_segment_lengths.size() << ") are different.");
      WBAssert(point_list.size() == plane_segment_angles.size(),
               "Internal error: The size of point_list (" << point_list.size()
               << ") and plane_segment_angles (" << plane_segment_angles.size() << ") are different.");
      WBAssert(point_list.size() == plane_segment_angles.size(),
               "Internal error: The size of point_list (" << point_list.size()
               << ") and global_x_list (" << global_x_list.size() << ") are different.");*/
      if (global_x_list.size() == 0)
        {
          // fill it
          global_x_list.resize(point_list.size());
          for (unsigned int i = 0; i < point_list.size(); ++i)
            global_x_list[i] = i;
        }
      WBAssertThrow(global_x_list.size() == point_list.size(), "The given global_x_list doesn't have "
                    "the same size as the point list. This is required.");

      double distance = INFINITY;
      double new_distance = INFINITY;
      double along_plane_distance = INFINITY;
      double new_along_plane_distance  = INFINITY;

      const CoordinateSystem natural_coordinate_system = coordinate_system->natural_coordinate_system();
      const bool bool_cartesian = natural_coordinate_system == cartesian;

      const Point<3> check_point_natural(coordinate_system->cartesian_to_natural_coordinates(check_point.get_array()),natural_coordinate_system);
      const Point<3> check_point_surface(bool_cartesian ? check_point_natural[0] : start_radius,
                                         check_point_natural[1],
                                         bool_cartesian ? start_radius           : check_point_natural[2],
                                         natural_coordinate_system);
      const Point<2> check_point_surface_2d(bool_cartesian ? check_point_natural[0] : check_point_natural[1],
                                            bool_cartesian ? check_point_natural[1] : check_point_natural[2],
                                            natural_coordinate_system);

      // The section which is checked.
      double section = 0.0;

      // The 'horizontal' fraction between the points at the surface.
      double section_fraction = 0.0;

      // What segment the point on the line is in.
      unsigned int segment = 0;

      // The 'vertical' fraction, indicates how far in the current segment the
      // point on the line is.
      double segment_fraction = 0.0;
      double total_average_angle = 0.0;

      const DepthMethod depth_method = coordinate_system->depth_method();
      WBAssertThrow(depth_method == DepthMethod::none
                    || depth_method == DepthMethod::angle_at_starting_point_with_surface
                    || depth_method == DepthMethod::angle_at_begin_segment_with_surface,
                    "Only the depth methods none, angle_at_starting_point_with_surface and "
                    "angle_at_begin_segment_with_surface are implemented");

      // loop over all the planes to find out which one is closest to the point.

      for (unsigned int i_section=0; i_section < point_list.size()-1; ++i_section)
        {
          const unsigned int current_section = i_section;
          const unsigned int next_section = i_section+1;
          // translate to orignal coordinates current and next section
          const unsigned int original_current_section = (unsigned int)std::floor(global_x_list[i_section]);
          const unsigned int original_next_section = original_current_section + 1;
          // see on what side the line P1P2 reference point is. This is based on the determinant
          const double reference_on_side_of_line = (point_list[next_section][0] - point_list[current_section][0])
                                                   * (reference_point[1] - point_list[current_section][1])
                                                   - (point_list[next_section][1] - point_list[current_section][1])
                                                   * (reference_point[0] - point_list[current_section][0])
                                                   < 0 ? 1 : -1;



          const Point<2> P1(point_list[current_section]);

          const Point<2> P2(point_list[next_section]);

          const Point<2> P1P2 = P2 - P1;
          const Point<2> P1PC = check_point_surface_2d - P1;


          // Compute the closest point on the line P1 to P2 from the check
          // point at the surface. We do this in natural coordinates on
          // purpose, because in spherical coordinates it is more accurate.
          Point<2> closest_point_on_line_2d = P1 + ((P1PC * P1P2) / (P1P2 * P1P2)) * P1P2;


          // compute what fraction of the distance between P1 and P2 the
          // closest point lies.
          const Point<2> P1CPL = closest_point_on_line_2d - P1;

          // This determines where the check point is between the coordinates
          // in the coordinate list.
          const double fraction_CPL_P1P2_strict = (P1CPL * P1P2 <= 0 ? -1.0 : 1.0)
                                                  * (1 - (P1P2.norm() - P1CPL.norm()) / P1P2.norm());

          // If the point on the line does not lay between point P1 and P2
          // then ignore it. Otherwise continue.
          if (fraction_CPL_P1P2_strict >= 0 && fraction_CPL_P1P2_strict <= 1.0)
            {
              // now figure out where the point is in relation with the user
              // defined coordinates
              const double fraction_CPL_P1P2 = global_x_list[i_section] - (int)global_x_list[i_section]
                                               + (global_x_list[i_section+1]-global_x_list[i_section]) * fraction_CPL_P1P2_strict;

              const Point<2> unit_normal_to_plane_spherical = P1P2 / P1P2.norm();
              const Point<2> closest_point_on_line_plus_normal_to_plane_spherical = closest_point_on_line_2d + 1e-8 * (closest_point_on_line_2d.norm() > 1.0 ? closest_point_on_line_2d.norm() : 1.0) * unit_normal_to_plane_spherical;

              WBAssert(closest_point_on_line_plus_normal_to_plane_spherical.norm() != 0.0,
                       "Internal error: The norm of variable 'closest_point_on_line_plus_normal_to_plane_spherical' "
                       "is  zero, while this may not happen.");

              // We now need 3d points from this point on, so make them.
              // The order of a Cartesian coordinate is x,y,z and the order of
              // a spherical coordinate it radius, long, lat (in rad).
              const Point<3> closest_point_on_line_surface(bool_cartesian ? closest_point_on_line_2d[0] : start_radius,
                                                           bool_cartesian ? closest_point_on_line_2d[1] : closest_point_on_line_2d[0],
                                                           bool_cartesian ? start_radius : closest_point_on_line_2d[1],
                                                           natural_coordinate_system);

              Point<3> closest_point_on_line_bottom = closest_point_on_line_surface;
              closest_point_on_line_bottom[bool_cartesian ? 2 : 0] = 0;

              const Point<3> closest_point_on_line_plus_normal_to_plane_surface_spherical(bool_cartesian ? closest_point_on_line_plus_normal_to_plane_spherical[0] : start_radius,
                                                                                          bool_cartesian ? closest_point_on_line_plus_normal_to_plane_spherical[1] : closest_point_on_line_plus_normal_to_plane_spherical[0],
                                                                                          bool_cartesian ? start_radius : closest_point_on_line_plus_normal_to_plane_spherical[1],
                                                                                          natural_coordinate_system);

              // Now that we have both the check point and the
              // closest_point_on_line, we need to push them to cartesian.
              const Point<3> check_point_cartesian(check_point);
              const Point<3> check_point_surface_cartesian(coordinate_system->natural_to_cartesian_coordinates(check_point_surface.get_array()),cartesian);
              const Point<3> closest_point_on_line_cartesian(coordinate_system->natural_to_cartesian_coordinates(closest_point_on_line_surface.get_array()),cartesian);
              const Point<3> closest_point_on_line_bottom_cartesian(coordinate_system->natural_to_cartesian_coordinates(closest_point_on_line_bottom.get_array()),cartesian);
              const Point<3> closest_point_on_line_plus_normal_to_plane_cartesian(coordinate_system->natural_to_cartesian_coordinates(closest_point_on_line_plus_normal_to_plane_surface_spherical.get_array()),cartesian);


              // if the two points are the same, we don't need to search any further
              if (std::fabs((check_point_cartesian - closest_point_on_line_cartesian).norm()) < 2e-14)
                {
                  distance = 0.0;
                  along_plane_distance = 0.0;
                  section = current_section;
                  section_fraction = fraction_CPL_P1P2;
                  segment = 0;
                  segment_fraction = 0.0;
                  total_average_angle = plane_segment_angles[original_current_section][0][0]
                                        + fraction_CPL_P1P2 * (plane_segment_angles[original_next_section][0][0]
                                                               - plane_segment_angles[original_current_section][0][0]);
                  break;
                }

              Point<3> normal_to_plane = closest_point_on_line_plus_normal_to_plane_cartesian - closest_point_on_line_cartesian;
              normal_to_plane = normal_to_plane / normal_to_plane.norm();

              // The y-axis is from the bottom/center to the closest_point_on_line,
              // the x-axis is 90 degrees rotated from that, so we rotate around
              // the line P1P2.
              // Todo: Assert that the norm of the axis are not equal to zero.
              Point<3> y_axis = closest_point_on_line_cartesian - closest_point_on_line_bottom_cartesian;

              WBAssert(y_axis.norm() != 0,
                       "Internal error: The y_axis.norm() is zero. Y_axis is " << y_axis[0] << ":" << y_axis[1] << ":" << y_axis[2]
                       << ". closest_point_on_line_cartesian = " << closest_point_on_line_cartesian[0] << ":" << closest_point_on_line_cartesian[1] << ":" << closest_point_on_line_cartesian[2]
                       << ", closest_point_on_line_bottom_cartesian = " << closest_point_on_line_bottom_cartesian[0] << ":" << closest_point_on_line_bottom_cartesian[1] << ":" << closest_point_on_line_bottom_cartesian[2]);

              WBAssert(!std::isnan(y_axis[0]),
                       "Internal error: The y_axis variable is not a number: " << y_axis[0]);
              WBAssert(!std::isnan(y_axis[1]),
                       "Internal error: The y_axis variable is not a number: " << y_axis[1]);
              WBAssert(!std::isnan(y_axis[2]),
                       "Internal error: The y_axis variable is not a number: " << y_axis[2]);


              y_axis = y_axis / y_axis.norm();


              WBAssert(!std::isnan(y_axis[0]),
                       "Internal error: The y_axis variable is not a number: " << y_axis[0]);
              WBAssert(!std::isnan(y_axis[1]),
                       "Internal error: The y_axis variable is not a number: " << y_axis[1]);
              WBAssert(!std::isnan(y_axis[2]),
                       "Internal error: The y_axis variable is not a number: " << y_axis[2]);


              // shorthand notation for computing the x_axis
              double vx = y_axis[0];
              double vy = y_axis[1];
              double vz = y_axis[2];
              double ux = normal_to_plane[0];
              double uy = normal_to_plane[1];
              double uz = normal_to_plane[2];

              Point<3> x_axis(ux*ux*vx + ux*uy*vy - uz*vy + uy*uz*vz + uy*vz,
                              uy*ux*vx + uz*vx + uy*uy*vy + uy*uz*vz - ux*vz,
                              uz*ux*vx - uy*vx + uz*uy*vy + ux*vy + uz*uz*vz,
                              cartesian);

              WBAssert(!std::isnan(x_axis[0]),
                       "Internal error: The x_axis variable is not a number: " << x_axis[0]);
              WBAssert(!std::isnan(x_axis[1]),
                       "Internal error: The x_axis variable is not a number: " << x_axis[1]);
              WBAssert(!std::isnan(x_axis[2]),
                       "Internal error: The x_axis variable is not a number: " << x_axis[2]);

              x_axis = x_axis *(reference_on_side_of_line / x_axis.norm());


              WBAssert(!std::isnan(x_axis[0]),
                       "Internal error: The x_axis variable is not a number: " << x_axis[0]);
              WBAssert(!std::isnan(x_axis[1]),
                       "Internal error: The x_axis variable is not a number: " << x_axis[1]);
              WBAssert(!std::isnan(x_axis[2]),
                       "Internal error: The x_axis variable is not a number: " << x_axis[2]);

              Point<2> check_point_2d(x_axis * (check_point_cartesian - closest_point_on_line_bottom_cartesian),
                                      y_axis * (check_point_cartesian - closest_point_on_line_bottom_cartesian),
                                      cartesian);


              Point<2> begin_segment(x_axis * (closest_point_on_line_cartesian - closest_point_on_line_bottom_cartesian),
                                     y_axis * (closest_point_on_line_cartesian - closest_point_on_line_bottom_cartesian),
                                     cartesian);


              WBAssert(!std::isnan(check_point_2d[0]),
                       "Internal error: The check_point_2d variable is not a number: " << check_point_2d[0]);
              WBAssert(!std::isnan(check_point_2d[1]),
                       "Internal error: The check_point_2d variable is not a number: " << check_point_2d[1]);


              WBAssert(!std::isnan(begin_segment[0]),
                       "Internal error: The begin_segment variable is not a number: " << begin_segment[0]);
              WBAssert(!std::isnan(begin_segment[1]),
                       "Internal error: The begin_segment variable is not a number: " << begin_segment[1]);

              Point<2> end_segment = begin_segment;


              double total_length = 0.0;
              double add_angle = 0.0;
              double average_angle = 0.0;
              for (unsigned int i_segment = 0; i_segment < plane_segment_lengths[original_current_section].size(); i_segment++)
                {
                  const unsigned int current_segment = i_segment;

                  // compute the angle between the the previous begin and end if
                  // the depth method is angle_at_begin_segment_with_surface.
                  if (i_segment != 0 && depth_method == DepthMethod::angle_at_begin_segment_with_surface)
                    {
                      const double add_angle_inner = (begin_segment * end_segment) / (begin_segment.norm() * end_segment.norm());

                      WBAssert(!std::isnan(add_angle_inner),
                               "Internal error: The add_angle_inner variable is not a number: " << add_angle_inner
                               << ". Variables: begin_segment = " << begin_segment[0] << ":" << begin_segment[1]
                               << ", end_segment = " << end_segment[0] << ":" << end_segment[1]
                               << ", begin_segment * end_segment / (begin_segment.norm() * end_segment.norm()) = "
                               << std::setprecision(32) << begin_segment * end_segment / (begin_segment.norm() * end_segment.norm())
                               << ".");

                      // there could be round of error problems here is the inner part is close to one
                      WBAssert(add_angle_inner >= 0 && add_angle_inner <= 1,
                               "Internal error: The variable add_angle_inner is smaller than zero or larger then one,"
                               "which causes the std::acos to return nan. If it is only a little bit larger then one, "
                               "this is probably caused by that begin and end segment are the same and round off error. "
                               "The value of add_angle_inner = " << add_angle_inner);

                      add_angle += std::acos(add_angle_inner);

                      WBAssert(!std::isnan(add_angle),
                               "Internal error: The add_angle variable is not a number: " << add_angle
                               << ". Variables: begin_segment = " << begin_segment[0] << ":" << begin_segment[1]
                               << ", end_segment = " << end_segment[0] << ":" << end_segment[1]
                               << ", begin_segment * end_segment / (begin_segment.norm() * end_segment.norm()) = "
                               << std::setprecision(32) << begin_segment * end_segment / (begin_segment.norm() * end_segment.norm())
                               << ", std::acos(begin_segment * end_segment / (begin_segment.norm() * end_segment.norm())) = "
                               << std::acos(begin_segment * end_segment / (begin_segment.norm() * end_segment.norm())));
                    }




                  begin_segment = end_segment;

                  WBAssert(!std::isnan(begin_segment[0]),
                           "Internal error: The begin_segment variable is not a number: " << begin_segment[0]);
                  WBAssert(!std::isnan(begin_segment[1]),
                           "Internal error: The begin_segment variable is not a number: " << begin_segment[1]);


                  // This interpolates different properties between P1 and P2 (the
                  // points of the plane at the surface)
                  const double degree_90_to_rad = 0.5 * const_pi;

                  WBAssert(plane_segment_angles.size() > original_next_section,
                           "Error: original_next_section = " << original_next_section
                           << ", and plane_segment_angles.size() = " << plane_segment_angles.size());


                  WBAssert(plane_segment_angles[original_next_section].size() > current_segment,
                           "Error: current_segment = "  << current_segment
                           << ", and current_segment.size() = " << plane_segment_angles[original_next_section].size());

                  /*std::cout << "plane_segment_angles = " << plane_segment_angles[original_current_section][current_segment][1] << std::endl;
                  std::cout << "fraction_CPL_P1P2 = " << fraction_CPL_P1P2 << std::endl;
                  std::cout << "original_next_section = " << original_next_section << std::endl;
                  std::cout << "plane_segment_angles.size() = " << plane_segment_angles.size() << std::endl;
                  std::cout << "current_segment = " << current_segment << std::endl;
                  std::cout << "plane_segment_angles[original_next_section].size() = " << plane_segment_angles[original_next_section].size() << std::endl;
                  std::cout << "add_angle = " << add_angle << std::endl;
                  */

                  const double interpolated_angle_top    = plane_segment_angles[original_current_section][current_segment][0]
                                                           + fraction_CPL_P1P2 * (plane_segment_angles[original_next_section][current_segment][0]
                                                                                  - plane_segment_angles[original_current_section][current_segment][0])
                                                           + add_angle;

                  const double interpolated_angle_bottom = plane_segment_angles[original_current_section][current_segment][1]
                                                           + fraction_CPL_P1P2 * (plane_segment_angles[original_next_section][current_segment][1]
                                                                                  - plane_segment_angles[original_current_section][current_segment][1])
                                                           + add_angle;


                  double interpolated_segment_length     = plane_segment_lengths[original_current_section][current_segment]
                                                           + fraction_CPL_P1P2 * (plane_segment_lengths[original_next_section][current_segment]
                                                                                  - plane_segment_lengths[original_current_section][current_segment]);
                  WBAssert(!std::isnan(interpolated_angle_top),
                           "Internal error: The interpolated_angle_top variable is not a number: " << interpolated_angle_top);

                  // We want to know where the end point of this segment is (and
                  // the start of the next segment). There are two cases which we
                  // will deal with separately. The first one is if the angle is
                  // constant. The second one is if the angle changes.
                  const double difference_in_angle_along_segment = interpolated_angle_top - interpolated_angle_bottom;

                  if (std::fabs(difference_in_angle_along_segment) < 1e-8)
                    {
                      // The angle is constant. It is easy find find the end of
                      // this segment and the distance.
                      if (interpolated_segment_length != 0)
                        {
                          end_segment[0] += interpolated_segment_length * std::sin(degree_90_to_rad - interpolated_angle_top);
                          end_segment[1] -= interpolated_segment_length * std::cos(degree_90_to_rad - interpolated_angle_top);

                          Point<2> begin_end_segment = end_segment - begin_segment;
                          Point<2> normal_2d_plane(-begin_end_segment[0],begin_end_segment[1], cartesian);
                          WBAssert(normal_2d_plane.norm() != 0, "Internal Error: normal_2d_plane.norm() is zero, which should not happen. "
                                   << "Extra info: begin_end_segment[0] = " << begin_end_segment[0]
                                   << ", begin_end_segment[1] = " << begin_end_segment[1]
                                   << ", end_segment: [" << end_segment[0] << "," << end_segment[1] << "]"
                                   << ", begin_segment: [" << begin_segment[0] << "," << begin_segment[1] << "]"
                                  );
                          normal_2d_plane /= normal_2d_plane.norm();

                          // Now find the distance of a point to this line.
                          // Based on http://geomalgorithms.com/a02-_lines.html.
                          const Point<2> BSP_ESP = end_segment - begin_segment;
                          const Point<2> BSP_CP = check_point_2d - begin_segment;

                          const double c1 = BSP_ESP * BSP_CP;
                          const double c2 = BSP_ESP * BSP_ESP;

                          if (c1 < 0 || c2 < c1)
                            {
                              new_distance = INFINITY;
                              new_along_plane_distance = INFINITY;
                            }
                          else
                            {
                              const Point<2> Pb = begin_segment + (c1/c2) * BSP_ESP;
                              const double side_of_line =  (begin_segment[0] - end_segment[0]) * (check_point_2d[1] - begin_segment[1])
                                                           - (begin_segment[1] - end_segment[1]) * (check_point_2d[0] - begin_segment[0])
                                                           < 0 ? -1.0 : 1.0;

                              new_distance = side_of_line * (check_point_2d - Pb).norm();
                              new_along_plane_distance = (begin_segment - Pb).norm();
                            }

                        }
                    }
                  else
                    {
                      // The angle is not constant. This means that we need to
                      // define a circle. First find the center of the circle.
                      const double radius_angle_circle = std::fabs(interpolated_segment_length/difference_in_angle_along_segment);

                      WBAssert(!std::isnan(radius_angle_circle),
                               "Internal error: The radius_angle_circle variable is not a number: " << radius_angle_circle
                               << ". interpolated_segment_length = " << interpolated_segment_length
                               << ", difference_in_angle_along_segment = " << difference_in_angle_along_segment);

                      const double cos_angle_top = std::cos(interpolated_angle_top);

                      WBAssert(!std::isnan(cos_angle_top),
                               "Internal error: The radius_angle_circle variable is not a number: " << cos_angle_top
                               << ". interpolated_angle_top = " << interpolated_angle_top);

                      Point<2> center_circle(cartesian);
                      if (std::fabs(interpolated_angle_top - 0.5 * const_pi) < 1e-8)
                        {
                          // if interpolated_angle_top is 90 degrees, the tan function
                          // is undefined (1/0). What we really want in this case is
                          // set the center to the correct location which is x = the x
                          //begin point + radius and y = the y begin point.
                          center_circle[0] = difference_in_angle_along_segment > 0 ? begin_segment[0] + radius_angle_circle : begin_segment[0] - radius_angle_circle;
                          center_circle[1] = begin_segment[1];
                        }
                      else if (std::fabs(interpolated_angle_top - 1.5 * const_pi) < 1e-8)
                        {
                          // if interpolated_angle_top is 270 degrees, the tan function
                          // is undefined (-1/0). What we really want in this case is
                          // set the center to the correct location which is x = the x
                          //begin point - radius and y = the y begin point.
                          center_circle[0] = difference_in_angle_along_segment > 0 ? begin_segment[0] - radius_angle_circle : begin_segment[0] + radius_angle_circle;
                          center_circle[1] = begin_segment[1];
                        }
                      else
                        {
                          double tan_angle_top = std::tan(interpolated_angle_top);

                          WBAssert(!std::isnan(tan_angle_top),
                                   "Internal error: The tan_angle_top variable is not a number: " << tan_angle_top);
                          const double center_circle_y = difference_in_angle_along_segment < 0 ?
                                                         begin_segment[1] - radius_angle_circle * cos_angle_top
                                                         : begin_segment[1] + radius_angle_circle * cos_angle_top;

                          WBAssert(!std::isnan(center_circle_y),
                                   "Internal error: The center_circle_y variable is not a number: " << center_circle_y
                                   << ". begin_segment[1] = " << begin_segment[1]
                                   << ", radius_angle_circle = " << radius_angle_circle
                                   << ", cos_angle_top = " << cos_angle_top);

                          // to prevent round off errors becomming dominant, we check
                          // whether center_circle_y - begin_segment[1] should be zero.
                          // TODO: improve this to some kind of relative difference.
                          const double CCYBS = center_circle_y - begin_segment[1];

                          WBAssert(!std::isnan(CCYBS),
                                   "Internal error: The CCYBS variable is not a number: " << CCYBS);



                          center_circle[0] = begin_segment[0] + tan_angle_top * (CCYBS);
                          center_circle[1] = center_circle_y;
                        }

                      WBAssert(!std::isnan(center_circle[0]) || !std::isnan(center_circle[1]),
                               "Internal error: The center variable contains not a number: " << center_circle[0] << ":" << center_circle[0]);
                      WBAssert(std::fabs((begin_segment-center_circle).norm() - std::fabs(radius_angle_circle))
                               < 1e-8 * std::fabs((begin_segment-center_circle).norm() + std::fabs(radius_angle_circle)),
                               "Internal error: The center of the circle is not a radius away from the begin point. " << std::endl
                               << "The center is located at " << center_circle[0] << ":" << center_circle[1] << std::endl
                               << "The begin point is located at " << begin_segment[0] << ":" << begin_segment[1] << std::endl
                               << "The computed radius is " << std::fabs((begin_segment-center_circle).norm())
                               << ", and it should be " << radius_angle_circle << ".");


                      // Now compute the location of the end of the segment by
                      // rotating P1 around the center_circle
                      Point<2> BSPC = begin_segment - center_circle;
                      const double sin_angle_diff = sin(difference_in_angle_along_segment);
                      const double cos_angle_diff = cos(difference_in_angle_along_segment);
                      end_segment[0] = cos_angle_diff * BSPC[0] - sin_angle_diff * BSPC[1] + center_circle[0];
                      end_segment[1] = sin_angle_diff * BSPC[0] + cos_angle_diff * BSPC[1] + center_circle[1];



                      WBAssert(std::fabs((end_segment-center_circle).norm() - std::fabs(radius_angle_circle))
                               < 1e-8 * std::fabs((end_segment-center_circle).norm() + std::fabs(radius_angle_circle)) ,
                               "Internal error: The center of the circle is not a radius away from the end point. " << std::endl
                               << "The center is located at " << center_circle[0] << ":" << center_circle[1] << std::endl
                               << "The end point is located at " << end_segment[0] << ":" << end_segment[1] << std::endl
                               << "The computed radius is " << std::fabs((end_segment-center_circle).norm())
                               << ", and it should be " << radius_angle_circle << ".");

                      // Now check if the angle of the check point in this circle
                      // is larger then the angle of P1 and smaller then P1 + angle
                      // difference. If that is the case then the distance from the
                      // plane is radius - (center - check_point).norm(). Otherwise
                      // it is infinity.
                      // The angle of the check point is computed with the help of
                      // dot product. But before that we need to adjust the check
                      // point 2d.
                      const Point<2> CPCR = check_point_2d - center_circle;
                      const double CPCR_norm = CPCR.norm();

                      const double dot_product = CPCR * Point<2>(0, radius_angle_circle, cartesian);
                      // If the x of the check point is larger then the x of center
                      // the circle, the angle is more than 180 degree, but the dot
                      // product will decrease instead of increase from 180 degrees.
                      // To fix this we make a special case for this.
                      // Furthermore, when the check point is at the same location as
                      // the center of the circle, we count that point as belonging
                      // to the top of the top segment (0 degree).
                      double check_point_angle = CPCR_norm == 0 ? 2.0 * const_pi : (check_point_2d[0] <= center_circle[0]
                                                                                    ? std::acos(dot_product/(CPCR_norm * radius_angle_circle))
                                                                                    : 2.0 * const_pi - std::acos(dot_product/(CPCR_norm * radius_angle_circle)));
                      check_point_angle = difference_in_angle_along_segment >= 0 ? const_pi - check_point_angle : 2.0 * const_pi - check_point_angle;

                      // In the case that it is exactly 2 * pi, bring it back to zero
                      check_point_angle = (std::fabs(check_point_angle - 2 * const_pi) < 1e-14 ? 0 : check_point_angle);

                      if ((difference_in_angle_along_segment > 0 && (check_point_angle <= interpolated_angle_top || std::fabs(check_point_angle - interpolated_angle_top) < 1e-12)
                           && (check_point_angle >= interpolated_angle_bottom || std::fabs(check_point_angle - interpolated_angle_bottom) < 1e-12))
                          || (difference_in_angle_along_segment < 0 && (check_point_angle >= interpolated_angle_top || std::fabs(check_point_angle - interpolated_angle_top) < 1e-12)
                              && (check_point_angle <= interpolated_angle_bottom || std::fabs(check_point_angle - interpolated_angle_bottom) < 1e-12)))
                        {
                          new_distance = (radius_angle_circle - CPCR_norm) * (difference_in_angle_along_segment < 0 ? 1 : -1);
                          new_along_plane_distance = (radius_angle_circle * check_point_angle - radius_angle_circle * interpolated_angle_top) * (difference_in_angle_along_segment < 0 ? 1 : -1);
                        }

                    }

                  // Now we need to see whether we need to update the information
                  // based on whether this segment is the closest one to the point
                  // up to now. To do this we first look whether the point falls
                  // within the bound of the segment and if it is actually closer.
                  // TODO: find out whether the fabs() are needed.
                  if (new_along_plane_distance >= -1e-10 &&
                      new_along_plane_distance <= std::fabs(interpolated_segment_length) &&
                      std::fabs(new_distance) < std::fabs(distance))
                    {
                      // There are two specific cases we are concerned with. The
                      // first case is that we want to have both the positive and
                      // negative distances (above and below the line). The second
                      // case is that we only want positive distances.
                      distance = only_positive ? std::fabs(new_distance) : new_distance;
                      along_plane_distance = new_along_plane_distance + total_length;
                      section = current_section;
                      section_fraction = fraction_CPL_P1P2;
                      segment = i_segment;
                      segment_fraction = new_along_plane_distance / interpolated_segment_length;
                      total_average_angle = (average_angle * total_length
                                             + 0.5 * (interpolated_angle_top + interpolated_angle_bottom  - 2 * add_angle) * new_along_plane_distance);
                      total_average_angle = (total_average_angle == 0 ? 0 : total_average_angle /
                                             (total_length + new_along_plane_distance));
                    }

                  // increase average angle
                  average_angle = (average_angle * total_length +
                                   0.5 * (interpolated_angle_top + interpolated_angle_bottom  - 2 * add_angle) * interpolated_segment_length);
                  average_angle = (average_angle == 0 ? 0 : average_angle /
                                   (total_length + interpolated_segment_length));
                  // increase the total length for the next segment.
                  total_length += interpolated_segment_length;
                }
            }
        }
      std::map<std::string, double> return_values;
      return_values["distanceFromPlane"] = distance;
      return_values["distanceAlongPlane"] = along_plane_distance;
      return_values["sectionFraction"] = section_fraction;
      return_values["segmentFraction"] = segment_fraction;
      return_values["section"] = section;
      return_values["segment"] = segment;
      return_values["averageAngle"] = total_average_angle;
      return return_values;
    }

    void interpolation::set_points(const std::vector<double> &x,
                                   const std::vector<double> &y,
                                   bool monotone_spline)
    {
      assert(x.size() == y.size());
      m_x = x;
      m_y = y;
      const unsigned int n = x.size();
      for (unsigned int i = 0; i < n-1; i++)
        {
          assert(m_x[i] < m_x[i+1]);
        }

      if (monotone_spline == true)
        {
          /**
           * This monotone spline algorithm is based on the javascript version
           * at https://en.wikipedia.org/wiki/Monotone_cubic_interpolation. The
           * parameters from this algorithm prevent overshooting in the
           * interpolation spline.
           */
          std::vector<double> dys(n-1), dxs(n-1), ms(n-1);
          for (unsigned int i=0; i < n-1; i++)
            {
              dxs[i] = x[i+1]-x[i];
              dys[i] = y[i+1]-y[i];
              ms[i] = dys[i]/dxs[i];
            }

          // get m_a parameter
          m_c.resize(n);
          m_c[0] = 0;

          for (unsigned int i = 0; i < n-2; i++)
            {
              const double m0 = ms[i];
              const double m1 = ms[i+1];

              if (m0 * m1 <= 0)
                {
                  m_c[i+1] = 0;
                }
              else
                {
                  const double dx0 = dxs[i];
                  const double dx1 = dxs[i+1];
                  const double common = dx0 + dx1;
                  m_c[i+1] = 3*common/((common + dx0)/m0 + (common + dx1)/m1);
                }
            }
          m_c[n-1] = ms[n-2];

          // Get b and c coefficients
          m_a.resize(n);
          m_b.resize(n);
          for (unsigned int i = 0; i < m_c.size()-1; i++)
            {
              const double c1 = m_c[i];
              const double m0 = ms[i];

              const double invDx = 1/dxs[i];
              const double common0 = c1 + m_c[i+1] - m0 - m0;
              m_b[i] = (m0 - c1 - common0) * invDx;
              m_a[i] = common0 * invDx * invDx;
            }
        }
      else     // linear interpolation
        {
          m_a.resize(n);
          m_b.resize(n);
          m_c.resize(n);
          for (unsigned int i = 0; i<n-1; i++)
            {
              m_a[i] = 0.0;
              m_b[i] = 0.0;
              m_c[i] = (m_y[i+1]-m_y[i])/(m_x[i+1]-m_x[i]);
            }
        }

      // for the right boundary we define
      // f_{n-1}(x) = b*(x-x_{n-1})^2 + c*(x-x_{n-1}) + y_{n-1}
      double h = x[n-1]-x[n-2];
      // m_b[n-1] is determined by the boundary condition
      if (!monotone_spline)
        {
          m_a[n-1] = 0.0;
          m_c[n-1] = 3.0*m_a[n-2]*h*h+2.0*m_b[n-2]*h+m_c[n-2];   // = f'_{n-2}(x_{n-1})
        }
    }

    double interpolation::operator() (double x) const
    {
      size_t n = m_x.size();
      // find the closest point m_x[idx] < x, idx=0 even if x<m_x[0]
      std::vector<double>::const_iterator it;
      it = std::lower_bound(m_x.begin(),m_x.end(),x);
      int idx = std::max( int(it-m_x.begin())-1, 0);

      double h = x-m_x[idx];
      double interpol;
      if (x<m_x[0])
        {
          // extrapolation to the left
          interpol = ((m_b[0])*h + m_c[0])*h + m_y[0];
        }
      else if (x>m_x[n-1])
        {
          // extrapolation to the right
          interpol = ((m_b[n-1])*h + m_c[n-1])*h + m_y[n-1];
        }
      else
        {
          // interpolation
          interpol = ((m_a[idx]*h + m_b[idx])*h + m_c[idx])*h + m_y[idx];
        }
      return interpol;
    }

    template const std::array<double,2> convert_point_to_array<2>(const Point<2> &point_);
    template const std::array<double,3> convert_point_to_array<3>(const Point<3> &point_);
  }
}



