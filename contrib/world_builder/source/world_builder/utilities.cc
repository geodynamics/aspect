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
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>

#ifdef WB_WITH_MPI
// we don't need the c++ MPI wrappers
#define OMPI_SKIP_MPICXX 1
#define MPICH_SKIP_MPICXX
#include <mpi.h>
#endif

#include "world_builder/nan.h"
#include "world_builder/utilities.h"
#include <world_builder/coordinate_system.h>


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
          other_point[0] += point[0] < 0 ? 2.0 * Consts::PI : -2.0 * Consts::PI;

          return (polygon_contains_point_implementation(point_list, point) ||
                  polygon_contains_point_implementation(point_list, other_point));
        }

      return polygon_contains_point_implementation(point_list, point);

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
      const size_t pointNo = point_list.size();
      size_t    wn = 0;    // the  winding number counter
      size_t   j=pointNo-1;

      // loop through all edges of the polygon
      for (size_t i=0; i<pointNo; i++)
        {
          // edge from V[i] to  V[i+1]
          if (point_list[j][1] <= point[1])
            {
              // first check if a point is directly on a line (within epsilon)
              if (approx(point_list[i][0],point[0]) && approx(point_list[i][1],point[1]))
                return true;
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
                  else if ( std::abs(is_left) < std::numeric_limits<double>::epsilon())
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
                  else if (std::abs(is_left) < std::numeric_limits<double>::epsilon())
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
    fraction_from_ellipse_center(const Point<2> &ellipse_center,
                                 const double semi_major_axis,
                                 const double eccentricity,
                                 const double theta,
                                 const Point<2> &point)
    {
      const double x_rotated = (point[0] - ellipse_center[0]) * std::cos(theta) + (point[1] - ellipse_center[1])* std::sin(theta);
      const double y_rotated = -(point[0] - ellipse_center[0]) * std::sin(theta) + (point[1] - ellipse_center[1])* std::cos(theta);

      const double semi_minor_axis = semi_major_axis * std::sqrt(1 - std::pow(eccentricity, 2));

      if (semi_major_axis < 10 * std::numeric_limits<double>::min()
          || semi_minor_axis < 10 * std::numeric_limits<double>::min())
        return false;

      const double ellipse = std::pow((x_rotated), 2) / std::pow(semi_major_axis, 2)
                             + std::pow((y_rotated), 2) / std::pow(semi_minor_axis, 2);

      return ellipse;
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

      const size_t n_poly_points = point_list.size();
      WBAssertThrow(n_poly_points >= 3, "Not enough polygon points were specified.");

      // Initialize a vector of distances for each point of the polygon with a very large distance
      std::vector<double> distances(n_poly_points, 1e23);

      // Create another polygon but with all points shifted 1 position to the right
      std::vector<Point<2> > shifted_point_list(n_poly_points, Point<2>(point.get_coordinate_system()));
      shifted_point_list[0] = point_list[n_poly_points-1];

      for (size_t i = 0; i < n_poly_points-1; ++i)
        shifted_point_list[i+1] = point_list[i];

      for (size_t i = 0; i < n_poly_points; ++i)
        {
          // Create vector along the polygon line segment
          const Point<2> vector_segment = shifted_point_list[i] - point_list[i];
          // Create vector from point to the second segment point
          const Point<2> vector_point_segment = point - point_list[i];

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


    std::array<double,3>
    cartesian_to_spherical_coordinates(const Point<3> &position)
    {
      std::array<double,3> scoord;

      scoord[0] = position.norm(); // R
      scoord[1] = std::atan2(position[1],position[0]); // Phi/long -> The result is always between -180 and 180 degrees: [-pi,pi]
      //if (scoord[1] < 0.0)
      //scoord[1] += 2.0*Consts::PI; // correct phi to [0,2*pi]

      //lat
      if (scoord[0] > std::numeric_limits<double>::min())
        scoord[2] = 0.5 * Consts::PI - std::acos(position[2]/scoord[0]);
      else
        scoord[2] = 0.0;

      return scoord;
    }

    Point<3>
    spherical_to_cartesian_coordinates(const std::array<double,3> &scoord)
    {
      const double cos_lat = scoord[0] * std::sin(0.5 * Consts::PI - scoord[2]);

      return Point<3>(cos_lat * std::cos(scoord[1]), // X
                      cos_lat * std::sin(scoord[1]), // Y
                      scoord[0] * std::cos(0.5 * Consts::PI - scoord[2]), // Z
                      cartesian);;
    }



    CoordinateSystem
    string_to_coordinate_system(const std::string &coordinate_system)
    {
      if (coordinate_system == "cartesian")
        return CoordinateSystem::cartesian;
      if (coordinate_system == "spherical")
        return CoordinateSystem::spherical;
      WBAssertThrow(false, "Coordinate system not implemented.");

      return invalid;
    }


    template<unsigned int dim>
    std::array<double,dim>
    convert_point_to_array(const Point<dim> &point_)
    {
      std::array<double,dim> array;
      for (size_t i = 0; i < dim; ++i)
        array[i] = point_[i];
      return array;
    }

    double
    string_to_double(const std::string &string)
    {
      // trim whitespace on either side of the text if necessary
      std::string s = string;
      while ((!s.empty()) && (s[0] == ' '))
        s.erase(s.begin());
      while ((!s.empty()) && (s[s.size() - 1] == ' '))
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
      while ((!s.empty()) && (s[0] == ' '))
        s.erase(s.begin());
      while ((!s.empty()) && (s[s.size() - 1] == ' '))
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
      while ((!s.empty()) && (s[0] == ' '))
        s.erase(s.begin());
      while ((!s.empty()) && (s[s.size() - 1] == ' '))
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

    PointDistanceFromCurvedPlanes
    distance_point_from_curved_planes(const Point<3> &check_point, // cartesian point in cartesian and spherical system
                                      const Objects::NaturalCoordinate &natural_coordinate, // cartesian point cartesian system, spherical point in spherical system
                                      const Point<2> &reference_point, // in (rad) spherical coordinates in spherical system
                                      const std::vector<Point<2> > &point_list, // in  (rad) spherical coordinates in spherical system
                                      const std::vector<std::vector<double> > &plane_segment_lengths,
                                      const std::vector<std::vector<Point<2> > > &plane_segment_angles,
                                      const double start_radius,
                                      const std::unique_ptr<CoordinateSystems::Interface> &coordinate_system,
                                      const bool only_positive,
                                      const Objects::BezierCurve &bezier_curve)
    {
      double distance = std::numeric_limits<double>::infinity();
      double new_distance = std::numeric_limits<double>::infinity();
      double along_plane_distance = std::numeric_limits<double>::infinity();
      double new_along_plane_distance  = std::numeric_limits<double>::infinity();
      double new_depth_reference_surface = std::numeric_limits<double>::infinity();

      const CoordinateSystem natural_coordinate_system = coordinate_system->natural_coordinate_system();
      const bool bool_cartesian = natural_coordinate_system == cartesian;

      const std::array<double,3> &check_point_surface_2d_array = natural_coordinate.get_coordinates();
      const Point<3> check_point_surface(bool_cartesian ? check_point_surface_2d_array[0] : start_radius,
                                         check_point_surface_2d_array[1],
                                         bool_cartesian ? start_radius : check_point_surface_2d_array[2],
                                         natural_coordinate_system);
      const Point<2> check_point_surface_2d(natural_coordinate.get_surface_coordinates(),
                                            natural_coordinate_system);

      // The section which is checked.
      size_t section = 0;

      // The 'horizontal' fraction between the points at the surface.
      double section_fraction = 0.0;

      // What segment the point on the line is in.
      size_t segment = 0;

      // The 'vertical' fraction, indicates how far in the current segment the
      // point on the line is.
      double segment_fraction = 0.0;
      double total_average_angle = 0.0;
      double depth_reference_surface = 0.0;

      const DepthMethod depth_method = coordinate_system->depth_method();

      size_t i_section_min_distance = 0;
      double fraction_CPL_P1P2 = std::numeric_limits<double>::infinity();

      // get an estimate for the closest point between P1 and P2.
      //constexpr double parts = 1;
      //constexpr double one_div_parts = 1./parts;
      //double minimum_distance_to_reference_point = std::numeric_limits<double>::infinity();
      //const size_t number_of_points = point_list.size();

      const Objects::ClosestPointOnCurve closest_point_on_curve = bezier_curve.closest_point_on_curve_segment(check_point_surface_2d);
      Point<2> closest_point_on_line_2d = closest_point_on_curve.point;

      // We now need 3d points from this point on, so make them.
      // The order of a Cartesian coordinate is x,y,z and the order of
      // a spherical coordinate it radius, long, lat (in rad).
      const Point<3> closest_point_on_line_surface(bool_cartesian ? closest_point_on_line_2d[0] : start_radius,
                                                   bool_cartesian ? closest_point_on_line_2d[1] : closest_point_on_line_2d[0],
                                                   bool_cartesian ? start_radius : closest_point_on_line_2d[1],
                                                   natural_coordinate_system);

      const Point<3> closest_point_on_line_cartesian(coordinate_system->natural_to_cartesian_coordinates(closest_point_on_line_surface.get_array()),cartesian);

      if (!std::isnan(closest_point_on_line_2d[0]))
        {
          i_section_min_distance = closest_point_on_curve.index;
          fraction_CPL_P1P2 = closest_point_on_curve.parametric_fraction;

          Point<3> closest_point_on_line_bottom = closest_point_on_line_surface;
          closest_point_on_line_bottom[bool_cartesian ? 2 : 0] = 0;

          WBAssert(!std::isnan(closest_point_on_line_bottom[0])
                   &&
                   !std::isnan(closest_point_on_line_bottom[1])
                   &&
                   !std::isnan(closest_point_on_line_bottom[2]),
                   "Internal error: The closest_point_on_line_bottom variables variable contains not a number: " << closest_point_on_line_bottom);

          // Now that we have both the check point and the
          // closest_point_on_line, we need to push them to cartesian.
          const Point<3> closest_point_on_line_bottom_cartesian(coordinate_system->natural_to_cartesian_coordinates(closest_point_on_line_bottom.get_array()),cartesian);
          const Point<3> check_point_surface_cartesian(coordinate_system->natural_to_cartesian_coordinates(check_point_surface.get_array()),cartesian);


          WBAssert(!std::isnan(closest_point_on_line_bottom_cartesian[0]),
                   "Internal error: The closest_point_on_line_bottom_cartesian[0] variable is not a number: " << closest_point_on_line_bottom_cartesian[0]);
          WBAssert(!std::isnan(closest_point_on_line_bottom_cartesian[1]),
                   "Internal error: The closest_point_on_line_bottom_cartesian[1] variable is not a number: " << closest_point_on_line_bottom_cartesian[1]);
          WBAssert(!std::isnan(closest_point_on_line_bottom_cartesian[2]),
                   "Internal error: The closest_point_on_line_bottom_cartesian[2] variable is not a number: " << closest_point_on_line_bottom_cartesian[2]);


          // translate to original coordinates current and next section
          const size_t original_current_section = i_section_min_distance;
          const size_t original_next_section = original_current_section + 1;


          // These are the mostly likely cases for the x and y axis, so initialize them to these values. They will be checked
          // in the else statement or replaced in the if statement.
          Point<3> y_axis = closest_point_on_line_cartesian - closest_point_on_line_bottom_cartesian;
          Point<3> x_axis = closest_point_on_line_cartesian - check_point_surface_cartesian;

          // This are accounting for corner cases.
          // If the point to check is exactly on or below the line, we can not compute the x-axis with this method.
          // We could use an other method where we use the two point before and after it, but we can also
          // just nudge it into a direction, which seems to work very well.
          if (std::fabs((check_point_surface - closest_point_on_line_surface).norm()) < 2e-14)
            {
              // If the point to check is on the line, we don't need to search any further, because we know the distance is zero.
              if (std::fabs((check_point - closest_point_on_line_cartesian).norm()) > 2e-14)
                {
                  const Point<2> &P1(point_list[i_section_min_distance]);
                  const Point<2> &P2(point_list[i_section_min_distance+1]);

                  const Point<2> P1P2 = P2 - P1;
                  const Point<2> unit_normal_to_plane_spherical = P1P2 / P1P2.norm();
                  const Point<2> closest_point_on_line_plus_normal_to_plane_spherical = closest_point_on_line_2d + 1e-8 * (closest_point_on_line_2d.norm() > 1.0 ? closest_point_on_line_2d.norm() : 1.0) * unit_normal_to_plane_spherical;

                  WBAssert(std::fabs(closest_point_on_line_plus_normal_to_plane_spherical.norm()) > std::numeric_limits<double>::epsilon(),
                           "Internal error: The norm of variable 'closest_point_on_line_plus_normal_to_plane_spherical' "
                           "is  zero, while this may not happen.");

                  const Point<3> closest_point_on_line_plus_normal_to_plane_surface_spherical(bool_cartesian ? closest_point_on_line_plus_normal_to_plane_spherical[0] : start_radius,
                                                                                              bool_cartesian ? closest_point_on_line_plus_normal_to_plane_spherical[1] : closest_point_on_line_plus_normal_to_plane_spherical[0],
                                                                                              bool_cartesian ? start_radius : closest_point_on_line_plus_normal_to_plane_spherical[1],
                                                                                              natural_coordinate_system);
                  const Point<3> closest_point_on_line_plus_normal_to_plane_cartesian(coordinate_system->natural_to_cartesian_coordinates(closest_point_on_line_plus_normal_to_plane_surface_spherical.get_array()),cartesian);
                  Point<3> normal_to_plane = closest_point_on_line_plus_normal_to_plane_cartesian - closest_point_on_line_cartesian;
                  normal_to_plane = normal_to_plane / normal_to_plane.norm();

                  // The y-axis is from the bottom/center to the closest_point_on_line,
                  // the x-axis is 90 degrees rotated from that, so we rotate around
                  // the line P1P2.
                  // Todo: Assert that the norm of the axis are not equal to zero.
                  y_axis = closest_point_on_line_cartesian - closest_point_on_line_bottom_cartesian;

                  WBAssert(std::abs(y_axis.norm()) > std::numeric_limits<double>::epsilon(),
                           "World Builder error: Cannot determine the up direction in the model. This is most likely due to the provided start radius being zero."
                           << " Technical details: The y_axis.norm() is zero. Y_axis is " << y_axis[0] << ':' << y_axis[1] << ':' << y_axis[2]
                           << ". closest_point_on_line_cartesian = " << closest_point_on_line_cartesian[0] << ':' << closest_point_on_line_cartesian[1] << ':' << closest_point_on_line_cartesian[2]
                           << ", closest_point_on_line_bottom_cartesian = " << closest_point_on_line_bottom_cartesian[0] << ':' << closest_point_on_line_bottom_cartesian[1] << ':' << closest_point_on_line_bottom_cartesian[2]);

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
                  const double vx = y_axis[0];
                  const double vy = y_axis[1];
                  const double vz = y_axis[2];
                  const double ux = normal_to_plane[0];
                  const double uy = normal_to_plane[1];
                  const double uz = normal_to_plane[2];

                  x_axis = Point<3>(ux*ux*vx + ux*uy*vy - uz*vy + uy*uz*vz + uy*vz,
                                    uy*ux*vx + uz*vx + uy*uy*vy + uy*uz*vz - ux*vz,
                                    uz*ux*vx - uy*vx + uz*uy*vy + ux*vy + uz*uz*vz,
                                    cartesian);

                  // see on what side the line P1P2 reference point is. This is based on the determinant
                  const Point<2> reference_p = ((closest_point_on_curve.normal-closest_point_on_line_2d)*1e2)+closest_point_on_line_2d;
                  const double reference_on_side_of_line =  (closest_point_on_line_2d-reference_p).norm_square() < (check_point_surface_2d-reference_p).norm_square() ? -1 : 1;

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
                }
              else
                {
                  total_average_angle = plane_segment_angles[original_current_section][0][0]
                                        + fraction_CPL_P1P2 * (plane_segment_angles[original_next_section][0][0]
                                                               - plane_segment_angles[original_current_section][0][0]);

                  PointDistanceFromCurvedPlanes return_values(natural_coordinate.get_coordinate_system());
                  return_values.distance_from_plane = 0.0;
                  return_values.distance_along_plane = 0.0;
                  return_values.fraction_of_section = fraction_CPL_P1P2;
                  return_values.fraction_of_segment = 0.0;
                  return_values.section = i_section_min_distance;
                  return_values.segment = 0;
                  return_values.average_angle = total_average_angle;
                  return_values.depth_reference_surface = 0.0;
                  return_values.closest_trench_point = closest_point_on_line_cartesian;
                  return return_values;
                }
            }
          else
            {
              WBAssert(std::abs(y_axis.norm()) > std::numeric_limits<double>::epsilon(),
                       "World Builder error: Cannot determine the up direction in the model. This is most likely due to the provided start radius being zero."
                       << " Technical details: The y_axis.norm() is zero. Y_axis is " << y_axis
                       << ". closest_point_on_line_cartesian = " << closest_point_on_line_cartesian
                       << ", closest_point_on_line_bottom_cartesian = " << closest_point_on_line_bottom_cartesian);

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


              Point<2> check_point_surface_2d_temp = check_point_surface_2d;

              if (!bool_cartesian)
                {
                  const double normal = std::fabs(point_list[i_section_min_distance+static_cast<size_t>(std::round(fraction_CPL_P1P2))][0]-check_point_surface_2d[0]);
                  const double plus   = std::fabs(point_list[i_section_min_distance+static_cast<size_t>(std::round(fraction_CPL_P1P2))][0]-(check_point_surface_2d[0]+2*Consts::PI));
                  const double min    = std::fabs(point_list[i_section_min_distance+static_cast<size_t>(std::round(fraction_CPL_P1P2))][0]-(check_point_surface_2d[0]-2*Consts::PI));

                  // find out whether the check point, checkpoint + 2pi or check point -2 pi is closest to the point list.
                  if (plus < normal)
                    {
                      check_point_surface_2d_temp[0]+= 2*Consts::PI;
                    }
                  else if (min < normal)
                    {
                      check_point_surface_2d_temp[0]-= 2*Consts::PI;
                    }
                }

              // check whether the check point and the reference point are on the same side, if not, change the side.
              const Point<2> AB_normal = closest_point_on_curve.normal*closest_point_on_line_2d.distance(reference_point);//*AB.norm();
              const Point<2> local_reference_point = AB_normal*1.+closest_point_on_line_2d;
              const bool reference_normal_on_side_of_line =  (closest_point_on_line_2d-local_reference_point).norm_square() < (check_point_surface_2d_temp-local_reference_point).norm_square();
              const bool reference_point_on_side_of_line =  (point_list[point_list.size()-1][0] - point_list[0][0])*(reference_point[1] - point_list[0][1]) - (reference_point[0] - point_list[0][0])*(point_list[point_list.size()-1][1] - point_list[0][1]) < 0.;
              const double reference_on_side_of_line =  reference_normal_on_side_of_line == reference_point_on_side_of_line ? 1 : -1;

              WBAssert(!std::isnan(x_axis[0]),
                       "Internal error: The x_axis variable is not a number: " << x_axis[0]);
              WBAssert(!std::isnan(x_axis[1]),
                       "Internal error: The x_axis variable is not a number: " << x_axis[1]);
              WBAssert(!std::isnan(x_axis[2]),
                       "Internal error: The x_axis variable is not a number: " << x_axis[2]);

              WBAssert(x_axis.norm() > 0.0, "x_axis norm is zero");

              x_axis = x_axis *(reference_on_side_of_line / x_axis.norm());

              WBAssert(!std::isnan(x_axis[0]),
                       "Internal error: The x_axis variable is not a number: " << x_axis[0]);
              WBAssert(!std::isnan(x_axis[1]),
                       "Internal error: The x_axis variable is not a number: " << x_axis[1]);
              WBAssert(!std::isnan(x_axis[2]),
                       "Internal error: The x_axis variable is not a number: " << x_axis[2]);

            }

          WBAssert(!std::isnan(x_axis[0]),
                   "Internal error: The x_axis[0] variable is not a number: " << x_axis[0] << ". Relevant values:  check_point = " << check_point << '.');
          WBAssert(!std::isnan(x_axis[1]),
                   "Internal error: The x_axis[1] variable is not a number: " << x_axis[1]);
          WBAssert(!std::isnan(x_axis[2]),
                   "Internal error: The x_axis[2] variable is not a number: " << x_axis[2]);

          // now that we have the x and y axes computed, convert the 3d check point into a 2d one.
          Point<2> check_point_2d(x_axis * (check_point - closest_point_on_line_bottom_cartesian),
                                  y_axis * (check_point - closest_point_on_line_bottom_cartesian),
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
          double add_angle_correction = 0.0;
          double average_angle = 0.0;
          for (size_t i_segment = 0; i_segment < plane_segment_lengths[original_current_section].size(); i_segment++)
            {
              const size_t current_segment = i_segment;

              // compute the angle between the previous begin and end if
              // the depth method is angle_at_begin_segment_with_surface.
              if (i_segment != 0
                  &&
                  (depth_method == DepthMethod::angle_at_begin_segment_with_surface
                   ||
                   depth_method == DepthMethod::angle_at_begin_segment_applied_to_end_segment_with_surface))
                {
                  double add_angle_inner = (begin_segment * end_segment) / (begin_segment.norm() * end_segment.norm());

                  WBAssert(!std::isnan(add_angle_inner),
                           "Internal error: The add_angle_inner variable is not a number: " << add_angle_inner
                           << ". Variables: begin_segment = " << begin_segment
                           << ", end_segment = " << end_segment
                           << ", begin_segment * end_segment / (begin_segment.norm() * end_segment.norm()) = "
                           << std::setprecision(32) << begin_segment * end_segment / (begin_segment.norm() * end_segment.norm())
                           << '.');

                  // there could be round of error problems here is the inner part is close to one
                  if (add_angle_inner < 0. && add_angle_inner >= -1e-14)
                    add_angle_inner = 0.;
                  if (add_angle_inner > 1. && add_angle_inner <= 1.+1e-14)
                    add_angle_inner = 1.;

                  WBAssert(add_angle_inner >= 0 && add_angle_inner <= 1,
                           "Internal error: The variable add_angle_inner is smaller than zero or larger then one,"
                           "which causes the std::acos to return nan. If it is only a little bit larger then one, "
                           "this is probably caused by that begin and end segment are the same and round off error. "
                           "The value of add_angle_inner = " << add_angle_inner);

                  add_angle_correction = std::acos(add_angle_inner);
                  add_angle += add_angle_correction;

                  WBAssert(!std::isnan(add_angle),
                           "Internal error: The add_angle variable is not a number: " << add_angle
                           << ". Variables: begin_segment = " << begin_segment
                           << ", end_segment = " << end_segment
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
              const double degree_90_to_rad = 0.5 * Consts::PI;

              WBAssert(plane_segment_angles.size() > original_next_section,
                       "Error: original_next_section = " << original_next_section
                       << ", and plane_segment_angles.size() = " << plane_segment_angles.size());


              WBAssert(plane_segment_angles[original_next_section].size() > current_segment,
                       "Error: current_segment = "  << current_segment
                       << ", and current_segment.size() = " << plane_segment_angles[original_next_section].size());

              const double interpolated_angle_top    = plane_segment_angles[original_current_section][current_segment][0]
                                                       + fraction_CPL_P1P2 * (plane_segment_angles[original_next_section][current_segment][0]
                                                                              - plane_segment_angles[original_current_section][current_segment][0])
                                                       + add_angle
                                                       + (depth_method == DepthMethod::angle_at_begin_segment_applied_to_end_segment_with_surface
                                                          && i_segment != 0 ? -add_angle_correction: 0);

              const double interpolated_angle_bottom = plane_segment_angles[original_current_section][current_segment][1]
                                                       + fraction_CPL_P1P2 * (plane_segment_angles[original_next_section][current_segment][1]
                                                                              - plane_segment_angles[original_current_section][current_segment][1])
                                                       + add_angle;


              const double interpolated_segment_length     = plane_segment_lengths[original_current_section][current_segment]
                                                             + fraction_CPL_P1P2 * (plane_segment_lengths[original_next_section][current_segment]
                                                                                    - plane_segment_lengths[original_current_section][current_segment]);

              if (interpolated_segment_length < 1e-14)
                continue;

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
                  if (std::fabs(interpolated_segment_length) > std::numeric_limits<double>::epsilon())
                    {
                      end_segment[0] += interpolated_segment_length * std::sin(degree_90_to_rad - interpolated_angle_top);
                      end_segment[1] -= interpolated_segment_length * std::cos(degree_90_to_rad - interpolated_angle_top);

                      Point<2> begin_end_segment = end_segment - begin_segment;
                      Point<2> normal_2d_plane(-begin_end_segment[0],begin_end_segment[1], cartesian);
                      WBAssert(std::fabs(normal_2d_plane.norm()) > std::numeric_limits<double>::epsilon(), "Internal Error: normal_2d_plane.norm() is zero, which should not happen. "
                               << "Extra info: begin_end_segment[0] = " << begin_end_segment[0]
                               << ", begin_end_segment[1] = " << begin_end_segment[1]
                               << ", end_segment: [" << end_segment[0] << ',' << end_segment[1] << ']'
                               << ", begin_segment: [" << begin_segment[0] << ',' << begin_segment[1] << ']'
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
                          new_distance = std::numeric_limits<double>::infinity();
                          new_along_plane_distance = std::numeric_limits<double>::infinity();
                          new_depth_reference_surface = std::numeric_limits<double>::infinity();
                        }
                      else
                        {
                          const Point<2> Pb = begin_segment + (c1/c2) * BSP_ESP;
                          const double side_of_line =  (begin_segment[0] - end_segment[0]) * (check_point_2d[1] - begin_segment[1])
                                                       - (begin_segment[1] - end_segment[1]) * (check_point_2d[0] - begin_segment[0])
                                                       < 0 ? -1.0 : 1.0;

                          new_distance = side_of_line * (check_point_2d - Pb).norm();
                          new_along_plane_distance = (begin_segment - Pb).norm();
                          new_depth_reference_surface = start_radius - Pb[1];

                          WBAssert(!std::isnan(new_depth_reference_surface),
                                   "new_depth_reference_surface is not a number: " << new_depth_reference_surface << ". "
                                   << "start_radius = " << start_radius << ",Pb[1] = " << Pb[1] << ".");
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
                  if (std::fabs(interpolated_angle_top - 0.5 * Consts::PI) < 1e-8)
                    {
                      // if interpolated_angle_top is 90 degrees, the tan function
                      // is undefined (1/0). What we really want in this case is
                      // set the center to the correct location which is x = the x
                      //begin point + radius and y = the y begin point.
                      center_circle[0] = difference_in_angle_along_segment > 0 ? begin_segment[0] + radius_angle_circle : begin_segment[0] - radius_angle_circle;
                      center_circle[1] = begin_segment[1];
                    }
                  else if (std::fabs(interpolated_angle_top - 1.5 * Consts::PI) < 1e-8)
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
                      const double tan_angle_top = std::tan(interpolated_angle_top);

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

                      // to prevent round off errors becoming dominant, we check
                      // whether center_circle_y - begin_segment[1] should be zero.
                      // TODO: improve this to some kind of relative difference.
                      const double CCYBS = center_circle_y - begin_segment[1];

                      WBAssert(!std::isnan(CCYBS),
                               "Internal error: The CCYBS variable is not a number: " << CCYBS);



                      center_circle[0] = begin_segment[0] + tan_angle_top * (CCYBS);
                      center_circle[1] = center_circle_y;
                    }

                  WBAssert(!std::isnan(center_circle[0]) || !std::isnan(center_circle[1]),
                           "Internal error: The center variable contains not a number: " << center_circle[0] << ':' << center_circle[0]);
                  WBAssert(std::fabs((begin_segment-center_circle).norm() - std::fabs(radius_angle_circle))
                           < 1e-8 * std::fabs((begin_segment-center_circle).norm() + std::fabs(radius_angle_circle)),
                           "Internal error: The center of the circle is not a radius away from the begin point. " << std::endl
                           << "The center is located at " << center_circle[0] << ':' << center_circle[1] << std::endl
                           << "The begin point is located at " << begin_segment[0] << ':' << begin_segment[1] << std::endl
                           << "The computed radius is " << std::fabs((begin_segment-center_circle).norm())
                           << ", and it should be " << radius_angle_circle << '.');


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
                           << "The center is located at " << center_circle[0] << ':' << center_circle[1] << std::endl
                           << "The end point is located at " << end_segment[0] << ':' << end_segment[1] << std::endl
                           << "The computed radius is " << std::fabs((end_segment-center_circle).norm())
                           << ", and it should be " << radius_angle_circle << '.');

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
                  double check_point_angle = std::fabs(CPCR_norm) < std::numeric_limits<double>::epsilon() ? 2.0 * Consts::PI : (check_point_2d[0] <= center_circle[0]
                                             ? std::acos(dot_product/(CPCR_norm * radius_angle_circle))
                                             : 2.0 * Consts::PI - std::acos(dot_product/(CPCR_norm * radius_angle_circle)));
                  check_point_angle = difference_in_angle_along_segment >= 0 ? Consts::PI - check_point_angle : 2.0 * Consts::PI - check_point_angle;

                  // In the case that it is exactly 2 * pi, bring it back to zero
                  check_point_angle = (std::fabs(check_point_angle - 2 * Consts::PI) < 1e-14 ? 0 : check_point_angle);

                  if ((difference_in_angle_along_segment > 0 && (check_point_angle <= interpolated_angle_top || std::fabs(check_point_angle - interpolated_angle_top) < 1e-12)
                       && (check_point_angle >= interpolated_angle_bottom || std::fabs(check_point_angle - interpolated_angle_bottom) < 1e-12))
                      || (difference_in_angle_along_segment < 0 && (check_point_angle >= interpolated_angle_top || std::fabs(check_point_angle - interpolated_angle_top) < 1e-12)
                          && (check_point_angle <= interpolated_angle_bottom || std::fabs(check_point_angle - interpolated_angle_bottom) < 1e-12)))
                    {
                      new_distance = (radius_angle_circle - CPCR_norm) * (difference_in_angle_along_segment < 0 ? 1 : -1);
                      new_along_plane_distance = (radius_angle_circle * check_point_angle - radius_angle_circle * interpolated_angle_top) * (difference_in_angle_along_segment < 0 ? 1 : -1);
                      // compute the new depth by rotating the begin point to the check point location.
                      new_depth_reference_surface = start_radius-(sin(check_point_angle + interpolated_angle_top) * BSPC[0] + cos(check_point_angle + interpolated_angle_top) * BSPC[1] + center_circle[1]);

                      WBAssert(!std::isnan(new_depth_reference_surface),
                               "new_depth_reference_surface is not a number: " << new_depth_reference_surface << ". "
                               << "start_radius = " << start_radius << ", check_point_angle = " << check_point_angle << ", interpolated_angle_top = " << interpolated_angle_top
                               << ", BSPC[0] = " << BSPC[0] << ".");
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
                  section = i_section_min_distance;
                  section_fraction = fraction_CPL_P1P2;
                  segment = i_segment;
                  segment_fraction = new_along_plane_distance / interpolated_segment_length;
                  total_average_angle = (average_angle * total_length
                                         + 0.5 * (interpolated_angle_top + interpolated_angle_bottom  - 2 * add_angle) * new_along_plane_distance);
                  total_average_angle = (std::fabs(total_average_angle) < std::numeric_limits<double>::epsilon() ? 0 : total_average_angle /
                                         (total_length + new_along_plane_distance));
                  depth_reference_surface = new_depth_reference_surface;
                }

              // increase average angle
              average_angle = (average_angle * total_length +
                               0.5 * (interpolated_angle_top + interpolated_angle_bottom  - 2 * add_angle) * interpolated_segment_length);
              average_angle = (std::fabs(average_angle) < std::numeric_limits<double>::epsilon() ? 0 : average_angle /
                               (total_length + interpolated_segment_length));
              // increase the total length for the next segment.
              total_length += interpolated_segment_length;
            }
        }

      WBAssert(!std::isnan(depth_reference_surface), "depth_reference_surface is not a number: " << depth_reference_surface << ".");

      PointDistanceFromCurvedPlanes return_values(natural_coordinate.get_coordinate_system());
      return_values.distance_from_plane = distance;
      return_values.distance_along_plane = along_plane_distance;
      return_values.fraction_of_section = section_fraction;
      return_values.fraction_of_segment = segment_fraction;
      return_values.section = section;
      return_values.segment = segment;
      return_values.average_angle = total_average_angle;
      return_values.depth_reference_surface = depth_reference_surface;
      return_values.closest_trench_point = closest_point_on_line_cartesian;
      return return_values;
    }

    void interpolation::set_points(const std::vector<double> &y)
    {
      const size_t n = y.size();
      mx_size_min = n;
      m.resize(n);
      for (unsigned int i = 0; i < n; ++i)
        {
          m[i][3] = y[i];
        }

      /**
       * This monotone spline algorithm is based on the javascript version
       * at https://en.wikipedia.org/wiki/Monotone_cubic_interpolation. The
       * parameters from this algorithm prevent overshooting in the
       * interpolation spline.
       */

      // get m_a parameter
      m[0][2] = 0;

      for (size_t i = 0; i < n-2; i++)
        {
          const double m0 = y[i+1]-y[i];
          const double m1 =  y[i+2]-y[i+1];

          if (m0 * m1 <= 0)
            {
              m[i+1][2] = 0;
            }
          else
            {
              m[i+1][2] = 2*m0*m1/(m0+m1);
            }
        }
      m[n-1][2] =  y[n-1]-y[n-2];

      // Get b and c coefficients
      //m_a.resize(n);
      //m_b.resize(n);
      for (size_t i = 0; i < n-1; i++)
        {
          const double c1 = m[i][2];
          const double m0 = y[i+1]-y[i];

          const double common0 = c1 + m[i+1][2] - m0 - m0;
          m[i][1] = (m0 - c1 - common0);
          m[i][0] = common0;
        }
    }

    double wrap_angle(const double angle)
    {
      return angle - 360.0*std::floor(angle/360.0);
    }

    double interpolate_angle_across_zero(const double angle_1,
                                         const double angle_2,
                                         const double fraction)
    {
      double theta_1 = angle_1;
      double theta_2 = angle_2;
      double rotation_angle;

      if (std::abs(theta_2 - theta_1) > Consts::PI)
        {
          if (theta_2 > theta_1)
            theta_1 += 2.*Consts::PI;
          else
            theta_2 += 2.*Consts::PI;
        }
      rotation_angle = (1-fraction) * theta_1 + fraction * theta_2;

      // make sure angle is between 0 and 360 degrees
      rotation_angle = rotation_angle - 2*Consts::PI*std::floor(rotation_angle/(2 * Consts::PI));

      return rotation_angle;
    }

    std::array<double,3>
    euler_angles_from_rotation_matrix(const std::array<std::array<double,3>,3> &rotation_matrix)
    {
      const double rad_to_degree = 180.0/Consts::PI;
      std::array<double,3> euler_angles;
      //const double s2 = std::sqrt(rotation_matrix[2][1] * rotation_matrix[2][1] + rotation_matrix[2][0] * rotation_matrix[2][0]);
      const std::ostringstream os;
      for (size_t i = 0; i < 3; i++)
        for (size_t j = 0; j < 3; j++)
          WBAssert(std::fabs(rotation_matrix[i][j]) <= 1.0,
                   "rotation_matrix[" + std::to_string(i) + "][" + std::to_string(j) +
                   "] is larger than one: " + std::to_string(rotation_matrix[i][j]) + ". rotation_matrix = \n"
                   + std::to_string(rotation_matrix[0][0]) + " " + std::to_string(rotation_matrix[0][1]) + " " + std::to_string(rotation_matrix[0][2]) + "\n"
                   + std::to_string(rotation_matrix[1][0]) + " " + std::to_string(rotation_matrix[1][1]) + " " + std::to_string(rotation_matrix[1][2]) + "\n"
                   + std::to_string(rotation_matrix[2][0]) + " " + std::to_string(rotation_matrix[2][1]) + " " + std::to_string(rotation_matrix[2][2]));


      const double theta = std::acos(rotation_matrix[2][2]);
      const double phi1  = std::atan2(rotation_matrix[2][0]/-sin(theta),rotation_matrix[2][1]/-sin(theta));
      const double phi2  = std::atan2(rotation_matrix[0][2]/-sin(theta),rotation_matrix[1][2]/sin(theta));

      euler_angles[0] = wrap_angle(phi1 * rad_to_degree);
      euler_angles[1] = wrap_angle(theta * rad_to_degree);
      euler_angles[2] = wrap_angle(phi2 * rad_to_degree);

      return euler_angles;
    }

    std::array<std::array<double,3>,3>
    euler_angles_to_rotation_matrix(double phi1_d, double theta_d, double phi2_d)
    {

      const double degree_to_rad = Consts::PI/180.0;
      const double phi1 = phi1_d * degree_to_rad;
      const double theta = theta_d * degree_to_rad;
      const double phi2 = phi2_d * degree_to_rad;
      std::array<std::array<double,3>,3> rot_matrix;


      rot_matrix[0][0] = cos(phi2)*cos(phi1) - cos(theta)*sin(phi1)*sin(phi2);
      rot_matrix[0][1] = -cos(phi2)*sin(phi1) - cos(theta)*cos(phi1)*sin(phi2);
      rot_matrix[0][2] = -sin(phi2)*sin(theta);

      rot_matrix[1][0] = sin(phi2)*cos(phi1) + cos(theta)*sin(phi1)*cos(phi2);
      rot_matrix[1][1] = -sin(phi2)*sin(phi1) + cos(theta)*cos(phi1)*cos(phi2);
      rot_matrix[1][2] = cos(phi2)*sin(theta);

      rot_matrix[2][0] = -sin(theta)*sin(phi1);
      rot_matrix[2][1] = -sin(theta)*cos(phi1);
      rot_matrix[2][2] = cos(theta);
      return rot_matrix;
    }

    std::string
    read_and_distribute_file_content(const std::string &filename)
    {
      std::string data_string;

#ifdef WB_WITH_MPI
      int mpi_initialized;
      MPI_Initialized(&mpi_initialized);
      if (mpi_initialized != 0)
        {
          const MPI_Comm comm = MPI_COMM_WORLD;
          int my_rank = 0;
          MPI_Comm_rank(comm, &my_rank);
          if (my_rank == 0)
            {
              int filesize;
              std::ifstream filestream;
              filestream.open(filename.c_str());
              WBAssertThrow (filestream.is_open(), std::string("Could not open file <") + filename + ">.");

              // Need to convert to unsigned int, because MPI_Bcast does not support
              // size_t or const unsigned int
              int invalid_file_size = -1;

              if (!filestream)
                {
                  // broadcast failure state, then throw
                  const int ierr = MPI_Bcast(&invalid_file_size, 1, MPI_UNSIGNED, 0, comm);
                  WBAssertThrow (ierr,
                                 std::string("Could not open file <") + filename + ">.");
                }

              // Read data from disk
              std::stringstream datastream;

              try
                {
                  datastream << filestream.rdbuf();
                }
              catch (const std::ios::failure &)
                {
                  // broadcast failure state, then throw
                  const int ierr = MPI_Bcast(&invalid_file_size, 1, MPI_UNSIGNED, 0, comm);
                  WBAssertThrow(ierr == 0, "MPI_Bcast failed.");
                  WBAssertThrow (false,
                                 std::string("Could not read file content from <") + filename + ">.");
                }

              data_string = datastream.str();
              WBAssertThrow(static_cast<long int>(data_string.size()) < std::numeric_limits<int>::max(),
                            "File is too large to be send with MPI.");
              filesize = static_cast<int>(data_string.size());

              // Distribute data_size and data across processes
              int ierr = MPI_Bcast(&filesize, 1, MPI_UNSIGNED, 0, comm);
              WBAssertThrow(ierr == 0, "MPI_Bcast failed.");

              // Receive and store data
              ierr = MPI_Bcast(&data_string[0],
                               filesize,
                               MPI_CHAR,
                               0,
                               comm);
              WBAssertThrow(ierr == 0, "MPI_Bcast failed.");
            }
          else
            {
              // Prepare for receiving data
              int filesize = 0;
              int ierr = MPI_Bcast(&filesize, 1, MPI_UNSIGNED, 0, comm);
              WBAssertThrow(ierr == 0, "MPI_Bcast failed.");
              WBAssertThrow(filesize != -1,
                            std::string("Could not open file <") + filename + ">.");

              data_string.resize(static_cast<size_t>(filesize));

              // Receive and store data
              ierr = MPI_Bcast(&data_string[0],
                               filesize,
                               MPI_CHAR,
                               0,
                               comm);
              WBAssertThrow(ierr == 0, "MPI_Bcast failed.");
            }
        }
      else
        {
          std::ifstream filestream;
          filestream.open(filename.c_str());
          if (!filestream)
            {
              WBAssertThrow (false,
                             std::string("Could not open file <") + filename + ">.");
            }
          std::stringstream datastream;
          datastream << filestream.rdbuf();
          data_string = datastream.str();
        }
#else
      std::ifstream filestream;
      filestream.open(filename.c_str());
      if (!filestream)
        {
          WBAssertThrow (false,
                         std::string("Could not open file <") + filename + ">.");
        }
      std::stringstream datastream;
      datastream << filestream.rdbuf();
      data_string = datastream.str();
#endif

      return data_string;
    }

    template std::array<double,2> convert_point_to_array<2>(const Point<2> &point_);
    template std::array<double,3> convert_point_to_array<3>(const Point<3> &point_);


    std::vector<double>
    calculate_ridge_distance_and_spreading(std::vector<std::vector<Point<2>>> mid_oceanic_ridges,
                                           std::vector<std::vector<double>> mid_oceanic_spreading_velocities,
                                           const std::unique_ptr<WorldBuilder::CoordinateSystems::Interface> &coordinate_system,
                                           const Objects::NaturalCoordinate &position_in_natural_coordinates_at_min_depth,
                                           const std::vector<std::vector<double>> &subducting_plate_velocities,
                                           const std::vector<double> &ridge_migration_times)
    {
      const double seconds_in_year = 60.0 * 60.0 * 24.0 * 365.25;  // sec/y

      double distance_ridge = std::numeric_limits<double>::max();
      double spreading_velocity_at_ridge = 0;
      double subducting_velocity_at_trench = 0;
      double ridge_migration_time = 0;

      // first find if the coordinate is on this side of a ridge
      unsigned int relevant_ridge = 0;
      const Point<2> check_point(position_in_natural_coordinates_at_min_depth.get_surface_coordinates(),
                                 position_in_natural_coordinates_at_min_depth.get_coordinate_system());

      Point<2> other_check_point = check_point;
      if (check_point.get_coordinate_system() == CoordinateSystem::spherical)
        other_check_point[0] += check_point[0] < 0 ? 2.0 * WorldBuilder::Consts::PI : -2.0 * WorldBuilder::Consts::PI;

      // if there is only one ridge, there is no transform
      if (mid_oceanic_ridges[0].size() > 1)
        {
          // There are more than one ridge, so there are transform faults
          // Find the first which is on the same side
          for (relevant_ridge = 0; relevant_ridge < mid_oceanic_ridges.size()-1; relevant_ridge++)
            {
              const Point<2> transform_point_0 = mid_oceanic_ridges[relevant_ridge+1][0];
              const Point<2> transform_point_1 = mid_oceanic_ridges[relevant_ridge][mid_oceanic_ridges[relevant_ridge].size()-1];
              const Point<2> reference_point   = mid_oceanic_ridges[relevant_ridge][0];
              const bool reference_on_side_of_line = (transform_point_1[0] - transform_point_0[0])
                                                     * (reference_point[1] - transform_point_0[1])
                                                     - (transform_point_1[1] - transform_point_0[1])
                                                     * (reference_point[0] - transform_point_0[0])
                                                     < 0;
              const bool checkpoint_on_side_of_line = (transform_point_1[0] - transform_point_0[0])
                                                      * (check_point[1] - transform_point_0[1])
                                                      - (transform_point_1[1] - transform_point_0[1])
                                                      * (check_point[0] - transform_point_0[0])
                                                      < 0;


              if (reference_on_side_of_line == checkpoint_on_side_of_line)
                {
                  break;
                }

            }
        }

      for (unsigned int i_coordinate = 0; i_coordinate < mid_oceanic_ridges[relevant_ridge].size() - 1; i_coordinate++)
        {
          const Point<2> segment_point0 = mid_oceanic_ridges[relevant_ridge][i_coordinate];
          const Point<2> segment_point1 = mid_oceanic_ridges[relevant_ridge][i_coordinate + 1];

          const double spreading_velocity_point0 = mid_oceanic_spreading_velocities[relevant_ridge][i_coordinate];
          const double spreading_velocity_point1 = mid_oceanic_spreading_velocities[relevant_ridge][i_coordinate + 1];

          // When subducting_velocities is not input by the user, default value is 0, which
          // results in subducting velocity == spreading_velocity. When a single value is
          // input by the user, subducting velocity != spreading_velocity, but
          // subducting velocity is spatially constant.
          double subducting_velocity_point0 = subducting_plate_velocities[0][0];
          double subducting_velocity_point1 = subducting_plate_velocities[0][0];

          // When subducting_velocities is input as an array, spatial variation
          if (subducting_plate_velocities[0].size() > 1)
            {
              WBAssert(subducting_plate_velocities.size() == mid_oceanic_ridges.size() && \
                       subducting_plate_velocities[relevant_ridge].size() == mid_oceanic_ridges[relevant_ridge].size(),
                       "subducting velocity and ridge coordinates must be the same dimension");
              WBAssert(ridge_migration_times.size() == mid_oceanic_ridges.size(),
                       "the times for ridge migration specified in 'spreading velocity' must be the same dimension "
                       "as ridge coordinates.");
              subducting_velocity_point0 = subducting_plate_velocities[relevant_ridge][i_coordinate];
              subducting_velocity_point1 = subducting_plate_velocities[relevant_ridge][i_coordinate + 1];
              ridge_migration_time = ridge_migration_times[relevant_ridge];
            }

          {
            // based on http://geomalgorithms.com/a02-_lines.html
            const Point<2> v = segment_point1 - segment_point0;
            const Point<2> w1 = check_point - segment_point0;
            const Point<2> w2 = other_check_point - segment_point0;

            const double c1 = (w1[0] * v[0] + w1[1] * v[1]);
            const double c = (v[0] * v[0] + v[1] * v[1]);
            const double c2 = (w2[0] * v[0] + w2[1] * v[1]);


            Point<2> Pb1(coordinate_system->natural_coordinate_system());
            // This part is needed when we want to consider segments instead of lines
            // If you want to have infinite lines, use only the else statement.

            // First, compare the results from the two compare points
            double spreading_velocity_at_ridge_pt1 = 0.0;
            double subducting_velocity_at_trench_pt1 = 0.0;
            double spreading_velocity_at_ridge_pt2 = 0.0;
            double subducting_velocity_at_trench_pt2 = 0.0;

            if (c1 <= 0)
              {
                Pb1=segment_point0;
                spreading_velocity_at_ridge_pt1 = spreading_velocity_point0;
                subducting_velocity_at_trench_pt1 = subducting_velocity_point0;
              }
            else if (c <= c1)
              {
                Pb1=segment_point1;
                spreading_velocity_at_ridge_pt1 = spreading_velocity_point1;
                subducting_velocity_at_trench_pt1 = subducting_velocity_point1;
              }
            else
              {
                Pb1=segment_point0 + (c1 / c) * v;
                spreading_velocity_at_ridge_pt1 = spreading_velocity_point0 + (spreading_velocity_point1 - spreading_velocity_point0) * (c1 / c);
                subducting_velocity_at_trench_pt1 = subducting_velocity_point0 + (subducting_velocity_point1 - subducting_velocity_point0) * (c1 / c);
              }

            Point<2> Pb2(coordinate_system->natural_coordinate_system());
            if (c2 <= 0)
              {
                Pb2=segment_point0;
                spreading_velocity_at_ridge_pt2 = spreading_velocity_point0;
                subducting_velocity_at_trench_pt2 = subducting_velocity_point0;
              }
            else if (c <= c2)
              {
                Pb2=segment_point1;
                spreading_velocity_at_ridge_pt2 = spreading_velocity_point1;
                subducting_velocity_at_trench_pt2 = spreading_velocity_point1;
              }
            else
              {
                Pb2=segment_point0 + (c2 / c) * v;
                spreading_velocity_at_ridge_pt2 = spreading_velocity_point0 + (spreading_velocity_point1 - spreading_velocity_point0) * (c2 / c);
                subducting_velocity_at_trench_pt2 = subducting_velocity_point0 + (subducting_velocity_point1 - subducting_velocity_point0) * (c2 / c);
              }

            Point<3> compare_point1(coordinate_system->natural_coordinate_system());
            Point<3> compare_point2(coordinate_system->natural_coordinate_system());

            compare_point1[0] = coordinate_system->natural_coordinate_system() == cartesian ? Pb1[0] :  position_in_natural_coordinates_at_min_depth.get_depth_coordinate();
            compare_point1[1] = coordinate_system->natural_coordinate_system() == cartesian ? Pb1[1] : Pb1[0];
            compare_point1[2] = coordinate_system->natural_coordinate_system() == cartesian ? position_in_natural_coordinates_at_min_depth.get_depth_coordinate() : Pb1[1];

            compare_point2[0] = coordinate_system->natural_coordinate_system() == cartesian ? Pb2[0] :  position_in_natural_coordinates_at_min_depth.get_depth_coordinate();
            compare_point2[1] = coordinate_system->natural_coordinate_system() == cartesian ? Pb2[1] : Pb2[0];
            compare_point2[2] = coordinate_system->natural_coordinate_system() == cartesian ? position_in_natural_coordinates_at_min_depth.get_depth_coordinate() : Pb2[1];

            const double compare_distance1 = coordinate_system->distance_between_points_at_same_depth(Point<3>(position_in_natural_coordinates_at_min_depth.get_coordinates(),
                                             position_in_natural_coordinates_at_min_depth.get_coordinate_system()),
                                             compare_point1);

            const double compare_distance2 = coordinate_system->distance_between_points_at_same_depth(Point<3>(position_in_natural_coordinates_at_min_depth.get_coordinates(),
                                             position_in_natural_coordinates_at_min_depth.get_coordinate_system()),
                                             compare_point2);

            double compare_distance = compare_distance1;
            double spreading_velocity_at_ridge_pt = spreading_velocity_at_ridge_pt1;
            double subducting_velocity_at_trench_pt = subducting_velocity_at_trench_pt1;

            // This is required in spherical coordinates to ensure that the distance
            // returned is the shortest distance around the sphere.
            if (compare_distance2 < compare_distance1)
              {
                compare_distance = compare_distance2;
                spreading_velocity_at_ridge_pt = spreading_velocity_at_ridge_pt2;
                subducting_velocity_at_trench_pt = subducting_velocity_at_trench_pt2;
              }

            // Then, the distance and velocities are taken from the nearest point on the ridge
            if (i_coordinate == 0 || compare_distance < distance_ridge)
              {
                distance_ridge = compare_distance;
                spreading_velocity_at_ridge = spreading_velocity_at_ridge_pt;
                subducting_velocity_at_trench = subducting_velocity_at_trench_pt;
              }
          }
        }
      std::vector<double> result;
      result.push_back(spreading_velocity_at_ridge / seconds_in_year); // m/s
      result.push_back(distance_ridge);
      result.push_back(subducting_velocity_at_trench / seconds_in_year); // m/s
      result.push_back(ridge_migration_time);
      return result;
    }

    // TODO: implement method for modifying the age of the slab based on ridge/trench migration.
    std::vector<double>
    calculate_effective_trench_and_plate_ages(std::vector<double> ridge_parameters, double distance_along_plane)
    {
      WBAssert(ridge_parameters.size() == 4, "Internal error: ridge_parameters have the wrong size: " << ridge_parameters.size() << " instead of 4.");
      const double seconds_in_year = 60.0 * 60.0 * 24.0 * 365.25;  // sec/y
      const double spreading_velocity = ridge_parameters[0] * seconds_in_year; // m/yr
      const double distance_ridge = ridge_parameters[1];
      const double subducting_velocity = ridge_parameters[2] * seconds_in_year; // m/yr

      WBAssertThrow(subducting_velocity >= 0, "The subducting velocity is less than 0. "
                    "Subducting velocity: " << subducting_velocity);

      // Plate age increases with distance along the slab in the mantle
      double effective_plate_age = (distance_ridge + distance_along_plane) / spreading_velocity; // m/(m/y) = yr

      // Age of trench when the query point was at the trench
      const double age_at_trench = effective_plate_age - distance_along_plane / subducting_velocity; // m/(m/y) = yr
      WBAssertThrow(age_at_trench >= 0, "The age of trench at subducting initiation is less than 0. "
                    "Age at trench: " << age_at_trench);

      std::vector<double> result;
      result.push_back(age_at_trench);
      result.push_back(effective_plate_age);
      return result;

    }

    std::array<std::array<double,3>,3>
    multiply_3x3_matrices(const std::array<std::array<double,3>,3> mat1, const std::array<std::array<double,3>,3> mat2)
    {
      std::array<std::array<double,3>,3> result;
      for (size_t i = 0; i < 3; i++)
        {
          for (size_t j = 0; j < 3; j++)
            {
              result[i][j] = 0;
              for (size_t k = 0; k < 3; k++)
                {
                  result[i][j] += mat1[i][k] * mat2[k][j];
                }
            }
        }

      return result;
    }
  } // namespace Utilities
} // namespace WorldBuilder



