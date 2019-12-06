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

#define CATCH_CONFIG_MAIN

#include <iostream>
#include <memory>

#include <catch2.h>

#include <world_builder/config.h>
#include <world_builder/coordinate_systems/interface.h>

#include <world_builder/features/interface.h>
#include <world_builder/features/continental_plate.h>
#include <world_builder/features/fault_models/temperature/uniform.h>
#include <world_builder/features/fault_models/composition/uniform.h>

#include <world_builder/point.h>

#include <world_builder/types/array.h>
#include <world_builder/types/bool.h>
#include <world_builder/types/double.h>
#include <world_builder/types/plugin_system.h>
#include <world_builder/types/point.h>
#include <world_builder/types/object.h>
#include <world_builder/types/segment.h>
#include <world_builder/types/string.h>
#include <world_builder/types/unsigned_int.h>

#include <world_builder/utilities.h>
#include <world_builder/wrapper_c.h>
using namespace WorldBuilder;

using Catch::Matchers::Contains;

/**
 * Compare the given two std::vector<double> entries with an epsilon (using Catch::Approx)
 */
inline void compare_vectors_approx(
  const std::vector<double> &computed,
  const std::vector<double> &expected)
{
  CHECK(computed.size() == expected.size());
  for (unsigned int i=0; i< computed.size(); ++i)
    {
      INFO("vector index i=" << i << ": ");
      CHECK(computed[i] == Approx(expected[i]));
    }
}

/**
 * Compare the given two std::array<double,3> entries with an epsilon (using Catch::Approx)
 */
inline void compare_3d_arrays_approx(
  const std::array<double,3> &computed,
  const std::array<double,3> &expected)
{
  CHECK(computed.size() == expected.size());
  for (unsigned int i=0; i< computed.size(); ++i)
    {
      INFO("vector index i=" << i << ": ");
      CHECK(computed[i] == Approx(expected[i]));
    }
}


TEST_CASE("WorldBuilder Point: Testing initialize and operators")
{
  // Test initialization of the Point class
  Point<2> p2(cartesian);
  Point<3> p3(cartesian);

  CHECK(p2.get_array() == std::array<double,2> {{0,0}});
  CHECK(p3.get_array() == std::array<double,3> {{0,0,0}});

  const Point<2> p2_array(std::array<double,2> {{1,2}},cartesian);
  const Point<3> p3_array(std::array<double,3> {{1,2,3}},cartesian);

  CHECK(p2_array.get_array() == std::array<double,2> {{1,2}});
  CHECK(p3_array.get_array() == std::array<double,3> {{1,2,3}});

  const Point<2> p2_point(p2_array);
  const Point<3> p3_point(p3_array);

  CHECK(p2_point.get_array() == std::array<double,2> {{1,2}});
  CHECK(p3_point.get_array() == std::array<double,3> {{1,2,3}});

  const Point<2> p2_explicit(3,4,cartesian);
  const Point<3> p3_explicit(4,5,6,cartesian);

  CHECK(p2_explicit.get_array() == std::array<double,2> {{3,4}});
  CHECK(p3_explicit.get_array() == std::array<double,3> {{4,5,6}});


  // Test Point operators

  // Test assign operator
  p2 = p2_array;
  p3 = p3_array;

  CHECK(p2.get_array() == std::array<double,2> {{1,2}});
  CHECK(p3.get_array() == std::array<double,3> {{1,2,3}});

  // Test multiply operator
  p2 = 2 * p2 * 1.0;
  p3 = 2 * p3 * 1.0;

  CHECK(p2.get_array() == std::array<double,2> {{2,4}});
  CHECK(p3.get_array() == std::array<double,3> {{2,4,6}});

  p2 *= 2;
  p3 *= 2;

  CHECK(p2.get_array() == std::array<double,2> {{4,8}});
  CHECK(p3.get_array() == std::array<double,3> {{4,8,12}});

  // Test dot operator
  CHECK(p2_array * p2_explicit == 11);
  CHECK(p3_array * p3_explicit == 32);

  // Test add operator
  p2 = p2 + p2;
  p3 = p3 + p3;

  CHECK(p2.get_array() == std::array<double,2> {{8,16}});
  CHECK(p3.get_array() == std::array<double,3> {{8,16,24}});

  p2 += p2;
  p3 += p3;

  CHECK(p2.get_array() == std::array<double,2> {{16,32}});
  CHECK(p3.get_array() == std::array<double,3> {{16,32,48}});

  // Test subtract operator
  p2 = p2 - (0.5 * p2);
  p3 = p3 - (0.5 * p3);

  CHECK(p2.get_array() == std::array<double,2> {{8,16}});
  CHECK(p3.get_array() == std::array<double,3> {{8,16,24}});

  p2 -=  (0.5 * p2);
  p3 -=  (0.5 * p3);

  CHECK(p2.get_array() == std::array<double,2> {{4,8}});
  CHECK(p3.get_array() == std::array<double,3> {{4,8,12}});

  // Test coordinate system
  //CHECK(p2.get_coordinate_system() == CoordinateSystem::cartesian);
  //CHECK(p3.get_coordinate_system() == CoordinateSystem::cartesian);

  // Test norm and norm_square
  CHECK(p2.norm_square() == 80);
  CHECK(p3.norm_square() == 224);

  CHECK(p2.norm() == std::sqrt(80));
  CHECK(p3.norm() == std::sqrt(224));

  // Test Point utility classes
  std::array<double,2> an2 = Utilities::convert_point_to_array(p2_point);
  std::array<double,3> an3 = Utilities::convert_point_to_array(p3_point);

  CHECK(an2 == std::array<double,2> {{1,2}});
  CHECK(an3 == std::array<double,3> {{1,2,3}});

  CHECK_THROWS_WITH(Point<2>(1,2,3,cartesian),Contains("Can't use the 3d constructor in 2d."));
  CHECK_THROWS_WITH(Point<3>(1,2,cartesian),Contains("Can't use the 2d constructor in 3d."));



}


TEST_CASE("WorldBuilder Utilities: string to conversions")
{
  // Test string to number conversion
  CHECK(Utilities::string_to_double("1") == 1.0);
  CHECK(Utilities::string_to_double(" 1 ") == 1.0);
  CHECK(Utilities::string_to_double(" 1.01 ") == 1.01);

  CHECK_THROWS_WITH(Utilities::string_to_double("1a"),
                    Contains("Could not convert \"1a\" to a double."));
  CHECK_THROWS_WITH(Utilities::string_to_double("a1"),
                    Contains("Could not convert \"a1\" to a double."));
  CHECK_THROWS_WITH(Utilities::string_to_double("a"),
                    Contains("Could not convert \"a\" to a double."));

  CHECK(Utilities::string_to_int("2") == 2);
  CHECK(Utilities::string_to_int(" 2 ") == 2);

  CHECK_THROWS_WITH(Utilities::string_to_int(" 2.02 "),
                    Contains("Could not convert \"2.02\" to an int."));
  CHECK_THROWS_WITH(Utilities::string_to_int("2b"),
                    Contains("Could not convert \"2b\" to an int."));
  CHECK_THROWS_WITH(Utilities::string_to_int("b2"),
                    Contains("Could not convert \"b2\" to an int."));
  CHECK_THROWS_WITH(Utilities::string_to_int("b"),
                    Contains("Could not convert \"b\" to an int."));

  CHECK(Utilities::string_to_unsigned_int("3") == 3);
  CHECK(Utilities::string_to_unsigned_int(" 3 ") == 3);

  CHECK_THROWS_WITH(Utilities::string_to_unsigned_int(" 3.03 "),
                    Contains("Could not convert \"3.03\" to an unsigned int."));
  CHECK_THROWS_WITH(Utilities::string_to_unsigned_int("3c"),
                    Contains("Could not convert \"3c\" to an unsigned int."));
  CHECK_THROWS_WITH(Utilities::string_to_unsigned_int("c3"),
                    Contains("Could not convert \"c3\" to an unsigned int."));
  CHECK_THROWS_WITH(Utilities::string_to_unsigned_int("c"),
                    Contains("Could not convert \"c\" to an unsigned int."));

  // Test point to array conversion
  const Point<2> p2(1,2,cartesian);
  const Point<3> p3(1,2,3,cartesian);

  CHECK(Utilities::convert_point_to_array(p2) == std::array<double,2> {{1,2}});
  CHECK(Utilities::convert_point_to_array(p3) == std::array<double,3> {{1,2,3}});

  // Test coordinate system
  CHECK(Utilities::string_to_coordinate_system("cartesian") == CoordinateSystem::cartesian);
  CHECK(Utilities::string_to_coordinate_system("spherical") == CoordinateSystem::spherical);
  CHECK_THROWS_WITH(Utilities::string_to_coordinate_system("other"), Contains("Coordinate system not implemented."));
}


TEST_CASE("WorldBuilder Utilities: interpolation")
{
  Utilities::interpolation linear;
  std::vector<double> x = {{0,1,2,6}};
  std::vector<double> y = {{10,5,5,35}};
  linear.set_points(x,y,false);
  CHECK(linear(-1) == Approx(15.0));
  CHECK(linear(-0.9) == Approx(14.5));
  CHECK(linear(-0.5) == Approx(12.5));
  CHECK(linear(-0.1) == Approx(10.5));
  CHECK(linear(0) == Approx(10.0));
  CHECK(linear(0.1) == Approx(9.5));
  CHECK(linear(0.5) == Approx(7.5));
  CHECK(linear(0.9) == Approx(5.5));
  CHECK(linear(1) == Approx(5.0));
  CHECK(linear(1.1) == Approx(5.0));
  CHECK(linear(1.5) == Approx(5.0));
  CHECK(linear(1.9) == Approx(5.0));
  CHECK(linear(2) == Approx(5.0));
  CHECK(linear(2.1) == Approx(5.75));
  CHECK(linear(2.5) == Approx(8.75));
  CHECK(linear(2.9) == Approx(11.75));
  CHECK(linear(3) == Approx(12.5));
  CHECK(linear(3.1) == Approx(13.25));
  CHECK(linear(3.5) == Approx(16.25));
  CHECK(linear(3.9) == Approx(19.25));
  CHECK(linear(4) == Approx(20));
  CHECK(linear(5) == Approx(27.5));
  CHECK(linear(6) == Approx(35));
  CHECK(linear(7) == Approx(42.5));

  Utilities::interpolation monotone_cubic_spline;
  monotone_cubic_spline.set_points(x,y,true);

  CHECK(monotone_cubic_spline(-1) == Approx(-5));
  CHECK(monotone_cubic_spline(-0.9) == Approx(-2.15));
  CHECK(monotone_cubic_spline(-0.7) == Approx(2.65));
  CHECK(monotone_cubic_spline(-0.5) == Approx(6.25));
  CHECK(monotone_cubic_spline(-0.3) == Approx(8.65));
  CHECK(monotone_cubic_spline(-0.1) == Approx(9.85));
  CHECK(monotone_cubic_spline(0) == Approx(10.0));
  CHECK(monotone_cubic_spline(0.1) == Approx(9.86));
  CHECK(monotone_cubic_spline(0.3) == Approx(8.92));
  CHECK(monotone_cubic_spline(0.5) == Approx(7.5));
  CHECK(monotone_cubic_spline(0.7) == Approx(6.08));
  CHECK(monotone_cubic_spline(0.9) == Approx(5.14));
  CHECK(monotone_cubic_spline(1) == Approx(5.0));
  CHECK(monotone_cubic_spline(1.1) == Approx(5.0));
  CHECK(monotone_cubic_spline(1.3) == Approx(5.0));
  CHECK(monotone_cubic_spline(1.5) == Approx(5.0));
  CHECK(monotone_cubic_spline(1.7) == Approx(5.0));
  CHECK(monotone_cubic_spline(1.9) == Approx(5.0));
  CHECK(monotone_cubic_spline(2) == Approx(5.0));
  CHECK(monotone_cubic_spline(2.1) == Approx(5.03703125));
  CHECK(monotone_cubic_spline(2.3) == Approx(5.32484375));
  CHECK(monotone_cubic_spline(2.5) == Approx(5.87890625));
  CHECK(monotone_cubic_spline(2.7) == Approx(6.67671875));
  CHECK(monotone_cubic_spline(2.9) == Approx(7.69578125));
  CHECK(monotone_cubic_spline(3) == Approx(8.28125));
  CHECK(monotone_cubic_spline(3.1) == Approx(8.91359375));
  CHECK(monotone_cubic_spline(3.3) == Approx(10.30765625));
  CHECK(monotone_cubic_spline(3.5) == Approx(11.85546875));
  CHECK(monotone_cubic_spline(3.7) == Approx(13.53453125));
  CHECK(monotone_cubic_spline(3.9) == Approx(15.32234375));
  CHECK(monotone_cubic_spline(4) == Approx(16.25));
  CHECK(monotone_cubic_spline(4.5) == Approx(21.11328125));
  CHECK(monotone_cubic_spline(5) == Approx(26.09375));
  CHECK(monotone_cubic_spline(5.5) == Approx(30.83984375));
  CHECK(monotone_cubic_spline(6) == Approx(35));
  CHECK(monotone_cubic_spline(6.5) == Approx(38.75));
  CHECK(monotone_cubic_spline(7) == Approx(42.5));

  Utilities::interpolation monotone_cubic_spline2;
  y[1] = -5;
  y[3] = -35;
  monotone_cubic_spline2.set_points(x,y,true);
  CHECK(monotone_cubic_spline2(-1) == Approx(-35));
  CHECK(monotone_cubic_spline2(-0.5) == Approx(-1.25));
  CHECK(monotone_cubic_spline2(0) == Approx(10));
  CHECK(monotone_cubic_spline2(0.5) == Approx(2.5));
  CHECK(monotone_cubic_spline2(1) == Approx(-5.0));
  CHECK(monotone_cubic_spline2(1.5) == Approx(0.0));
  CHECK(monotone_cubic_spline2(2) == Approx(5.0));
  CHECK(monotone_cubic_spline2(2.5) == Approx(3.828125));
  CHECK(monotone_cubic_spline2(3) == Approx(0.625));
  CHECK(monotone_cubic_spline2(3.5) == Approx(-4.140625));
  CHECK(monotone_cubic_spline2(4) == Approx(-10));
  CHECK(monotone_cubic_spline2(4.5) == Approx(-16.484375));
  CHECK(monotone_cubic_spline2(5) == Approx(-23.125));
  CHECK(monotone_cubic_spline2(5.5) == Approx(-29.453125));
  CHECK(monotone_cubic_spline2(6) == Approx(-35.0));
  CHECK(monotone_cubic_spline2(6.5) == Approx(-40));
  CHECK(monotone_cubic_spline2(7) == Approx(-45));

  Utilities::interpolation monotone_cubic_spline3;
  y[0] = 10;
  y[1] = -5;
  y[2] = -10;
  y[3] = -35;
  monotone_cubic_spline3.set_points(x,y,true);
  CHECK(monotone_cubic_spline3(-1) == Approx(-27.5));
  CHECK(monotone_cubic_spline3(-0.5) == Approx(0.625));
  CHECK(monotone_cubic_spline3(0) == Approx(10.0));
  CHECK(monotone_cubic_spline3(0.5) == Approx(3.4375));
  CHECK(monotone_cubic_spline3(1) == Approx(-5.0));
  CHECK(monotone_cubic_spline3(1.5) == Approx(-7.7272727273));
  CHECK(monotone_cubic_spline3(2) == Approx(-10.0));
  CHECK(monotone_cubic_spline3(2.5) == Approx(-12.9074928977));
  CHECK(monotone_cubic_spline3(3) == Approx(-15.9303977273));
  CHECK(monotone_cubic_spline3(3.5) == Approx(-19.0420809659));
  CHECK(monotone_cubic_spline3(4) == Approx(-22.2159090909));
  CHECK(monotone_cubic_spline3(4.5) == Approx(-25.4252485795));
  CHECK(monotone_cubic_spline3(5) == Approx(-28.6434659091));
  CHECK(monotone_cubic_spline3(5.5) == Approx(-31.8439275568));
  CHECK(monotone_cubic_spline3(6) == Approx(-35.0));
  CHECK(monotone_cubic_spline3(6.5) == Approx(-38.125));
  CHECK(monotone_cubic_spline3(7) == Approx(-41.25));

  // bi monotone cubic spline
  Utilities::interpolation monotone_cubic_spline_x;
  Utilities::interpolation monotone_cubic_spline_y;
  x[0] = 0;
  x[1] = 1;
  x[2] = 2;
  x[3] = 3;
  y[0] = 10;
  y[1] = 10;
  y[2] = 5;
  y[3] = 0;
  monotone_cubic_spline_x.set_points(x,y,true);
  y[0] = 0;
  y[1] = 5;
  y[2] = 10;
  y[3] = 10;
  monotone_cubic_spline_y.set_points(x,y,true);

  CHECK(monotone_cubic_spline_x(0) == Approx(10));
  CHECK(monotone_cubic_spline_x(0.25) == Approx(10));
  CHECK(monotone_cubic_spline_x(0.5) == Approx(10));
  CHECK(monotone_cubic_spline_x(0.75) == Approx(10));
  CHECK(monotone_cubic_spline_x(1) == Approx(10.0));
  CHECK(monotone_cubic_spline_x(1.25) == Approx(9.453125));
  CHECK(monotone_cubic_spline_x(1.5) == Approx(8.125));
  CHECK(monotone_cubic_spline_x(1.75) == Approx(6.484375));
  CHECK(monotone_cubic_spline_x(2) == Approx(5.0));
  CHECK(monotone_cubic_spline_x(2.25) == Approx(3.75));
  CHECK(monotone_cubic_spline_x(2.5) == Approx(2.5));
  CHECK(monotone_cubic_spline_x(2.75) == Approx(1.25));
  CHECK(monotone_cubic_spline_x(3) == Approx(0));

  CHECK(monotone_cubic_spline_y(0) == Approx(0));
  CHECK(monotone_cubic_spline_y(0.25) == Approx(0.546875));
  CHECK(monotone_cubic_spline_y(0.5) == Approx(1.875));
  CHECK(monotone_cubic_spline_y(0.75) == Approx(3.515625));
  CHECK(monotone_cubic_spline_y(1) == Approx(5.0));
  CHECK(monotone_cubic_spline_y(1.25) == Approx(6.484375));
  CHECK(monotone_cubic_spline_y(1.5) == Approx(8.125));
  CHECK(monotone_cubic_spline_y(1.75) == Approx(9.453125));
  CHECK(monotone_cubic_spline_y(2) == Approx(10.0));
  CHECK(monotone_cubic_spline_y(2.25) == Approx(10.0));
  CHECK(monotone_cubic_spline_y(2.5) == Approx(10.0));
  CHECK(monotone_cubic_spline_y(2.75) == Approx(10.0));
  CHECK(monotone_cubic_spline_y(3) == Approx(10.0));
}

TEST_CASE("WorldBuilder Utilities: Point in polygon")
{
  std::vector<Point<2> > point_list_4_elements(4, Point<2>(cartesian));
  point_list_4_elements[0] = Point<2>(0,0,cartesian);
  point_list_4_elements[1] = Point<2>(5,0,cartesian);
  point_list_4_elements[2] = Point<2>(5,5,cartesian);
  point_list_4_elements[3] = Point<2>(0,5,cartesian);

  std::vector<Point<2> > point_list_3_elements(3, Point<2>(cartesian));
  point_list_3_elements[0] = Point<2>(10,10,cartesian);
  point_list_3_elements[1] = Point<2>(10,15,cartesian);
  point_list_3_elements[2] = Point<2>(15,15,cartesian);

  std::vector<Point<2> > check_points(9, Point<2>(cartesian));
  check_points[0] = Point<2>(-1,-1,cartesian);
  check_points[1] = Point<2>(0,0,cartesian);
  check_points[2] = Point<2>(0,5,cartesian);
  check_points[3] = Point<2>(5,0,cartesian);
  check_points[4] = Point<2>(5,5,cartesian);
  check_points[5] = Point<2>(5,5.01,cartesian);
  check_points[6] = Point<2>(1,1,cartesian);
  check_points[7] = Point<2>(12.5,12,cartesian);
  check_points[8] = Point<2>(11.5,12,cartesian);

  std::vector<std::array<bool,2> > awnsers(9);
  awnsers[0] = {{false,false}};
  awnsers[1] = {{true,false}};
  awnsers[2] = {{true,false}};
  awnsers[3] = {{true,false}};
  awnsers[4] = {{true,false}};
  awnsers[5] = {{false,false}};
  awnsers[6] = {{true,false}};
  awnsers[7] = {{false,false}};
  awnsers[8] = {{false,true}};

  std::vector<std::array<double,2> > awnsers_signed_distance(9);
  awnsers_signed_distance[0] = {{-std::sqrt(2), -std::sqrt(11 * 11 + 11 * 11)}};
  awnsers_signed_distance[1] = {{0,-std::sqrt(10 * 10 + 10 * 10)}};
  awnsers_signed_distance[2] = {{0,-std::sqrt(125)}};
  awnsers_signed_distance[3] = {{0,-std::sqrt(125)}};
  awnsers_signed_distance[4] = {{0,-std::sqrt(50)}};
  awnsers_signed_distance[5] = {{-std::sqrt(0.01 * 0.01),-std::sqrt(5 * 5 + 4.99 * 4.99)}};
  awnsers_signed_distance[6] = {{1,-std::sqrt(9 * 9 + 9 * 9)}};
  awnsers_signed_distance[7] = {{-10.2591422643,-0.3535533906}};
  awnsers_signed_distance[8] = {{-9.5524865873,0.3535533906}};

  for (unsigned int i = 0; i < check_points.size(); ++i)
    {
      INFO("checking point " << i << " = (" << check_points[i][0] << ":" << check_points[i][1] << ")");
      CHECK(Utilities::polygon_contains_point(point_list_4_elements,check_points[i]) == awnsers[i][0]);
      CHECK(Utilities::polygon_contains_point(point_list_3_elements,check_points[i]) == awnsers[i][1]);
      CHECK(Utilities::signed_distance_to_polygon(point_list_4_elements,check_points[i]) == Approx(awnsers_signed_distance[i][0]));
      CHECK(Utilities::signed_distance_to_polygon(point_list_3_elements,check_points[i]) == Approx(awnsers_signed_distance[i][1]));
    }

  std::vector<Point<2> > point_list_2_elements(2, Point<2>(cartesian));
  CHECK_THROWS_WITH(Utilities::signed_distance_to_polygon(point_list_2_elements,check_points[0]),
                    Contains("Not enough polygon points were specified."));

  std::vector<Point<2> > point_list_1_elements(1, Point<2>(cartesian));
  CHECK_THROWS_WITH(Utilities::signed_distance_to_polygon(point_list_1_elements,check_points[0]),
                    Contains("Not enough polygon points were specified."));

  std::vector<Point<2> > point_list_0_elements(0, Point<2>(cartesian));
  CHECK_THROWS_WITH(Utilities::signed_distance_to_polygon(point_list_0_elements,check_points[0]),
                    Contains("Not enough polygon points were specified."));
}


TEST_CASE("WorldBuilder Utilities: Natural Coordinate")
{
  // Cartesian
  unique_ptr<CoordinateSystems::Interface> cartesian(CoordinateSystems::Interface::create("cartesian",NULL));

  // Test the natural coordinate system
  Utilities::NaturalCoordinate nca1(std::array<double,3> {{1,2,3}},*cartesian);
  CHECK(nca1.get_coordinates() == std::array<double,3> {{1,2,3}});
  CHECK(nca1.get_surface_coordinates() == std::array<double,2> {{1,2}});
  CHECK(nca1.get_depth_coordinate() == 3);

  Utilities::NaturalCoordinate ncp1(Point<3>(1,2,3,CoordinateSystem::cartesian),*cartesian);
  CHECK(ncp1.get_coordinates() == std::array<double,3> {{1,2,3}});
  CHECK(ncp1.get_surface_coordinates() == std::array<double,2> {{1,2}});
  CHECK(ncp1.get_depth_coordinate() == 3);


  unique_ptr<CoordinateSystems::Interface> spherical(CoordinateSystems::Interface::create("spherical",NULL));

  // Test the natural coordinate system
  Utilities::NaturalCoordinate nsa1(std::array<double,3> {{1,2,3}},*spherical);
  std::array<double,3> nsa1_array = nsa1.get_coordinates();
  CHECK(nsa1_array[0] == Approx(std::sqrt(1.0 * 1.0 + 2.0 * 2.0 + 3.0 * 3.0)));
  CHECK(nsa1_array[1] == Approx(1.1071487178));
  CHECK(nsa1_array[2] == Approx(0.9302740141));
  std::array<double,2> nsa1_surface_array = nsa1.get_surface_coordinates();
  CHECK(nsa1_surface_array[0] == Approx(1.1071487178));
  CHECK(nsa1_surface_array[1] == Approx(0.9302740141));
  CHECK(nsa1.get_depth_coordinate() == Approx(std::sqrt(1.0 * 1.0 + 2.0 * 2.0 + 3.0 * 3.0)));


  Utilities::NaturalCoordinate nsp1(Point<3>(1,2,3,CoordinateSystem::spherical),*spherical);
  std::array<double,3> nsp1_array = nsp1.get_coordinates();
  CHECK(nsp1_array[0] == Approx(std::sqrt(1.0 * 1.0 + 2.0 * 2.0 + 3.0 * 3.0)));
  CHECK(nsp1_array[1] == Approx(1.1071487178));
  CHECK(nsp1_array[2] == Approx(0.9302740141));
  std::array<double,2> nsp1_surface_array = nsp1.get_surface_coordinates();
  CHECK(nsp1_surface_array[0] == Approx(1.1071487178));
  CHECK(nsp1_surface_array[1] == Approx(0.9302740141));
  CHECK(nsp1.get_depth_coordinate() == Approx(std::sqrt(1.0 * 1.0 + 2.0 * 2.0 + 3.0 * 3.0)));

}

TEST_CASE("WorldBuilder Utilities: Coordinate systems transformations")
{
  // Test coordinate system transformation
  {
    Point<3> cartesian(3,4,5,CoordinateSystem::cartesian);

    Point<3> spherical(Utilities::cartesian_to_spherical_coordinates(Point<3>(cartesian.get_array(),CoordinateSystem::cartesian)), CoordinateSystem::spherical);

    compare_vectors_approx(std::vector<double>(std::begin(spherical.get_array()), std::end(spherical.get_array())),
    std::vector<double> {{std::sqrt(3*3+4*4+5*5),0.927295218001613,0.7853982}});

    Point<3> cartesian_back(Utilities::spherical_to_cartesian_coordinates(spherical.get_array()), CoordinateSystem::cartesian);

    compare_vectors_approx(std::vector<double>(std::begin(cartesian_back.get_array()), std::end(cartesian_back.get_array())),
    std::vector<double> {{3,4,5}});
  }

  {
    Point<3> cartesian(-2,-1,6,CoordinateSystem::cartesian);

    Point<3> spherical(Utilities::cartesian_to_spherical_coordinates(Point<3>(cartesian.get_array(),CoordinateSystem::cartesian)), CoordinateSystem::spherical);

    compare_vectors_approx(std::vector<double>(std::begin(spherical.get_array()), std::end(spherical.get_array())),
    std::vector<double> {{std::sqrt(2*2+1*1+6*6),-2.6779450446,1.2140629383}});

    Point<3> cartesian_back(Utilities::spherical_to_cartesian_coordinates(spherical.get_array()), CoordinateSystem::cartesian);

    compare_vectors_approx(std::vector<double>(std::begin(cartesian_back.get_array()), std::end(cartesian_back.get_array())),
    std::vector<double> {{-2,-1,6}});
  }


}

TEST_CASE("WorldBuilder Utilities: cross product")
{
  const Point<3> unit_x(1,0,0,cartesian);
  const Point<3> unit_y(0,1,0,cartesian);
  const Point<3> unit_z(0,0,1,cartesian);

  compare_3d_arrays_approx(Utilities::cross_product(unit_x, unit_x).get_array(), std::array<double,3> {{0,0,0}});
  compare_3d_arrays_approx(Utilities::cross_product(unit_x, unit_y).get_array(), std::array<double,3> {{0,0,1}});
  compare_3d_arrays_approx(Utilities::cross_product(unit_x, unit_z).get_array(), std::array<double,3> {{0,-1,0}});

  compare_3d_arrays_approx(Utilities::cross_product(unit_y, unit_x).get_array(), std::array<double,3> {{0,0,-1}});
  compare_3d_arrays_approx(Utilities::cross_product(unit_y, unit_y).get_array(), std::array<double,3> {{0,0,0}});
  compare_3d_arrays_approx(Utilities::cross_product(unit_y, unit_z).get_array(), std::array<double,3> {{1,0,0}});

  compare_3d_arrays_approx(Utilities::cross_product(unit_z, unit_x).get_array(), std::array<double,3> {{0,1,0}});
  compare_3d_arrays_approx(Utilities::cross_product(unit_z, unit_y).get_array(), std::array<double,3> {{-1,0,0}});
  compare_3d_arrays_approx(Utilities::cross_product(unit_z, unit_z).get_array(), std::array<double,3> {{0,0,0}});


  const double sqrt2 = sqrt(0.5);
  const Point<3> sqrt2_x(sqrt2,0,0,cartesian);
  const Point<3> sqrt2_y(0,sqrt2,0,cartesian);
  const Point<3> sqrt2_z(0,0,sqrt2,cartesian);

  const Point<3> unit_xy(sqrt2,sqrt2,0,cartesian);
  const Point<3> unit_xz(sqrt2,0,sqrt2,cartesian);
  const Point<3> unit_yz(0,sqrt2,sqrt2,cartesian);

  compare_3d_arrays_approx(Utilities::cross_product(unit_xy, sqrt2_x).get_array(), std::array<double,3> {{0,0,-0.5}});
  compare_3d_arrays_approx(Utilities::cross_product(unit_xy, sqrt2_y).get_array(), std::array<double,3> {{0,0,0.5}});
  compare_3d_arrays_approx(Utilities::cross_product(unit_xy, sqrt2_z).get_array(), std::array<double,3> {{0.5,-0.5,0}});

  compare_3d_arrays_approx(Utilities::cross_product(unit_xz, sqrt2_x).get_array(), std::array<double,3> {{0,0.5,0}});
  compare_3d_arrays_approx(Utilities::cross_product(unit_xz, sqrt2_y).get_array(), std::array<double,3> {{-0.5,0,0.5}});
  compare_3d_arrays_approx(Utilities::cross_product(unit_xz, sqrt2_z).get_array(), std::array<double,3> {{0,-0.5,0}});

  compare_3d_arrays_approx(Utilities::cross_product(unit_yz, sqrt2_x).get_array(), std::array<double,3> {{0,0.5,-0.5}});
  compare_3d_arrays_approx(Utilities::cross_product(unit_yz, sqrt2_y).get_array(), std::array<double,3> {{-0.5,0,0}});
  compare_3d_arrays_approx(Utilities::cross_product(unit_yz, sqrt2_z).get_array(), std::array<double,3> {{0.5,0,0}});

  const Point<3> point1(2,3,4,cartesian);
  const Point<3> point2(5,6,7,cartesian);

  compare_3d_arrays_approx(Utilities::cross_product(point1, point2).get_array(), std::array<double,3> {{-3,6,-3}});
  compare_3d_arrays_approx(Utilities::cross_product(point2, point1).get_array(), std::array<double,3> {{3,-6,3}});
}

TEST_CASE("WorldBuilder C wrapper")
{
  // First test a world builder file with a cross section defined
  std::string file = WorldBuilder::Data::WORLD_BUILDER_SOURCE_DIR + "/tests/data/simple_wb1.json";
  void *ptr_world = NULL;
  void **ptr_ptr_world = &ptr_world;
  const char *world_builder_file = file.c_str();

  create_world(ptr_ptr_world, world_builder_file);

  double temperature = 0;

  temperature_2d(*ptr_ptr_world, 1, 2, 0, 10, &temperature);
  CHECK(temperature == Approx(1600));
  temperature_3d(*ptr_ptr_world, 1, 2, 3, 0, 10, &temperature);
  CHECK(temperature == Approx(1600));
  temperature_2d(*ptr_ptr_world, 550e3, 0, 0, 10, &temperature);
  CHECK(temperature == Approx(150));
  temperature_3d(*ptr_ptr_world, 120e3, 500e3, 0, 0, 10, &temperature);
  CHECK(temperature == Approx(150));

  // Test the compositions
  double composition = 0.0;

  composition_2d(*ptr_ptr_world, 1, 2, 0, 2, &composition);
  CHECK(composition == 0.0);
  composition_3d(*ptr_ptr_world, 1, 2, 3, 0, 2, &composition);
  CHECK(composition == 0.0);
  composition_2d(*ptr_ptr_world,  550e3, 0, 0, 3, &composition);
  CHECK(composition == 1.0);
  composition_3d(*ptr_ptr_world, 120e3, 500e3, 0, 0, 3, &composition);
  CHECK(composition == 1.0);

  release_world(*ptr_ptr_world);

  // Now test a world builder file without a cross section defined
  file = WorldBuilder::Data::WORLD_BUILDER_SOURCE_DIR + "/tests/data/simple_wb2.json";
  ptr_world = NULL;
  ptr_ptr_world = &ptr_world;
  const char *world_builder_file2 = file.c_str();

  create_world(ptr_ptr_world, world_builder_file2);


  CHECK_THROWS_WITH(temperature_2d(*ptr_ptr_world, 1, 2, 0, 10, &temperature),
                    Contains("This function can only be called when the cross section "
                             "variable in the world builder file has been set. Dim is 3."));
  temperature_3d(*ptr_ptr_world, 1, 2, 3, 0, 10, &temperature);
  CHECK(temperature == Approx(1600));
  temperature_3d(*ptr_ptr_world, 120e3, 500e3, 0, 0, 10, &temperature);
  CHECK(temperature == Approx(150));

  // Test the compositions
  CHECK_THROWS_WITH(composition_2d(*ptr_ptr_world, 1, 2, 0, 2, &composition),
                    Contains("This function can only be called when the cross section "
                             "variable in the world builder file has been set. Dim is 3."));

  composition_3d(*ptr_ptr_world, 1, 2, 3, 0, 2, &composition);
  CHECK(composition == 0.0);
  composition_3d(*ptr_ptr_world, 120e3, 500e3, 0, 0, 3, &composition);
  CHECK(composition == 1.0);

  release_world(*ptr_ptr_world);
}


TEST_CASE("WorldBuilder Coordinate Systems: Interface")
{
  std::string file_name = WorldBuilder::Data::WORLD_BUILDER_SOURCE_DIR + "/tests/data/oceanic_plate_spherical.wb";
  WorldBuilder::World world(file_name);

  CHECK_THROWS_WITH(CoordinateSystems::Interface::create("!not_implemented_coordinate_system!",&world),
                    Contains("Internal error: Plugin with name '!not_implemented_coordinate_system!' is not found. "
                             "The size of factories is "));

  unique_ptr<CoordinateSystems::Interface> interface(CoordinateSystems::Interface::create("cartesian",&world));

  interface->declare_entries(world.parameters, "", {});

  CHECK(interface->cartesian_to_natural_coordinates(std::array<double,3> {{1,2,3}}) == std::array<double,3> {{1,2,3}});
  CHECK(interface->natural_to_cartesian_coordinates(std::array<double,3> {{1,2,3}}) == std::array<double,3> {{1,2,3}});

  CHECK(interface->natural_coordinate_system() == CoordinateSystem::cartesian);

}

TEST_CASE("WorldBuilder Coordinate Systems: Cartesian")
{
  unique_ptr<CoordinateSystems::Interface> cartesian(CoordinateSystems::Interface::create("cartesian",NULL));

  //todo:fix
  //cartesian->declare_entries();

  CHECK(cartesian->cartesian_to_natural_coordinates(std::array<double,3> {{1,2,3}}) == std::array<double,3> {{1,2,3}});
  CHECK(cartesian->natural_to_cartesian_coordinates(std::array<double,3> {{1,2,3}}) == std::array<double,3> {{1,2,3}});

  CHECK(cartesian->natural_coordinate_system() == CoordinateSystem::cartesian);

  // distance between two points at the same depth
  Point<3> point_1(0.0,0.0,10.0, CoordinateSystem::cartesian);
  Point<3> point_2(1.0,2.0,10.0, CoordinateSystem::cartesian);
  Point<3> point_3(3.0,2.0,10.0, CoordinateSystem::cartesian);
  Point<3> point_4(3.0,3.0,10.0, CoordinateSystem::cartesian);

  CHECK(cartesian->distance_between_points_at_same_depth(point_1, point_2) == Approx(std::sqrt(1 + 2 * 2)));
  CHECK(cartesian->distance_between_points_at_same_depth(point_2, point_3) == Approx(2.0));
  CHECK(cartesian->distance_between_points_at_same_depth(point_2, point_4) == Approx(std::sqrt(2 * 2 + 1)));

}

TEST_CASE("WorldBuilder Coordinate Systems: Spherical")
{
  // TODO: make test where a cartesian wb file is loaded into a spherical coordinate system.
  std::string file_name = WorldBuilder::Data::WORLD_BUILDER_SOURCE_DIR + "/tests/data/oceanic_plate_spherical.wb";

  WorldBuilder::World world(file_name);

  unique_ptr<CoordinateSystems::Interface> spherical(CoordinateSystems::Interface::create("spherical", &world));

  world.parameters.enter_subsection("coordinate system");
  {
    world.parameters.enter_subsection("spherical");
    {
      spherical->declare_entries(world.parameters,"", {});
    }
    world.parameters.leave_subsection();
  }
  world.parameters.leave_subsection();

  std::array<double,3> spherical_array = spherical->cartesian_to_natural_coordinates(std::array<double,3> {{1,2,3}});
  CHECK(spherical_array[0] == Approx(std::sqrt(1.0 * 1.0 + 2.0 * 2.0 + 3.0 * 3.0)));
  CHECK(spherical_array[1] == Approx(1.1071487178));
  CHECK(spherical_array[2] == Approx(0.9302740141));
  std::array<double,3> cartesian_array = spherical->natural_to_cartesian_coordinates(std::array<double,3> {{std::sqrt(1.0 * 1.0 + 2.0 * 2.0 + 3.0 * 3.0),1.1071487178,0.9302740141}});
  CHECK(cartesian_array[0] == Approx(1));
  CHECK(cartesian_array[1] == Approx(2));
  CHECK(cartesian_array[2] == Approx(3));

  CHECK(spherical->natural_coordinate_system() == CoordinateSystem::spherical);

  // distance between two points at the same depth
  double dtr = Utilities::const_pi / 180.0;
  // first check unit radius, this the central angle
  Point<3> unit_point_1(1.0, 0.0 * dtr, 0.0 * dtr, CoordinateSystem::spherical);
  Point<3> unit_point_2(1.0, 1.0 * dtr, 0.0 * dtr, CoordinateSystem::spherical);
  Point<3> unit_point_3(1.0, 0.0 * dtr, 1.0 * dtr, CoordinateSystem::spherical);
  Point<3> unit_point_4(1.0, 1.0 * dtr, 1.0 * dtr, CoordinateSystem::spherical);
  Point<3> unit_point_5(1.0, 90.0 * dtr, 90.0 * dtr, CoordinateSystem::spherical);
  Point<3> unit_point_6(1.0, -90.0 * dtr, 0.0 * dtr, CoordinateSystem::spherical);
  Point<3> unit_point_7(1.0, 90.0 * dtr, 180.0 * dtr, CoordinateSystem::spherical);

  CHECK(spherical->distance_between_points_at_same_depth(unit_point_1, unit_point_2) == Approx(dtr));
  CHECK(spherical->distance_between_points_at_same_depth(unit_point_1, unit_point_3) == Approx(dtr));
  CHECK(spherical->distance_between_points_at_same_depth(unit_point_1, unit_point_4) ==
        Approx(std::acos(std::sin(0) * std::sin(1*dtr) +
                         std::cos(0) * std::cos(1*dtr) * std::cos(1*dtr))));
  CHECK(spherical->distance_between_points_at_same_depth(unit_point_1, unit_point_5) == Approx(0.5 * Utilities::const_pi));
  CHECK(spherical->distance_between_points_at_same_depth(unit_point_6, unit_point_7) == Approx(Utilities::const_pi));

  // secondly check non-unit radius
  Point<3> point_1(10.0, 0.0 * dtr, 0.0 * dtr, CoordinateSystem::spherical);
  Point<3> point_2(10.0, 1.0 * dtr, 0.0 * dtr, CoordinateSystem::spherical);
  Point<3> point_3(10.0, 0.0 * dtr, 1.0 * dtr, CoordinateSystem::spherical);
  Point<3> point_4(10.0, 1.0 * dtr, 1.0 * dtr, CoordinateSystem::spherical);
  Point<3> point_5(10.0, 90.0 * dtr, 90.0 * dtr, CoordinateSystem::spherical);
  Point<3> point_6(10.0, -90.0 * dtr, 0.0 * dtr, CoordinateSystem::spherical);
  Point<3> point_7(10.0, 90.0 * dtr, 180.0 * dtr, CoordinateSystem::spherical);

  CHECK(spherical->distance_between_points_at_same_depth(point_1, point_2) == Approx(10 * dtr));
  CHECK(spherical->distance_between_points_at_same_depth(point_1, point_3) == Approx(10 * dtr));
  CHECK(spherical->distance_between_points_at_same_depth(point_1, point_4) ==
        Approx(10 * std::acos(std::sin(0) * std::sin(1*dtr) +
                              std::cos(0) * std::cos(1*dtr) * std::cos(1*dtr))));
  CHECK(spherical->distance_between_points_at_same_depth(point_1, point_5) == Approx(10 * 0.5 * Utilities::const_pi));
  CHECK(spherical->distance_between_points_at_same_depth(point_6, point_7) == Approx(10 * Utilities::const_pi));

}

TEST_CASE("WorldBuilder Features: Interface")
{
  std::string file_name = WorldBuilder::Data::WORLD_BUILDER_SOURCE_DIR + "/tests/data/simple_wb1.json";

  WorldBuilder::World world(file_name);
  CHECK_THROWS_WITH(Features::Interface::create("!not_implemented_feature!", &world),
                    Contains("Internal error: Plugin with name '!not_implemented_feature!' is not found. "
                             "The size of factories is "));

  std::unique_ptr<Features::Interface> interface = Features::Interface::create("continental plate", &world);

}

TEST_CASE("WorldBuilder Features: Continental Plate")
{
  std::string file_name = WorldBuilder::Data::WORLD_BUILDER_SOURCE_DIR + "/tests/data/continental_plate.wb";
  WorldBuilder::World world1(file_name);

  // Check continental plate directly
  std::unique_ptr<Features::Interface> continental_plate = Features::Interface::create("continental plate", &world1);

  // Check continental plate through the world
  std::array<double,3> position = {{0,0,0}};
  CHECK(world1.temperature(position, 0, 10) == Approx(1600));

  position = {{250e3,500e3,0}};
  CHECK(world1.temperature(position, 0, 10) == Approx(150));
  CHECK(world1.temperature(position, 74e3, 10) == Approx(150));
  CHECK(world1.temperature(position, 76e3, 10) == Approx(100));
  CHECK(world1.temperature(position, 149e3, 10) == Approx(100));
  CHECK(world1.temperature(position, 151e3, 10) == Approx(50));
  CHECK(world1.temperature(position, 224e3, 10) == Approx(50));
  CHECK(world1.temperature(position, 226e3, 10) == Approx(1704.5201415988));
  CHECK(world1.temperature(position, 240e3, 10) == Approx(1711.2149738521));
  CHECK(world1.temperature(position, 260e3, 10) == Approx(1720.8246597128));

  CHECK(world1.composition(position, 0, 0) == 0.0);
  CHECK(world1.composition(position, 0, 1) == 0.0);
  CHECK(world1.composition(position, 0, 2) == 0.0);
  CHECK(world1.composition(position, 0, 3) == 1.0);
  CHECK(world1.composition(position, 0, 4) == 0.0);
  CHECK(world1.composition(position, 0, 5) == 0.0);
  CHECK(world1.composition(position, 240e3, 0) == 0.0);
  CHECK(world1.composition(position, 240e3, 1) == 0.0);
  CHECK(world1.composition(position, 240e3, 2) == 0.0);
  CHECK(world1.composition(position, 240e3, 3) == 1.0);
  CHECK(world1.composition(position, 240e3, 4) == 0.0);
  CHECK(world1.composition(position, 240e3, 5) == 0.0);
  CHECK(world1.composition(position, 260e3, 0) == 0.0);
  CHECK(world1.composition(position, 260e3, 1) == 0.0);
  CHECK(world1.composition(position, 260e3, 2) == 0.0);
  CHECK(world1.composition(position, 260e3, 3) == 0.0);
  CHECK(world1.composition(position, 260e3, 4) == 0.0);
  CHECK(world1.composition(position, 260e3, 5) == 0.0);

  position = {{1500e3,1500e3,0}};
  CHECK(world1.temperature(position, 0, 10) == Approx(20));
  CHECK(world1.temperature(position, 240e3, 10) == Approx(21.3901871732));
  CHECK(world1.temperature(position, 260e3, 10) == Approx(1720.8246597128));

  CHECK(world1.composition(position, 0, 0) == 0.0);
  CHECK(world1.composition(position, 0, 1) == 0.0);
  CHECK(world1.composition(position, 0, 2) == 1.0);
  CHECK(world1.composition(position, 0, 3) == 0.0);
  CHECK(world1.composition(position, 0, 4) == 0.0);
  CHECK(world1.composition(position, 0, 5) == 0.0);
  CHECK(world1.composition(position, 240e3, 0) == 0.0);
  CHECK(world1.composition(position, 240e3, 1) == 0.0);
  CHECK(world1.composition(position, 240e3, 2) == 1.0);
  CHECK(world1.composition(position, 240e3, 3) == 0.0);
  CHECK(world1.composition(position, 240e3, 4) == 0.0);
  CHECK(world1.composition(position, 240e3, 5) == 0.0);
  CHECK(world1.composition(position, 260e3, 0) == 0.0);
  CHECK(world1.composition(position, 260e3, 1) == 0.0);
  CHECK(world1.composition(position, 260e3, 2) == 0.0);
  CHECK(world1.composition(position, 260e3, 3) == 0.0);
  CHECK(world1.composition(position, 260e3, 4) == 0.0);
  CHECK(world1.composition(position, 260e3, 5) == 0.0);

  position = {{250e3,1750e3,0}};
  CHECK(world1.temperature(position, 0, 10) == Approx(293.15));
  CHECK(world1.temperature(position, 240e3, 10) == Approx(1659.0985664065));
  CHECK(world1.temperature(position, 260e3, 10) == Approx(1720.8246597128));

  CHECK(world1.composition(position, 0, 0) == 0.0);
  CHECK(world1.composition(position, 0, 1) == 0.0);
  CHECK(world1.composition(position, 0, 2) == 0.0);
  CHECK(world1.composition(position, 0, 3) == 0.0);
  CHECK(world1.composition(position, 0, 4) == 1.0);
  CHECK(world1.composition(position, 0, 5) == 0.0);
  CHECK(world1.composition(position, 240e3, 0) == 0.0);
  CHECK(world1.composition(position, 240e3, 1) == 0.0);
  CHECK(world1.composition(position, 240e3, 2) == 0.0);
  CHECK(world1.composition(position, 240e3, 3) == 0.0);
  CHECK(world1.composition(position, 240e3, 4) == 1.0);
  CHECK(world1.composition(position, 240e3, 5) == 0.0);
  CHECK(world1.composition(position, 260e3, 0) == 0.0);
  CHECK(world1.composition(position, 260e3, 1) == 0.0);
  CHECK(world1.composition(position, 260e3, 2) == 0.0);
  CHECK(world1.composition(position, 260e3, 3) == 0.0);
  CHECK(world1.composition(position, 260e3, 4) == 0.0);
  CHECK(world1.composition(position, 260e3, 5) == 0.0);

  position = {{750e3,250e3,0}};
  CHECK(world1.temperature(position, 0, 10) == Approx(10));
  CHECK(world1.temperature(position, 240e3, 10) == Approx(48.4));
  CHECK(world1.temperature(position, 260e3, 10) == Approx(1720.8246597128));

  CHECK(world1.composition(position, 0, 0) == 0.0);
  CHECK(world1.composition(position, 0, 1) == 0.0);
  CHECK(world1.composition(position, 0, 2) == 0.0);
  CHECK(world1.composition(position, 0, 3) == 0.0);
  CHECK(world1.composition(position, 0, 4) == 0.0);
  CHECK(world1.composition(position, 0, 5) == 0.25);
  CHECK(world1.composition(position, 0, 6) == 0.75);
  CHECK(world1.composition(position, 240e3, 0) == 0.0);
  CHECK(world1.composition(position, 240e3, 1) == 0.0);
  CHECK(world1.composition(position, 240e3, 2) == 0.0);
  CHECK(world1.composition(position, 240e3, 3) == 0.0);
  CHECK(world1.composition(position, 240e3, 4) == 0.0);
  CHECK(world1.composition(position, 240e3, 5) == 0.25);
  CHECK(world1.composition(position, 240e3, 6) == 0.75);
  CHECK(world1.composition(position, 260e3, 0) == 0.0);
  CHECK(world1.composition(position, 260e3, 1) == 0.0);
  CHECK(world1.composition(position, 260e3, 2) == 0.0);
  CHECK(world1.composition(position, 260e3, 3) == 0.0);
  CHECK(world1.composition(position, 260e3, 4) == 0.0);
  CHECK(world1.composition(position, 260e3, 5) == 0.0);

  // the constant layers test
  position = {{1500e3,250e3,0}};
  CHECK(world1.temperature(position, 0, 10) == Approx(10));
  CHECK(world1.temperature(position, 240e3, 10) == Approx(48.4));
  CHECK(world1.temperature(position, 260e3, 10) == Approx(1720.8246597128));

  CHECK(world1.composition(position, 0, 0) == 0.0);
  CHECK(world1.composition(position, 0, 1) == 0.0);
  CHECK(world1.composition(position, 0, 2) == 0.0);
  CHECK(world1.composition(position, 0, 3) == 0.0);
  CHECK(world1.composition(position, 0, 4) == 0.0);
  CHECK(world1.composition(position, 0, 5) == 0.0);
  CHECK(world1.composition(position, 0, 6) == 0.25);
  CHECK(world1.composition(position, 0, 7) == 0.75);
  CHECK(world1.composition(position, 0, 8) == 0.0);
  CHECK(world1.composition(position, 75e3-1, 0) == 0.0);
  CHECK(world1.composition(position, 75e3-1, 1) == 0.0);
  CHECK(world1.composition(position, 75e3-1, 2) == 0.0);
  CHECK(world1.composition(position, 75e3-1, 3) == 0.0);
  CHECK(world1.composition(position, 75e3-1, 4) == 0.0);
  CHECK(world1.composition(position, 75e3-1, 5) == 0.0);
  CHECK(world1.composition(position, 75e3-1, 6) == 0.25);
  CHECK(world1.composition(position, 75e3-1, 7) == 0.75);
  CHECK(world1.composition(position, 75e3-1, 8) == 0.0);
  CHECK(world1.composition(position, 75e3+1, 0) == 0.0);
  CHECK(world1.composition(position, 75e3+1, 1) == 0.0);
  CHECK(world1.composition(position, 75e3+1, 2) == 0.0);
  CHECK(world1.composition(position, 75e3+1, 3) == 0.0);
  CHECK(world1.composition(position, 75e3+1, 4) == 0.0);
  CHECK(world1.composition(position, 75e3+1, 5) == 0.0);
  CHECK(world1.composition(position, 75e3+1, 6) == 0.0);
  CHECK(world1.composition(position, 75e3+1, 7) == 1.0);
  CHECK(world1.composition(position, 75e3+1, 8) == 0.0);
  CHECK(world1.composition(position, 150e3-1, 0) == 0.0);
  CHECK(world1.composition(position, 150e3-1, 1) == 0.0);
  CHECK(world1.composition(position, 150e3-1, 2) == 0.0);
  CHECK(world1.composition(position, 150e3-1, 3) == 0.0);
  CHECK(world1.composition(position, 150e3-1, 4) == 0.0);
  CHECK(world1.composition(position, 150e3-1, 5) == 0.0);
  CHECK(world1.composition(position, 150e3-1, 6) == 0.0);
  CHECK(world1.composition(position, 150e3-1, 7) == 1.0);
  CHECK(world1.composition(position, 150e3-1, 8) == 0.0);
  CHECK(world1.composition(position, 150e3+1, 0) == 0.0);
  CHECK(world1.composition(position, 150e3+1, 1) == 0.0);
  CHECK(world1.composition(position, 150e3+1, 2) == 0.0);
  CHECK(world1.composition(position, 150e3+1, 3) == 0.0);
  CHECK(world1.composition(position, 150e3+1, 4) == 0.0);
  CHECK(world1.composition(position, 150e3+1, 5) == 0.0);
  CHECK(world1.composition(position, 150e3+1, 6) == 0.0);
  CHECK(world1.composition(position, 150e3+1, 7) == 0.0);
  CHECK(world1.composition(position, 150e3+1, 8) == 1.0);
  CHECK(world1.composition(position, 240e3, 0) == 0.0);
  CHECK(world1.composition(position, 240e3, 1) == 0.0);
  CHECK(world1.composition(position, 240e3, 2) == 0.0);
  CHECK(world1.composition(position, 240e3, 3) == 0.0);
  CHECK(world1.composition(position, 240e3, 4) == 0.0);
  CHECK(world1.composition(position, 240e3, 5) == 0.0);
  CHECK(world1.composition(position, 240e3, 6) == 0.0);
  CHECK(world1.composition(position, 240e3, 7) == 0.0);
  CHECK(world1.composition(position, 240e3, 8) == 0.0);
  CHECK(world1.composition(position, 260e3, 0) == 0.0);
  CHECK(world1.composition(position, 260e3, 1) == 0.0);
  CHECK(world1.composition(position, 260e3, 2) == 0.0);
  CHECK(world1.composition(position, 260e3, 3) == 0.0);
  CHECK(world1.composition(position, 260e3, 4) == 0.0);
  CHECK(world1.composition(position, 260e3, 5) == 0.0);
  CHECK(world1.composition(position, 260e3, 6) == 0.0);
  CHECK(world1.composition(position, 260e3, 7) == 0.0);
  CHECK(world1.composition(position, 260e3, 8) == 0.0);
}

TEST_CASE("WorldBuilder Features: Mantle layer")
{
  std::string file_name = WorldBuilder::Data::WORLD_BUILDER_SOURCE_DIR + "/tests/data/mantle_layer_cartesian.wb";
  WorldBuilder::World world1(file_name);

  // Check continental plate directly
  std::unique_ptr<Features::Interface> mantle_layer = Features::Interface::create("mantle layer", &world1);

  // Check continental plate through the world
  std::array<double,3> position = {{0,0,0}};
  CHECK(world1.temperature(position, 0, 10) == Approx(1600));

  position = {{250e3,501e3,0}};
  CHECK(world1.temperature(position, 0+100e3, 10) == Approx(150));
  CHECK(world1.temperature(position, 240e3+100e3, 10) == Approx(150));
  CHECK(world1.temperature(position, 260e3+100e3, 10) == Approx(1769.6886536946));

  CHECK(world1.composition(position, 0+200e3, 0) == 0.0);
  CHECK(world1.composition(position, 0+200e3, 1) == 0.0);
  CHECK(world1.composition(position, 0+200e3, 2) == 0.0);
  CHECK(world1.composition(position, 0+200e3, 3) == 1.0);
  CHECK(world1.composition(position, 0+200e3, 4) == 0.0);
  CHECK(world1.composition(position, 0+200e3, 5) == 0.0);
  CHECK(world1.composition(position, 240e3+200e3, 0) == 0.0);
  CHECK(world1.composition(position, 240e3+200e3, 1) == 0.0);
  CHECK(world1.composition(position, 240e3+200e3, 2) == 0.0);
  CHECK(world1.composition(position, 240e3+200e3, 3) == 1.0);
  CHECK(world1.composition(position, 240e3+200e3, 4) == 0.0);
  CHECK(world1.composition(position, 240e3+200e3, 5) == 0.0);
  CHECK(world1.composition(position, 260e3+200e3, 0) == 0.0);
  CHECK(world1.composition(position, 260e3+200e3, 1) == 0.0);
  CHECK(world1.composition(position, 260e3+200e3, 2) == 0.0);
  CHECK(world1.composition(position, 260e3+200e3, 3) == 0.0);
  CHECK(world1.composition(position, 260e3+200e3, 4) == 0.0);
  CHECK(world1.composition(position, 260e3+200e3, 5) == 0.0);

  position = {{1500e3,1500e3,0}};
  CHECK(world1.temperature(position, 0+150e3, 10) == Approx(20));
  CHECK(world1.temperature(position, 240e3+150e3, 10) == Approx(20));
  CHECK(world1.temperature(position, 260e3+150e3, 10) == Approx(1794.6385365126));

  CHECK(world1.composition(position, 0+150e3, 0) == 0.0);
  CHECK(world1.composition(position, 0+150e3, 1) == 0.0);
  CHECK(world1.composition(position, 0+150e3, 2) == 1.0);
  CHECK(world1.composition(position, 0+150e3, 3) == 0.0);
  CHECK(world1.composition(position, 0+150e3, 4) == 0.0);
  CHECK(world1.composition(position, 0+150e3, 5) == 0.0);
  CHECK(world1.composition(position, 240e3+150e3, 0) == 0.0);
  CHECK(world1.composition(position, 240e3+150e3, 1) == 0.0);
  CHECK(world1.composition(position, 240e3+150e3, 2) == 1.0);
  CHECK(world1.composition(position, 240e3+150e3, 3) == 0.0);
  CHECK(world1.composition(position, 240e3+150e3, 4) == 0.0);
  CHECK(world1.composition(position, 240e3+150e3, 5) == 0.0);
  CHECK(world1.composition(position, 260e3+150e3, 0) == 0.0);
  CHECK(world1.composition(position, 260e3+150e3, 1) == 0.0);
  CHECK(world1.composition(position, 260e3+150e3, 2) == 0.0);
  CHECK(world1.composition(position, 260e3+150e3, 3) == 0.0);
  CHECK(world1.composition(position, 260e3+150e3, 4) == 0.0);
  CHECK(world1.composition(position, 260e3+150e3, 5) == 0.0);

  position = {{250e3,1750e3,0}};
  CHECK(world1.temperature(position, 0+250e3, 10) == Approx(293.15));
  CHECK(world1.temperature(position, 240e3+250e3, 10) == Approx(1778.5465550447));
  CHECK(world1.temperature(position, 260e3+250e3, 10) == Approx(1845.598526046));

  CHECK(world1.composition(position, 0+250e3, 0) == 0.0);
  CHECK(world1.composition(position, 0+250e3, 1) == 0.0);
  CHECK(world1.composition(position, 0+250e3, 2) == 0.0);
  CHECK(world1.composition(position, 0+250e3, 3) == 0.0);
  CHECK(world1.composition(position, 0+250e3, 4) == 1.0);
  CHECK(world1.composition(position, 0+250e3, 5) == 0.0);
  CHECK(world1.composition(position, 240e3+250e3, 0) == 0.0);
  CHECK(world1.composition(position, 240e3+250e3, 1) == 0.0);
  CHECK(world1.composition(position, 240e3+250e3, 2) == 0.0);
  CHECK(world1.composition(position, 240e3+250e3, 3) == 0.0);
  CHECK(world1.composition(position, 240e3+250e3, 4) == 1.0);
  CHECK(world1.composition(position, 240e3+250e3, 5) == 0.0);
  CHECK(world1.composition(position, 260e3+250e3, 0) == 0.0);
  CHECK(world1.composition(position, 260e3+250e3, 1) == 0.0);
  CHECK(world1.composition(position, 260e3+250e3, 2) == 0.0);
  CHECK(world1.composition(position, 260e3+250e3, 3) == 0.0);
  CHECK(world1.composition(position, 260e3+250e3, 4) == 0.0);
  CHECK(world1.composition(position, 260e3+250e3, 5) == 0.0);

  position = {{750e3,250e3,0}};
  CHECK(world1.temperature(position, 0+300e3, 10) == Approx(10));
  CHECK(world1.temperature(position, 240e3+300e3, 10) == Approx(48.4));
  CHECK(world1.temperature(position, 260e3+300e3, 10) == Approx(1871.6186210824));

  CHECK(world1.composition(position, 0+300e3, 0) == 0.0);
  CHECK(world1.composition(position, 0+300e3, 1) == 0.0);
  CHECK(world1.composition(position, 0+300e3, 2) == 0.0);
  CHECK(world1.composition(position, 0+300e3, 3) == 0.0);
  CHECK(world1.composition(position, 0+300e3, 4) == 0.0);
  CHECK(world1.composition(position, 0+300e3, 5) == 0.25);
  CHECK(world1.composition(position, 0+300e3, 6) == 0.75);
  CHECK(world1.composition(position, 240e3+300e3, 0) == 0.0);
  CHECK(world1.composition(position, 240e3+300e3, 1) == 0.0);
  CHECK(world1.composition(position, 240e3+300e3, 2) == 0.0);
  CHECK(world1.composition(position, 240e3+300e3, 3) == 0.0);
  CHECK(world1.composition(position, 240e3+300e3, 4) == 0.0);
  CHECK(world1.composition(position, 240e3+300e3, 5) == 0.25);
  CHECK(world1.composition(position, 240e3+300e3, 6) == 0.75);
  CHECK(world1.composition(position, 260e3+300e3, 0) == 0.0);
  CHECK(world1.composition(position, 260e3+300e3, 1) == 0.0);
  CHECK(world1.composition(position, 260e3+300e3, 2) == 0.0);
  CHECK(world1.composition(position, 260e3+300e3, 3) == 0.0);
  CHECK(world1.composition(position, 260e3+300e3, 4) == 0.0);
  CHECK(world1.composition(position, 260e3+300e3, 5) == 0.0);

  // the constant layers test
  position = {{1500e3,250e3,0}};
  CHECK(world1.temperature(position, 0+350e3, 10) == Approx(1764.7404561736));
  CHECK(world1.temperature(position, 240e3+350e3, 10) == Approx(1887.4064334793));
  CHECK(world1.temperature(position, 260e3+350e3, 10) == Approx(1898.0055593602));

  CHECK(world1.composition(position, 0+350e3, 0) == 0.0);
  CHECK(world1.composition(position, 0+350e3, 1) == 0.0);
  CHECK(world1.composition(position, 0+350e3, 2) == 0.0);
  CHECK(world1.composition(position, 0+350e3, 3) == 0.0);
  CHECK(world1.composition(position, 0+350e3, 4) == 0.0);
  CHECK(world1.composition(position, 0+350e3, 5) == 0.0);
  CHECK(world1.composition(position, 0+350e3, 6) == 1.0);
  CHECK(world1.composition(position, 0+350e3, 7) == 0.0);
  CHECK(world1.composition(position, 0+350e3, 8) == 0.0);
  CHECK(world1.composition(position, 0+350e3, 9) == 0.0);
  CHECK(world1.composition(position, 75e3-1+350e3, 0) == 0.0);
  CHECK(world1.composition(position, 75e3-1+350e3, 1) == 0.0);
  CHECK(world1.composition(position, 75e3-1+350e3, 2) == 0.0);
  CHECK(world1.composition(position, 75e3-1+350e3, 3) == 0.0);
  CHECK(world1.composition(position, 75e3-1+350e3, 4) == 0.0);
  CHECK(world1.composition(position, 75e3-1+350e3, 5) == 0.0);
  CHECK(world1.composition(position, 75e3-1+350e3, 6) == 1.0);
  CHECK(world1.composition(position, 75e3-1+350e3, 7) == 0.0);
  CHECK(world1.composition(position, 75e3-1+350e3, 8) == 0.0);
  CHECK(world1.composition(position, 75e3-1+350e3, 9) == 0.0);
  CHECK(world1.composition(position, 75e3+1+350e3, 0) == 0.0);
  CHECK(world1.composition(position, 75e3+1+350e3, 1) == 0.0);
  CHECK(world1.composition(position, 75e3+1+350e3, 2) == 0.0);
  CHECK(world1.composition(position, 75e3+1+350e3, 3) == 0.0);
  CHECK(world1.composition(position, 75e3+1+350e3, 4) == 0.0);
  CHECK(world1.composition(position, 75e3+1+350e3, 5) == 0.0);
  CHECK(world1.composition(position, 75e3+1+350e3, 6) == 0.0);
  CHECK(world1.composition(position, 75e3+1+350e3, 7) == 1.0);
  CHECK(world1.composition(position, 75e3+1+350e3, 8) == 0.0);
  CHECK(world1.composition(position, 75e3+1+350e3, 9) == 0.0);
  CHECK(world1.composition(position, 150e3-1+350e3, 0) == 0.0);
  CHECK(world1.composition(position, 150e3-1+350e3, 1) == 0.0);
  CHECK(world1.composition(position, 150e3-1+350e3, 2) == 0.0);
  CHECK(world1.composition(position, 150e3-1+350e3, 3) == 0.0);
  CHECK(world1.composition(position, 150e3-1+350e3, 4) == 0.0);
  CHECK(world1.composition(position, 150e3-1+350e3, 5) == 0.0);
  CHECK(world1.composition(position, 150e3-1+350e3, 6) == 0.0);
  CHECK(world1.composition(position, 150e3-1+350e3, 7) == 1.0);
  CHECK(world1.composition(position, 150e3-1+350e3, 8) == 0.0);
  CHECK(world1.composition(position, 150e3-1+350e3, 9) == 0.0);
  CHECK(world1.composition(position, 150e3+1+350e3, 0) == 0.0);
  CHECK(world1.composition(position, 150e3+1+350e3, 1) == 0.0);
  CHECK(world1.composition(position, 150e3+1+350e3, 2) == 0.0);
  CHECK(world1.composition(position, 150e3+1+350e3, 3) == 0.0);
  CHECK(world1.composition(position, 150e3+1+350e3, 4) == 0.0);
  CHECK(world1.composition(position, 150e3+1+350e3, 5) == 0.0);
  CHECK(world1.composition(position, 150e3+1+350e3, 6) == 0.0);
  CHECK(world1.composition(position, 150e3+1+350e3, 7) == 0.0);
  CHECK(world1.composition(position, 150e3+1+350e3, 8) == 0.25);
  CHECK(world1.composition(position, 150e3+1+350e3, 9) == 0.75);
  CHECK(world1.composition(position, 240e3+350e3, 0) == 0.0);
  CHECK(world1.composition(position, 240e3+350e3, 1) == 0.0);
  CHECK(world1.composition(position, 240e3+350e3, 2) == 0.0);
  CHECK(world1.composition(position, 240e3+350e3, 3) == 0.0);
  CHECK(world1.composition(position, 240e3+350e3, 4) == 0.0);
  CHECK(world1.composition(position, 240e3+350e3, 5) == 0.0);
  CHECK(world1.composition(position, 240e3+350e3, 6) == 0.0);
  CHECK(world1.composition(position, 240e3+350e3, 7) == 0.0);
  CHECK(world1.composition(position, 240e3+350e3, 8) == 0.0);
  CHECK(world1.composition(position, 240e3+350e3, 9) == 0.0);
  CHECK(world1.composition(position, 260e3+350e3, 0) == 0.0);
  CHECK(world1.composition(position, 260e3+350e3, 1) == 0.0);
  CHECK(world1.composition(position, 260e3+350e3, 2) == 0.0);
  CHECK(world1.composition(position, 260e3+350e3, 3) == 0.0);
  CHECK(world1.composition(position, 260e3+350e3, 4) == 0.0);
  CHECK(world1.composition(position, 260e3+350e3, 5) == 0.0);
  CHECK(world1.composition(position, 260e3+350e3, 6) == 0.0);
  CHECK(world1.composition(position, 260e3+350e3, 7) == 0.0);
  CHECK(world1.composition(position, 260e3+350e3, 8) == 0.0);
  CHECK(world1.composition(position, 260e3+350e3, 9) == 0.0);
}

TEST_CASE("WorldBuilder Features: Oceanic Plate")
{
  // Cartesian
  std::string file_name = WorldBuilder::Data::WORLD_BUILDER_SOURCE_DIR + "/tests/data/oceanic_plate_cartesian.wb";
  WorldBuilder::World world1(file_name);

  // Check continental plate directly
  std::unique_ptr<Features::Interface> continental_plate = Features::Interface::create("oceanic plate", &world1);

  // Check continental plate through the world
  // 2d
  std::array<double,2> position_2d = {{0,0}};
  CHECK(world1.temperature(position_2d, 0, 10) == Approx(1600));
  CHECK(world1.temperature(position_2d, 240e3, 10) == Approx(1711.2149738521));
  CHECK(world1.temperature(position_2d, 260e3, 10) == Approx(1720.8246597128));
  CHECK(world1.composition(position_2d, 0, 0) == 0.0);
  CHECK(world1.composition(position_2d, 0, 1) == 0.0);
  CHECK(world1.composition(position_2d, 0, 2) == 0.0);
  CHECK(world1.composition(position_2d, 0, 3) == 0.0);
  CHECK(world1.composition(position_2d, 0, 4) == 0.0);
  CHECK(world1.composition(position_2d, 0, 5) == 0.0);
  CHECK(world1.composition(position_2d, 0, 6) == 0.0);
  // 3d
  std::array<double,3> position = {{0,0,0}};
  CHECK(world1.temperature(position, 0, 10) == Approx(1600));
  CHECK(world1.temperature(position, 240e3, 10) == Approx(1711.2149738521));
  CHECK(world1.temperature(position, 260e3, 10) == Approx(1720.8246597128));
  CHECK(world1.composition(position, 0, 0) == 0.0);
  CHECK(world1.composition(position, 0, 1) == 0.0);
  CHECK(world1.composition(position, 0, 2) == 0.0);
  CHECK(world1.composition(position, 0, 3) == 0.0);
  CHECK(world1.composition(position, 0, 4) == 0.0);
  CHECK(world1.composition(position, 0, 5) == 0.0);
  CHECK(world1.composition(position, 0, 6) == 0.0);

  position = {{250e3,500e3,0}};
  CHECK(world1.temperature(position, 0, 10) == Approx(150));
  CHECK(world1.temperature(position, 240e3, 10) == Approx(150));
  CHECK(world1.temperature(position, 260e3, 10) == Approx(1720.8246597128));
  CHECK(world1.composition(position, 0, 0) == 0.0);
  CHECK(world1.composition(position, 0, 1) == 0.0);
  CHECK(world1.composition(position, 0, 2) == 0.0);
  CHECK(world1.composition(position, 0, 3) == 1.0);
  CHECK(world1.composition(position, 240e3, 3) == 1.0);
  CHECK(world1.composition(position, 260e3, 3) == 0.0);
  CHECK(world1.composition(position, 0, 4) == 0.0);
  CHECK(world1.composition(position, 0, 5) == 0.0);
  CHECK(world1.composition(position, 0, 6) == 0.0);

  position = {{1500e3,1500e3,0}};
  CHECK(world1.temperature(position, 0, 10) == Approx(20));
  CHECK(world1.temperature(position, 240e3, 10) == Approx(20));
  CHECK(world1.temperature(position, 260e3, 10) == Approx(1720.8246597128));
  CHECK(world1.composition(position, 0, 0) == 0.0);
  CHECK(world1.composition(position, 0, 1) == 0.0);
  CHECK(world1.composition(position, 0, 2) == 1.0);
  CHECK(world1.composition(position, 240e3, 2) == 1.0);
  CHECK(world1.composition(position, 260e3, 2) == 0.0);
  CHECK(world1.composition(position, 0, 3) == 0.0);
  CHECK(world1.composition(position, 0, 4) == 0.0);
  CHECK(world1.composition(position, 0, 5) == 0.0);
  CHECK(world1.composition(position, 0, 6) == 0.0);

  position = {{250e3,1750e3,0}};
  CHECK(world1.temperature(position, 0, 10) == Approx(293.15));
  CHECK(world1.temperature(position, 240e3, 10) == Approx(1659.0985664065));
  CHECK(world1.temperature(position, 260e3, 10) == Approx(1720.8246597128));
  CHECK(world1.composition(position, 0, 0) == 0.0);
  CHECK(world1.composition(position, 0, 1) == 0.0);
  CHECK(world1.composition(position, 0, 2) == 0.0);
  CHECK(world1.composition(position, 0, 3) == 0.0);
  CHECK(world1.composition(position, 0, 4) == 1.0);
  CHECK(world1.composition(position, 240e3, 4) == 1.0);
  CHECK(world1.composition(position, 260e3, 4) == 0.0);
  CHECK(world1.composition(position, 0, 5) == 0.0);
  CHECK(world1.composition(position, 0, 6) == 0.0);

  position = {{750e3,250e3,0}};
  CHECK(world1.temperature(position, 0, 10) == Approx(10));
  CHECK(world1.temperature(position, 240e3, 10) == Approx(48.4));
  CHECK(world1.temperature(position, 260e3, 10) == Approx(1720.8246597128));
  CHECK(world1.composition(position, 0, 0) == 0.0);
  CHECK(world1.composition(position, 0, 1) == 0.0);
  CHECK(world1.composition(position, 0, 2) == 0.0);
  CHECK(world1.composition(position, 0, 3) == 0.0);
  CHECK(world1.composition(position, 0, 4) == 0.0);
  CHECK(world1.composition(position, 0, 5) == 0.25);
  CHECK(world1.composition(position, 0, 6) == 0.75);
  CHECK(world1.composition(position, 240e3, 5) == 0.25);
  CHECK(world1.composition(position, 240e3, 6) == 0.75);
  CHECK(world1.composition(position, 260e3, 5) == 0.0);
  CHECK(world1.composition(position, 260e3, 6) == 0.0);

  position = {{1500e3, 0, 0}};
  CHECK(world1.temperature(position, 0, 10) == Approx(293.15));
  CHECK(world1.temperature(position, 10, 10) == Approx(303.6570169192));
  CHECK(world1.temperature(position, 240e3, 10) == Approx(1710.9310013));
  CHECK(world1.temperature(position, 260e3, 10) == Approx(1720.8246597128));
  CHECK(world1.composition(position, 0, 0) == 0.0);
  CHECK(world1.composition(position, 0, 1) == 0.0);
  CHECK(world1.composition(position, 0, 2) == 0.0);
  CHECK(world1.composition(position, 0, 3) == 0.0);
  CHECK(world1.composition(position, 0, 4) == 0.0);
  CHECK(world1.composition(position, 0, 5) == 0.0);
  CHECK(world1.composition(position, 0, 6) == 0.0);

  // test symmetry
  position = {{1600e3, 0, 0}};
  CHECK(world1.temperature(position, 0, 10) == Approx(293.15));
  CHECK(world1.temperature(position, 10, 10) == Approx(293.6215585565));
  CHECK(world1.temperature(position, 240e3, 10) == Approx(1711.2149738521));
  CHECK(world1.temperature(position, 260e3, 10) == Approx(1720.8246597128));

  position = {{1400e3, 0, 0}};
  CHECK(world1.temperature(position, 0, 10) == Approx(293.15));
  CHECK(world1.temperature(position, 10, 10) == Approx(293.6215585565));
  CHECK(world1.temperature(position, 240e3, 10) == Approx(1711.2149738521));
  CHECK(world1.temperature(position, 260e3, 10) == Approx(1720.8246597128));

  // the constant layers test
  position = {{200e3,200e3,0}};
  CHECK(world1.temperature(position, 0, 10) == Approx(293.15));
  CHECK(world1.temperature(position, 240e3, 10) == Approx(1708.6787610897));
  CHECK(world1.temperature(position, 260e3, 10) == Approx(1720.8246597128));

  CHECK(world1.composition(position, 0, 0) == 0.0);
  CHECK(world1.composition(position, 0, 1) == 0.0);
  CHECK(world1.composition(position, 0, 2) == 0.0);
  CHECK(world1.composition(position, 0, 3) == 0.0);
  CHECK(world1.composition(position, 0, 4) == 0.0);
  CHECK(world1.composition(position, 0, 5) == 0.0);
  CHECK(world1.composition(position, 0, 6) == 1.0);
  CHECK(world1.composition(position, 0, 7) == 0.0);
  CHECK(world1.composition(position, 0, 8) == 0.0);
  CHECK(world1.composition(position, 0, 9) == 0.0);
  CHECK(world1.composition(position, 75e3-1, 0) == 0.0);
  CHECK(world1.composition(position, 75e3-1, 1) == 0.0);
  CHECK(world1.composition(position, 75e3-1, 2) == 0.0);
  CHECK(world1.composition(position, 75e3-1, 3) == 0.0);
  CHECK(world1.composition(position, 75e3-1, 4) == 0.0);
  CHECK(world1.composition(position, 75e3-1, 5) == 0.0);
  CHECK(world1.composition(position, 75e3-1, 6) == 1.0);
  CHECK(world1.composition(position, 75e3-1, 7) == 0.0);
  CHECK(world1.composition(position, 75e3-1, 8) == 0.0);
  CHECK(world1.composition(position, 75e3-1, 9) == 0.0);
  CHECK(world1.composition(position, 75e3+1, 0) == 0.0);
  CHECK(world1.composition(position, 75e3+1, 1) == 0.0);
  CHECK(world1.composition(position, 75e3+1, 2) == 0.0);
  CHECK(world1.composition(position, 75e3+1, 3) == 0.0);
  CHECK(world1.composition(position, 75e3+1, 4) == 0.0);
  CHECK(world1.composition(position, 75e3+1, 5) == 0.0);
  CHECK(world1.composition(position, 75e3+1, 6) == 0.0);
  CHECK(world1.composition(position, 75e3+1, 7) == 1.0);
  CHECK(world1.composition(position, 75e3+1, 8) == 0.0);
  CHECK(world1.composition(position, 75e3+1, 9) == 0.0);
  CHECK(world1.composition(position, 150e3-1, 0) == 0.0);
  CHECK(world1.composition(position, 150e3-1, 1) == 0.0);
  CHECK(world1.composition(position, 150e3-1, 2) == 0.0);
  CHECK(world1.composition(position, 150e3-1, 3) == 0.0);
  CHECK(world1.composition(position, 150e3-1, 4) == 0.0);
  CHECK(world1.composition(position, 150e3-1, 5) == 0.0);
  CHECK(world1.composition(position, 150e3-1, 6) == 0.0);
  CHECK(world1.composition(position, 150e3-1, 7) == 1.0);
  CHECK(world1.composition(position, 150e3-1, 8) == 0.0);
  CHECK(world1.composition(position, 150e3-1, 9) == 0.0);
  CHECK(world1.composition(position, 150e3+1, 0) == 0.0);
  CHECK(world1.composition(position, 150e3+1, 1) == 0.0);
  CHECK(world1.composition(position, 150e3+1, 2) == 0.0);
  CHECK(world1.composition(position, 150e3+1, 3) == 0.0);
  CHECK(world1.composition(position, 150e3+1, 4) == 0.0);
  CHECK(world1.composition(position, 150e3+1, 5) == 0.0);
  CHECK(world1.composition(position, 150e3+1, 6) == 0.0);
  CHECK(world1.composition(position, 150e3+1, 7) == 0.0);
  CHECK(world1.composition(position, 150e3+1, 8) == 0.25);
  CHECK(world1.composition(position, 150e3+1, 9) == 0.75);
  CHECK(world1.composition(position, 240e3, 0) == 0.0);
  CHECK(world1.composition(position, 240e3, 1) == 0.0);
  CHECK(world1.composition(position, 240e3, 2) == 0.0);
  CHECK(world1.composition(position, 240e3, 3) == 0.0);
  CHECK(world1.composition(position, 240e3, 4) == 0.0);
  CHECK(world1.composition(position, 240e3, 5) == 0.0);
  CHECK(world1.composition(position, 240e3, 6) == 0.0);
  CHECK(world1.composition(position, 240e3, 7) == 0.0);
  CHECK(world1.composition(position, 240e3, 8) == 0.0);
  CHECK(world1.composition(position, 240e3, 9) == 0.0);
  CHECK(world1.composition(position, 260e3, 0) == 0.0);
  CHECK(world1.composition(position, 260e3, 1) == 0.0);
  CHECK(world1.composition(position, 260e3, 2) == 0.0);
  CHECK(world1.composition(position, 260e3, 3) == 0.0);
  CHECK(world1.composition(position, 260e3, 4) == 0.0);
  CHECK(world1.composition(position, 260e3, 5) == 0.0);
  CHECK(world1.composition(position, 260e3, 6) == 0.0);
  CHECK(world1.composition(position, 260e3, 7) == 0.0);
  CHECK(world1.composition(position, 260e3, 8) == 0.0);
  CHECK(world1.composition(position, 260e3, 9) == 0.0);

  // spherical
  file_name = WorldBuilder::Data::WORLD_BUILDER_SOURCE_DIR + "/tests/data/oceanic_plate_spherical.wb";
  WorldBuilder::World world2(file_name);

  // Check continental plate directly
  std::unique_ptr<Features::Interface> oceanic_plate = Features::Interface::create("oceanic plate", &world2);

  // Check continental plate through the world
  double dtr = Utilities::const_pi / 180.0;
  std::unique_ptr<WorldBuilder::CoordinateSystems::Interface> &coordinate_system = world2.parameters.coordinate_system;

  // 2d
  position_2d = {{6371000,0}};
  CHECK(world2.temperature(position_2d, 0, 10) == Approx(1600));
  CHECK(world2.composition(position_2d, 0, 0) == 0.0);
  CHECK(world2.composition(position_2d, 0, 1) == 0.0);
  CHECK(world2.composition(position_2d, 0, 2) == 0.0);
  CHECK(world2.composition(position_2d, 0, 3) == 0.0);
  CHECK(world2.composition(position_2d, 0, 4) == 0.0);
  CHECK(world2.composition(position_2d, 0, 5) == 0.0);
  CHECK(world2.composition(position_2d, 0, 6) == 0.0);

  // 3d
  position = {{6371000,0,0}};
  CHECK(world2.temperature(position, 0, 10) == Approx(1600));
  CHECK(world2.composition(position, 0, 0) == 0.0);
  CHECK(world2.composition(position, 0, 1) == 0.0);
  CHECK(world2.composition(position, 0, 2) == 0.0);
  CHECK(world2.composition(position, 0, 3) == 0.0);
  CHECK(world2.composition(position, 0, 4) == 0.0);
  CHECK(world2.composition(position, 0, 5) == 0.0);
  CHECK(world2.composition(position, 0, 6) == 0.0);

  position = {{6371000, -5 * dtr,-5 * dtr}};
  position = coordinate_system->natural_to_cartesian_coordinates(position);
  CHECK(world2.temperature(position, 0, 10) == Approx(150));
  CHECK(world2.temperature(position, 240e3, 10) == Approx(150));
  CHECK(world2.temperature(position, 260e3, 10) == Approx(1720.8246597128));
  CHECK(world2.composition(position, 0, 0) == 0.0);
  CHECK(world2.composition(position, 0, 1) == 0.0);
  CHECK(world2.composition(position, 0, 2) == 0.0);
  CHECK(world2.composition(position, 0, 3) == 0.0);
  CHECK(world2.composition(position, 0, 4) == 0.0);
  CHECK(world2.composition(position, 0, 5) == 0.0);
  CHECK(world2.composition(position, 0, 6) == 0.0);

  position = {{6371000, 5 * dtr,-5 * dtr}};
  position = coordinate_system->natural_to_cartesian_coordinates(position);
  CHECK(world2.temperature(position, 0, 10) == Approx(20));
  CHECK(world2.temperature(position, 240e3, 10) == Approx(20));
  CHECK(world2.temperature(position, 260e3, 10) == Approx(1720.8246597128));
  CHECK(world2.composition(position, 0, 0) == 0.0);
  CHECK(world2.composition(position, 0, 1) == 0.0);
  CHECK(world2.composition(position, 0, 2) == 1.0);
  CHECK(world2.composition(position, 240e3, 2) == 1.0);
  CHECK(world2.composition(position, 260e3, 2) == 0.0);
  CHECK(world2.composition(position, 0, 3) == 0.0);
  CHECK(world2.composition(position, 0, 4) == 0.0);
  CHECK(world2.composition(position, 0, 5) == 0.0);
  CHECK(world2.composition(position, 0, 6) == 0.0);

  position = {{6371000, 5 * dtr,5 * dtr}};
  position = coordinate_system->natural_to_cartesian_coordinates(position);
  CHECK(world2.temperature(position, 0, 10) == Approx(293.15));
  CHECK(world2.temperature(position, 240e3, 10) == Approx(1659.0985664065));
  CHECK(world2.temperature(position, 260e3, 10) == Approx(1720.8246597128));
  CHECK(world2.composition(position, 0, 0) == 0.0);
  CHECK(world2.composition(position, 0, 1) == 0.0);
  CHECK(world2.composition(position, 0, 2) == 0.0);
  CHECK(world2.composition(position, 0, 3) == 0.0);
  CHECK(world2.composition(position, 0, 4) == 1.0);
  CHECK(world2.composition(position, 240e3, 4) == 1.0);
  CHECK(world2.composition(position, 260e3, 4) == 0.0);
  CHECK(world2.composition(position, 0, 5) == 0.0);
  CHECK(world2.composition(position, 0, 6) == 0.0);

  position = {{6371000, -15 * dtr, -15 * dtr}};
  position = coordinate_system->natural_to_cartesian_coordinates(position);
  CHECK(world2.temperature(position, 0, 10) == Approx(10));
  CHECK(world2.temperature(position, 240e3, 10) == Approx(48.4));
  CHECK(world2.temperature(position, 260e3, 10) == Approx(1720.8246597128));
  CHECK(world2.composition(position, 0, 0) == 0.0);
  CHECK(world2.composition(position, 0, 1) == 0.0);
  CHECK(world2.composition(position, 0, 2) == 0.0);
  CHECK(world2.composition(position, 0, 3) == 0.0);
  CHECK(world2.composition(position, 0, 4) == 0.0);
  CHECK(world2.composition(position, 0, 5) == 0.25);
  CHECK(world2.composition(position, 0, 6) == 0.75);
  CHECK(world2.composition(position, 240e3, 5) == 0.25);
  CHECK(world2.composition(position, 240e3, 6) == 0.75);
  CHECK(world2.composition(position, 260e3, 5) == 0.0);
  CHECK(world2.composition(position, 260e3, 6) == 0.0);

  position = {{6371000, 15 * dtr, -19 * dtr}};
  position = coordinate_system->natural_to_cartesian_coordinates(position);
  CHECK(world2.temperature(position, 0, 0) == Approx(293.15));
  CHECK(world2.temperature(position, 10, 10) == Approx(303.6570169192));
  CHECK(world2.temperature(position, 240e3, 10) == Approx(1710.9310013));
  CHECK(world2.temperature(position, 260e3, 10) == Approx(1720.8246597128));
  CHECK(world2.composition(position, 0, 0) == 0.0);
  CHECK(world2.composition(position, 0, 1) == 0.0);
  CHECK(world2.composition(position, 0, 2) == 0.0);
  CHECK(world2.composition(position, 0, 3) == 0.0);
  CHECK(world2.composition(position, 0, 4) == 0.0);
  CHECK(world2.composition(position, 0, 5) == 0.0);
  CHECK(world2.composition(position, 0, 6) == 1.0);
  CHECK(world2.composition(position, 240e3, 6) == 1.0);
  CHECK(world2.composition(position, 260e3, 6) == 0.0);

  // test symmetry
  position = {{6371000, 16 * dtr, -19 * dtr}};
  position = coordinate_system->natural_to_cartesian_coordinates(position);
  CHECK(world2.temperature(position, 0, 10) == Approx(293.15));
  CHECK(world2.temperature(position, 10, 10) == Approx(293.596373966));
  CHECK(world2.temperature(position, 240e3, 10) == Approx(1711.2149738521));
  CHECK(world2.temperature(position, 260e3, 10) == Approx(1720.8246597128));

  position = {{6371000, 14 * dtr, -19 * dtr}};
  position = coordinate_system->natural_to_cartesian_coordinates(position);
  CHECK(world2.temperature(position, 0, 10) == Approx(293.15));
  CHECK(world2.temperature(position, 10, 10) == Approx(293.596373966));
  CHECK(world2.temperature(position, 240e3, 10) == Approx(1711.2149738521));
  CHECK(world2.temperature(position, 260e3, 10) == Approx(1720.8246597128));

  // test bend
  position = {{6371000, 12.5 * dtr, -12.5 * dtr}};
  position = coordinate_system->natural_to_cartesian_coordinates(position);
  CHECK(world2.temperature(position, 0, 0) == Approx(293.15));
  CHECK(world2.temperature(position, 10, 10) == Approx(303.6570169192));
  CHECK(world2.temperature(position, 240e3, 10) == Approx(1710.9310013));
  CHECK(world2.temperature(position, 260e3, 10) == Approx(1720.8246597128));
}

TEST_CASE("WorldBuilder Features: Subducting Plate")
{
  // Cartesian
  std::string file_name = WorldBuilder::Data::WORLD_BUILDER_SOURCE_DIR + "/tests/data/subducting_plate_constant_angles_cartesian.wb";
  WorldBuilder::World world1(file_name);

  // Check continental plate directly (upper case should automatically turn into lower case).
  std::unique_ptr<Features::Interface> continental_plate = Features::Interface::create("Subducting Plate", &world1);

  // Check continental plate through the world
  std::array<double,3> position = {{0,0,800e3}};
  CHECK(world1.temperature(position, 0, 10) == Approx(1600.0));
  CHECK(world1.temperature(position, 240e3, 10) == Approx(1711.2149738521));
  CHECK(world1.temperature(position, 260e3, 10) == Approx(1720.8246597128));
  CHECK(world1.composition(position, 0, 0) == 0.0);
  CHECK(world1.composition(position, 0, 1) == 0.0);
  CHECK(world1.composition(position, 0, 2) == 0.0);
  CHECK(world1.composition(position, 0, 3) == 0.0);
  CHECK(world1.composition(position, 0, 4) == 0.0);
  CHECK(world1.composition(position, 0, 5) == 0.0);
  CHECK(world1.composition(position, 0, 6) == 0.0);

  position = {{250e3,500e3,800e3}};
  // results strongly dependent on the summation number of the McKenzie temperature.
  CHECK(world1.temperature(position, 0, 10) == Approx(1599.9999999994));
  CHECK(world1.temperature(position, 1, 10) == Approx(1590.6681048292)); // we are in the plate for sure (colder than anywhere in the mantle)
  CHECK(world1.temperature(position, 5, 10) == Approx(1554.3294579725)); // we are in the plate for sure (colder than anywhere in the mantle)
  CHECK(world1.temperature(position, 10, 10) == Approx(1511.0782994151)); // we are in the plate for sure (colder than anywhere in the mantle)
  CHECK(world1.temperature(position, 100, 10) == Approx(1050.2588653774)); // we are in the plate for sure (colder than anywhere in the mantle)
  CHECK(world1.temperature(position, 500, 10) == Approx(876.8637089978)); // we are in the plate for sure (colder than anywhere in the mantle)
  CHECK(world1.temperature(position, 1000, 10) == Approx(829.7197877204)); // we are in the plate for sure (colder than anywhere in the mantle)
  CHECK(world1.temperature(position, 5000, 10) == Approx(620.6827823798)); // we are in the plate for sure (colder than anywhere in the mantle)
  CHECK(world1.temperature(position, 10e3, 10) == Approx(525.727843761));
  CHECK(world1.temperature(position, 25e3, 10) == Approx(538.4605698191));
  CHECK(world1.temperature(position, 50e3, 10) == Approx(751.6318639633));
  CHECK(world1.temperature(position, 75e3, 10) == Approx(991.5852995579));
  CHECK(world1.temperature(position, 150e3, 10) == Approx(1668.6311660012));
  //CHECK(world1.temperature(position, std::sqrt(2) * 100e3 - 1, 10) == Approx(150.0));
  //CHECK(world1.temperature(position, std::sqrt(2) * 100e3 + 1, 10) == Approx(1664.6283561404));
  CHECK(world1.composition(position, 0, 0) == 1.0);
  CHECK(world1.composition(position, 0, 1) == 0.0);
  CHECK(world1.composition(position, 0, 2) == 0.0);
  CHECK(world1.composition(position, 0, 3) == 0.0);
  CHECK(world1.composition(position, 0, 4) == 0.0);
  CHECK(world1.composition(position, 10, 0) == 1.0);
  CHECK(world1.composition(position, 10, 1) == 0.0);
  CHECK(world1.composition(position, 10, 2) == 0.0);
  CHECK(world1.composition(position, 10, 3) == 0.0);
  CHECK(world1.composition(position, 10, 4) == 0.0);
  CHECK(world1.composition(position, std::sqrt(2) * 33e3 - 1, 0) == 1.0);
  CHECK(world1.composition(position, std::sqrt(2) * 33e3 + 1, 0) == 0.0);
  CHECK(world1.composition(position, std::sqrt(2) * 33e3 - 1, 1) == 0.0);
  CHECK(world1.composition(position, std::sqrt(2) * 33e3 + 1, 1) == 1.0);
  CHECK(world1.composition(position, std::sqrt(2) * 66e3 - 1, 1) == 1.0);
  CHECK(world1.composition(position, std::sqrt(2) * 66e3 + 1, 1) == 0.0);
  CHECK(world1.composition(position, std::sqrt(2) * 66e3 - 1, 2) == 0.0);
  CHECK(world1.composition(position, std::sqrt(2) * 66e3 + 1, 2) == 0.25);
  CHECK(world1.composition(position, std::sqrt(2) * 66e3 + 1, 3) == 0.75);
  CHECK(world1.composition(position, std::sqrt(2) * 99e3 - 1, 2) == 0.25);
  CHECK(world1.composition(position, std::sqrt(2) * 99e3 - 1, 3) == 0.75);
  CHECK(world1.composition(position, std::sqrt(2) * 99e3 + 1, 2) == 0.0);
  // this comes form the first subducting plate
  CHECK(world1.composition(position, std::sqrt(2) * 99e3 + 1, 3) == 1.0);
  CHECK(world1.composition(position, std::sqrt(2) * 100e3 - 1, 3) == 1.0);
  CHECK(world1.composition(position, std::sqrt(2) * 100e3 + 1, 3) == 0.0);
  CHECK(world1.composition(position, 0, 4) == 0.0);
  CHECK(world1.composition(position, 0, 5) == 0.0);
  CHECK(world1.composition(position, 0, 6) == 0.0);


  position = {{250e3,600e3,800e3}};
  CHECK(world1.temperature(position, 0, 10) == Approx(1600));
  CHECK(world1.temperature(position, 10, 10) == Approx(1600.0044800063));
  CHECK(world1.temperature(position, 100e3, 10) == Approx(1645.4330950743));
  CHECK(world1.temperature(position, 100e3+1, 10) == Approx(1.0000424264));
  CHECK(world1.temperature(position, 101e3, 10) == Approx(1.0424264069));
  CHECK(world1.temperature(position, 110e3, 10) == Approx(1.4242640687));
  CHECK(world1.temperature(position, 150e3, 10) == Approx(3.1213203436));
  CHECK(world1.temperature(position, 200e3, 10) == Approx(2.0));
  CHECK(world1.composition(position, 0, 0) == 0.0);
  CHECK(world1.composition(position, 0, 1) == 0.0);
  CHECK(world1.composition(position, 0, 2) == 0.0);
  CHECK(world1.composition(position, 0, 3) == 0.0);
  CHECK(world1.composition(position, 10, 0) == 0.0);
  CHECK(world1.composition(position, 10, 1) == 0.0);
  CHECK(world1.composition(position, 10, 2) == 0.0);
  CHECK(world1.composition(position, 10, 3) == 0.0);
  CHECK(world1.composition(position, 100e3, 0) == Approx(0.0));
  CHECK(world1.composition(position, 100e3, 1) == Approx(0.0));
  CHECK(world1.composition(position, 100e3, 2) == Approx(0.0));
  CHECK(world1.composition(position, 100e3, 3) == Approx(0.0));
  CHECK(world1.composition(position, 100e3, 4) == Approx(0.0));
  CHECK(world1.composition(position, 100e3+1, 0) == Approx(0.0));
  CHECK(world1.composition(position, 100e3+1, 1) == Approx(0.0));
  CHECK(world1.composition(position, 100e3+1, 2) == Approx(0.0));
  CHECK(world1.composition(position, 100e3+1, 3) == Approx(1.0));
  CHECK(world1.composition(position, 100e3+1, 4) == Approx(0.0));
  CHECK(world1.composition(position, 101e3, 0) == Approx(0.0));
  CHECK(world1.composition(position, 101e3+1, 1) == Approx(0.0));
  CHECK(world1.composition(position, 101e3+1, 2) == Approx(0.0));
  CHECK(world1.composition(position, 101e3+1, 3) == Approx(1.0));
  CHECK(world1.composition(position, 101e3+1, 4) == Approx(0.0));
  CHECK(world1.composition(position, 150e3, 0) == Approx(0.0));
  CHECK(world1.composition(position, 150e3, 1) == Approx(0.0));
  CHECK(world1.composition(position, 150e3, 2) == Approx(0.0));
  CHECK(world1.composition(position, 150e3, 3) == Approx(1.0));
  CHECK(world1.composition(position, 150e3, 4) == Approx(0.0));
  CHECK(world1.composition(position, 200e3, 0) == Approx(0.0));
  CHECK(world1.composition(position, 200e3, 1) == Approx(0.0));
  CHECK(world1.composition(position, 200e3, 2) == Approx(0.0));
  CHECK(world1.composition(position, 200e3, 3) == Approx(0.25));
  CHECK(world1.composition(position, 200e3, 4) == Approx(0.75));
  CHECK(world1.composition(position, 0, 4) == 0.0);
  CHECK(world1.composition(position, 0, 5) == 0.0);
  CHECK(world1.composition(position, 0, 6) == 0.0);

  position = {{650e3,650e3,800e3}};
  CHECK(world1.temperature(position, 0, 10) == Approx(1600));
  CHECK(world1.temperature(position, 10, 10) == Approx(1600.0044800063));
  CHECK(world1.temperature(position, 100e3, 10) == Approx(4.3590710784));
  CHECK(world1.temperature(position, 100e3+1, 10) == Approx(4.3590710784));
  CHECK(world1.temperature(position, 101e3, 10) == Approx(4.3590710784));
  CHECK(world1.temperature(position, 110e3, 10) == Approx(4.3590710784));
  CHECK(world1.temperature(position, 150e3, 10) == Approx(4.3590710784));
  CHECK(world1.temperature(position, 200e3, 10) == Approx(1692.1562939786));
  CHECK(world1.composition(position, 0, 0) == 0.0);
  CHECK(world1.composition(position, 0, 1) == 0.0);
  CHECK(world1.composition(position, 0, 2) == 0.0);
  CHECK(world1.composition(position, 0, 3) == 0.0);
  CHECK(world1.composition(position, 10, 0) == 0.0);
  CHECK(world1.composition(position, 10, 1) == 0.0);
  CHECK(world1.composition(position, 10, 2) == 0.0);
  CHECK(world1.composition(position, 10, 3) == 0.0);
  CHECK(world1.composition(position, 100e3, 0) == Approx(0.0));
  CHECK(world1.composition(position, 100e3, 1) == Approx(0.0));
  CHECK(world1.composition(position, 100e3, 2) == Approx(0.0));
  CHECK(world1.composition(position, 100e3, 3) == Approx(0.0897677696));
  CHECK(world1.composition(position, 100e3, 4) == Approx(0.2693033088));
  CHECK(world1.composition(position, 100e3+1, 0) == Approx(0.0));
  CHECK(world1.composition(position, 100e3+1, 1) == Approx(0.0));
  CHECK(world1.composition(position, 100e3+1, 2) == Approx(0.0));
  CHECK(world1.composition(position, 100e3+1, 3) == Approx(0.0897677696));
  CHECK(world1.composition(position, 100e3+1, 4) == Approx(0.2693033088));
  CHECK(world1.composition(position, 101e3, 0) == Approx(0.0));
  CHECK(world1.composition(position, 101e3+1, 1) == Approx(0.0));
  CHECK(world1.composition(position, 101e3+1, 2) == Approx(0.0));
  CHECK(world1.composition(position, 101e3+1, 3) == Approx(0.0897677696));
  CHECK(world1.composition(position, 101e3+1, 4) == Approx(0.2693033088));
  CHECK(world1.composition(position, 150e3, 0) == Approx(0.0));
  CHECK(world1.composition(position, 150e3, 1) == Approx(0.0));
  CHECK(world1.composition(position, 150e3, 2) == Approx(0.0));
  CHECK(world1.composition(position, 150e3, 3) == Approx(0.0897677696));
  CHECK(world1.composition(position, 150e3, 4) == Approx(0.2693033088));
  CHECK(world1.composition(position, 200e3, 0) == Approx(0.0));
  CHECK(world1.composition(position, 200e3, 1) == Approx(0.0));
  CHECK(world1.composition(position, 200e3, 2) == Approx(0.0));
  CHECK(world1.composition(position, 200e3, 3) == Approx(0.0));
  CHECK(world1.composition(position, 200e3, 4) == Approx(0.0));
  CHECK(world1.composition(position, 0, 4) == 0.0);
  CHECK(world1.composition(position, 0, 5) == 0.0);
  CHECK(world1.composition(position, 0, 6) == 0.0);

  position = {{700e3,675e3,800e3}};
  CHECK(world1.temperature(position, 0, 10) == Approx(1600));
  CHECK(world1.temperature(position, 10, 10) == Approx(1600.0044800063));
  CHECK(world1.temperature(position, 100e3, 10) == Approx(4.4592312705));
  CHECK(world1.temperature(position, 100e3+1, 10) == Approx(4.4592312705));
  CHECK(world1.temperature(position, 101e3, 10) == Approx(4.4592312705));
  CHECK(world1.temperature(position, 110e3, 10) == Approx(4.4592312705));
  CHECK(world1.temperature(position, 150e3, 10) == Approx(4.4592312705));
  CHECK(world1.temperature(position, 200e3, 10) == Approx(1692.1562939786));
  CHECK(world1.composition(position, 0, 0) == 0.0);
  CHECK(world1.composition(position, 0, 1) == 0.0);
  CHECK(world1.composition(position, 0, 2) == 0.0);
  CHECK(world1.composition(position, 0, 3) == 0.0);
  CHECK(world1.composition(position, 10, 0) == 0.0);
  CHECK(world1.composition(position, 10, 1) == 0.0);
  CHECK(world1.composition(position, 10, 2) == 0.0);
  CHECK(world1.composition(position, 10, 3) == 0.0);
  CHECK(world1.composition(position, 100e3, 0) == Approx(0.0));
  CHECK(world1.composition(position, 100e3, 1) == Approx(0.0));
  CHECK(world1.composition(position, 100e3, 2) == Approx(0.0));
  CHECK(world1.composition(position, 100e3, 3) == Approx(0.1148078176));
  CHECK(world1.composition(position, 100e3, 4) == Approx(0.3444234529));
  CHECK(world1.composition(position, 100e3+1, 0) == Approx(0.0));
  CHECK(world1.composition(position, 100e3+1, 1) == Approx(0.0));
  CHECK(world1.composition(position, 100e3+1, 2) == Approx(0.0));
  CHECK(world1.composition(position, 100e3+1, 3) == Approx(0.1148078176));
  CHECK(world1.composition(position, 100e3+1, 4) == Approx(0.3444234529));
  CHECK(world1.composition(position, 101e3, 0) == Approx(0.0));
  CHECK(world1.composition(position, 101e3+1, 1) == Approx(0.0));
  CHECK(world1.composition(position, 101e3+1, 2) == Approx(0.0));
  CHECK(world1.composition(position, 101e3+1, 3) == Approx(0.1148078176));
  CHECK(world1.composition(position, 101e3+1, 4) == Approx(0.3444234529));
  CHECK(world1.composition(position, 150e3, 0) == Approx(0.0));
  CHECK(world1.composition(position, 150e3, 1) == Approx(0.0));
  CHECK(world1.composition(position, 150e3, 2) == Approx(0.0));
  CHECK(world1.composition(position, 150e3, 3) == Approx(0.1148078176));
  CHECK(world1.composition(position, 150e3, 4) == Approx(0.3444234529));
  CHECK(world1.composition(position, 200e3, 0) == Approx(0.0));
  CHECK(world1.composition(position, 200e3, 1) == Approx(0.0));
  CHECK(world1.composition(position, 200e3, 2) == Approx(0.0));
  CHECK(world1.composition(position, 200e3, 3) == Approx(0.0));
  CHECK(world1.composition(position, 200e3, 4) == Approx(0.0));
  CHECK(world1.composition(position, 0, 4) == 0.0);
  CHECK(world1.composition(position, 0, 5) == 0.0);
  CHECK(world1.composition(position, 0, 6) == 0.0);

  std::string file_name2 = WorldBuilder::Data::WORLD_BUILDER_SOURCE_DIR + "/tests/data/subducting_plate_different_angles_cartesian.wb";
  WorldBuilder::World world2(file_name2);

  position = {{250e3,500e3,800e3}};
  CHECK(world2.temperature(position, 0, 10) == Approx(1599.9999999994));
  CHECK(world2.composition(position, 0, 0) == 1.0);
  CHECK(world2.temperature(position, 1, 10) == Approx(1586.7321267113));
  CHECK(world2.composition(position, 1, 0) == 1.0);
  CHECK(world2.temperature(position, 1e3, 10) == Approx(233.1298084984));
  CHECK(world2.composition(position, 1e3, 0) == 1.0);
  CHECK(world2.temperature(position, 10e3, 10) == Approx(412.3206666645));
  CHECK(world2.composition(position, 10e3, 0) == 1.0);
  CHECK(world2.temperature(position, 20e3, 10) == Approx(544.1584323159));
  CHECK(world2.composition(position, 20e3, 0) == 1.0);
  CHECK(world2.temperature(position, 40e3, 10) == Approx(814.1198929245));
  CHECK(world2.composition(position, 40e3, 0) == 1.0);
  CHECK(world2.temperature(position, 60e3, 10) == Approx(1087.9994157551));
  CHECK(world2.composition(position, 60e3, 0) == 1.0);
  CHECK(world2.temperature(position, 80e3, 10) == Approx(1365.1437340828));
  CHECK(world2.composition(position, 80e3, 0) == 1.0);
  CHECK(world2.temperature(position, 100e3, 10) == Approx(1645.4330950743));
  CHECK(world2.composition(position, 100e3, 0) == 1.0);


  file_name = WorldBuilder::Data::WORLD_BUILDER_SOURCE_DIR + "/tests/data/subducting_plate_different_angles_cartesian_2.wb";
  WorldBuilder::World world4(file_name);

  position = {{250e3,500e3,800e3}};
  CHECK(world4.temperature(position, 0, 10) == Approx(-1));
  CHECK(world4.composition(position, 0, 0) == 1.0);
  CHECK(world4.composition(position, 0, 1) == 0.0);
  CHECK(world4.composition(position, 0, 2) == 0.0);
  CHECK(world4.composition(position, 0, 3) == 0.0);
  CHECK(world4.temperature(position, 1, 10) == Approx(-1));
  CHECK(world4.composition(position, 1, 0) == 1.0);
  CHECK(world4.composition(position, 1, 1) == 0.0);
  CHECK(world4.composition(position, 1, 2) == 0.0);
  CHECK(world4.composition(position, 1, 3) == 0.0);
  CHECK(world4.temperature(position, 1e3, 10) == Approx(-1));
  CHECK(world4.composition(position, 1e3, 0) == 1.0);
  CHECK(world4.composition(position, 1e3, 1) == 0.0);
  CHECK(world4.composition(position, 1e3, 2) == 0.0);
  CHECK(world4.composition(position, 1e3, 3) == 0.0);
  CHECK(world4.temperature(position, 10e3, 10) == Approx(-1));
  CHECK(world4.composition(position, 10e3, 0) == 1.0);
  CHECK(world4.composition(position, 10e3, 1) == 0.0);
  CHECK(world4.composition(position, 10e3, 2) == 0.0);
  CHECK(world4.composition(position, 10e3, 3) == 0.0);
  CHECK(world4.temperature(position, 20e3, 10) == Approx(573.5111391079));
  CHECK(world4.composition(position, 20e3, 0) == 0.0);
  CHECK(world4.composition(position, 20e3, 1) == 1.0);
  CHECK(world4.composition(position, 20e3, 2) == 0.0);
  CHECK(world4.composition(position, 20e3, 3) == 0.0);
  CHECK(world4.temperature(position, 30e3, 10) == Approx(1374.5429651305));
  CHECK(world4.composition(position, 30e3, 0) == 0.0);
  CHECK(world4.composition(position, 30e3, 1) == 1.0);
  CHECK(world4.composition(position, 30e3, 2) == 0.0);
  CHECK(world4.composition(position, 30e3, 3) == 0.0);
  CHECK(world4.temperature(position, 35e3, 10) == Approx(-3));
  CHECK(world4.composition(position, 35e3, 0) == 0.0);
  CHECK(world4.composition(position, 35e3, 1) == 0.0);
  CHECK(world4.composition(position, 35e3, 2) == 0.25);
  CHECK(world4.composition(position, 35e3, 3) == 0.75);
  CHECK(world4.temperature(position, 40e3, 10) == Approx(-3));
  CHECK(world4.composition(position, 40e3, 0) == 0.0);
  CHECK(world4.composition(position, 40e3, 1) == 0.0);
  CHECK(world4.composition(position, 40e3, 2) == 0.25);
  CHECK(world4.composition(position, 40e3, 3) == 0.75);
  CHECK(world4.temperature(position, 45e3, 10) == Approx(-3));
  CHECK(world4.composition(position, 45e3, 0) == 0.0);
  CHECK(world4.composition(position, 45e3, 1) == 0.0);
  CHECK(world4.composition(position, 45e3, 2) == 0.25);
  CHECK(world4.composition(position, 45e3, 3) == 0.75);
  CHECK(world4.temperature(position, 50e3, 10) == Approx(1622.5575343016));
  CHECK(world4.composition(position, 50e3, 0) == 0.0);
  CHECK(world4.composition(position, 50e3, 1) == 0.0);
  CHECK(world4.composition(position, 50e3, 2) == 0.0);
  CHECK(world4.composition(position, 50e3, 3) == 0.0);
  CHECK(world4.temperature(position, 60e3, 10) == Approx(1627.1070617637));
  CHECK(world4.composition(position, 60e3, 0) == 0.0);
  CHECK(world4.composition(position, 60e3, 1) == 0.0);
  CHECK(world4.composition(position, 60e3, 2) == 0.0);
  CHECK(world4.composition(position, 60e3, 3) == 0.0);
  CHECK(world4.temperature(position, 80e3, 10) == Approx(1636.2444220394));
  CHECK(world4.composition(position, 80e3, 0) == 0.0);
  CHECK(world4.composition(position, 80e3, 1) == 0.0);
  CHECK(world4.composition(position, 80e3, 2) == 0.0);
  CHECK(world4.composition(position, 80e3, 3) == 0.0);
  CHECK(world4.temperature(position, 100e3, 10) == Approx(1645.4330950743));
  CHECK(world4.composition(position, 100e3, 0) == 0.0);
  CHECK(world4.composition(position, 100e3, 1) == 0.0);
  CHECK(world4.composition(position, 100e3, 2) == 0.0);
  CHECK(world4.composition(position, 100e3, 3) == 0.0);

}

TEST_CASE("WorldBuilder Features: Fault")
{
  // Cartesian
  std::string file_name = WorldBuilder::Data::WORLD_BUILDER_SOURCE_DIR + "/tests/data/fault_constant_angles_cartesian.wb";
  WorldBuilder::World world1(file_name);

  // Check continental plate directly (upper case should automatically turn into lower case).
  std::unique_ptr<Features::Interface> fault = Features::Interface::create("fault", &world1);

  // Check fault plate through the world
  std::array<double,3> position = {{0,0,800e3}};
  CHECK(world1.temperature(position, 0, 10) == Approx(1600));
  CHECK(world1.temperature(position, 220e3, 10) == Approx(1701.6589518333));
  CHECK(world1.temperature(position, 230e3, 10) == Approx(965.8099855202));
  CHECK(world1.composition(position, 0, 0) == 0.0);
  CHECK(world1.composition(position, 0, 1) == 0.0);
  CHECK(world1.composition(position, 0, 2) == 0.0);
  CHECK(world1.composition(position, 0, 3) == 0.0);
  CHECK(world1.composition(position, 0, 4) == 0.0);
  CHECK(world1.composition(position, 0, 5) == 0.0);
  CHECK(world1.composition(position, 0, 6) == 0.0);

  position = {{250e3,500e3,800e3}};
  CHECK(world1.temperature(position, 0, 10) == Approx(1600.0));
  CHECK(world1.temperature(position, 10, 10) == Approx(150));
  CHECK(world1.temperature(position, std::sqrt(2) * 50e3 - 1, 10) == Approx(150.0));
  CHECK(world1.temperature(position, std::sqrt(2) * 50e3 + 1, 10) == Approx(1631.9945206949));
  CHECK(world1.composition(position, 0, 0) == 0.0);
  CHECK(world1.composition(position, 0, 1) == 0.0);
  CHECK(world1.composition(position, 0, 2) == 0.0);
  CHECK(world1.composition(position, 0, 3) == 0.0);
  CHECK(world1.composition(position, 10, 0) == 0.0);
  CHECK(world1.composition(position, 10, 1) == 0.0);
  CHECK(world1.composition(position, 10, 2) == 0.0);
  CHECK(world1.composition(position, 10, 3) == 0.25);
  CHECK(world1.composition(position, 10, 4) == 0.75);
  CHECK(world1.composition(position, std::sqrt(2) * 50e3 - 1, 3) == 0.25);
  CHECK(world1.composition(position, std::sqrt(2) * 50e3 - 1, 4) == 0.75);
  CHECK(world1.composition(position, std::sqrt(2) * 50e3 + 1, 3) == 0.0);
  CHECK(world1.composition(position, std::sqrt(2) * 50e3 + 1, 4) == 0.0);
  CHECK(world1.composition(position, 0, 4) == 0.0);
  CHECK(world1.composition(position, 0, 5) == 0.0);
  CHECK(world1.composition(position, 0, 6) == 0.0);


  position = {{250e3,250e3,800e3}};
  CHECK(world1.temperature(position, 0, 10) == Approx(1600.0));
  CHECK(world1.temperature(position, 1, 10) == Approx(293.1595620855));
  CHECK(world1.temperature(position, 5, 10) == Approx(293.1978104273));
  CHECK(world1.temperature(position, 10, 10) == Approx(293.2456208547));
  CHECK(world1.temperature(position, 100, 10) == Approx(294.1062085466));
  CHECK(world1.temperature(position, 500, 10) == Approx(297.9310427331));
  CHECK(world1.temperature(position, 1000, 10) == Approx(302.7120854661));
  CHECK(world1.temperature(position, 5000, 10) == Approx(340.9604273305));
  CHECK(world1.temperature(position, std::sqrt(2) * 50e3/2, 10) == Approx(631.2207737686));
  CHECK(world1.temperature(position, std::sqrt(2) * 50e3 - 1, 10) == Approx(969.2819854517));
  CHECK(world1.temperature(position, std::sqrt(2) * 50e3 + 1, 10) == Approx(1631.9945206949));
  CHECK(world1.composition(position, 0, 0) == 0.0);
  CHECK(world1.composition(position, 0, 1) == 0.0);
  CHECK(world1.composition(position, 0, 2) == 0.0);
  CHECK(world1.composition(position, 0, 3) == 0.0);
  CHECK(world1.composition(position, 10, 0) == 1.0);
  CHECK(world1.composition(position, 10, 1) == 0.0);
  CHECK(world1.composition(position, 10, 2) == 0.0);
  CHECK(world1.composition(position, 10, 3) == 0.0);

  CHECK(world1.composition(position, std::sqrt(2) * 33e3 * 0.5 - 1, 0) == 1.0);
  CHECK(world1.composition(position, std::sqrt(2) * 33e3 * 0.5 - 1, 1) == 0.0);
  CHECK(world1.composition(position, std::sqrt(2) * 33e3 * 0.5 - 1, 2) == 0.0);
  CHECK(world1.composition(position, std::sqrt(2) * 33e3 * 0.5 + 1, 0) == 0.0);
  CHECK(world1.composition(position, std::sqrt(2) * 33e3 * 0.5 + 1, 1) == 1.0);
  CHECK(world1.composition(position, std::sqrt(2) * 33e3 * 0.5 + 1, 2) == 0.0);
  CHECK(world1.composition(position, std::sqrt(2) * 66e3 * 0.5 - 1, 1) == 1.0);
  CHECK(world1.composition(position, std::sqrt(2) * 66e3 * 0.5 + 1, 1) == 0.0);
  CHECK(world1.composition(position, std::sqrt(2) * 99e3 * 0.5 - 1, 2) == 0.25);
  CHECK(world1.composition(position, std::sqrt(2) * 99e3 * 0.5 + 1, 2) == 0.0);
  CHECK(world1.composition(position, std::sqrt(2) * 99e3 * 0.5 - 1, 3) == 0.75);
  CHECK(world1.composition(position, std::sqrt(2) * 99e3 * 0.5 + 1, 3) == 0.0);
  CHECK(world1.composition(position, std::sqrt(2) * 100e3 * 0.5 - 1, 3) == 0.0);
  CHECK(world1.composition(position, std::sqrt(2) * 100e3 * 0.5 + 1, 3) == 0.0);
  CHECK(world1.composition(position, 0, 4) == 0.0);
  CHECK(world1.composition(position, 0, 5) == 0.0);
  CHECK(world1.composition(position, 0, 6) == 0.0);

  position = {{250e3,250e3,800e3}};
  CHECK(world1.composition(position, 1, 0) == 1.0);
  CHECK(world1.composition(position, 1, 1) == 0.0);
  CHECK(world1.composition(position, 1, 2) == 0.0);
  position = {{250e3,250e3-std::sqrt(2) * 33e3 * 0.5 + 1, 800e3}};
  CHECK(world1.composition(position, 1, 0) == 1.0);
  CHECK(world1.composition(position, 1, 1) == 0.0);
  CHECK(world1.composition(position, 1, 2) == 0.0);
  position = {{250e3,250e3-std::sqrt(2) * 66e3 * 0.5 + 1, 800e3}};
  CHECK(world1.composition(position, 1, 0) == 0.0);
  CHECK(world1.composition(position, 1, 1) == 1.0);
  CHECK(world1.composition(position, 1, 2) == 0.0);
  position = {{250e3,250e3-std::sqrt(2) * 99e3 * 0.5 + 1, 800e3}};
  CHECK(world1.composition(position, 1, 0) == 0.0);
  CHECK(world1.composition(position, 1, 1) == 0.0);
  CHECK(world1.composition(position, 1, 2) == 0.25);
  CHECK(world1.composition(position, 1, 3) == 0.75);

  std::string file_name2 = WorldBuilder::Data::WORLD_BUILDER_SOURCE_DIR + "/tests/data/fault_constant_angles_cartesian_force_temp.wb";
  WorldBuilder::World world2(file_name2);

  // Check fault plate through the world
  position = {{0,0,800e3}};
  CHECK(world2.temperature(position, 0, 10) == Approx(293.15));
  CHECK(world2.temperature(position, 220e3, 10) == Approx(1701.6589518333));
  CHECK(world2.temperature(position, 230e3, 10) == Approx(100.0));
  CHECK(world2.composition(position, 0, 0) == 0.0);
  CHECK(world2.composition(position, 0, 1) == 0.0);
  CHECK(world2.composition(position, 0, 2) == 0.0);
  CHECK(world2.composition(position, 0, 3) == 0.0);
  CHECK(world2.composition(position, 0, 4) == 0.0);
  CHECK(world2.composition(position, 0, 5) == 0.0);
  CHECK(world2.composition(position, 0, 6) == 0.0);

  position = {{250e3,500e3,800e3}};
  CHECK(world2.temperature(position, 0, 10) == Approx(293.15));
  CHECK(world2.temperature(position, 10, 10) == Approx(150));
  CHECK(world2.temperature(position, std::sqrt(2) * 50e3 - 1, 10) == Approx(150.0));
  CHECK(world2.temperature(position, std::sqrt(2) * 50e3 + 1, 10) == Approx(1631.9945206949));
  CHECK(world2.composition(position, 0, 0) == 0.0);
  CHECK(world2.composition(position, 0, 1) == 0.0);
  CHECK(world2.composition(position, 0, 2) == 0.0);
  CHECK(world2.composition(position, 0, 3) == 0.0);
  CHECK(world2.composition(position, 10, 0) == 0.0);
  CHECK(world2.composition(position, 10, 1) == 0.0);
  CHECK(world2.composition(position, 10, 2) == 0.0);
  CHECK(world1.composition(position, 10, 3) == 0.25);
  CHECK(world1.composition(position, 10, 4) == 0.75);
  CHECK(world2.composition(position, std::sqrt(2) * 50e3 - 1, 3) == 0.25);
  CHECK(world2.composition(position, std::sqrt(2) * 50e3 - 1, 4) == 0.75);
  CHECK(world2.composition(position, std::sqrt(2) * 50e3 + 1, 3) == 0.0);
  CHECK(world2.composition(position, std::sqrt(2) * 50e3 + 1, 4) == 0.0);
  CHECK(world2.composition(position, 0, 4) == 0.0);
  CHECK(world2.composition(position, 0, 5) == 0.0);
  CHECK(world2.composition(position, 0, 6) == 0.0);


  position = {{250e3,250e3,800e3}};
  CHECK(world2.temperature(position, 0, 10) == Approx(293.15));
  CHECK(world2.temperature(position, 1, 10) == Approx(100.0));
  CHECK(world2.temperature(position, 5, 10) == Approx(100.0));
  CHECK(world2.temperature(position, 10, 10) == Approx(100.0));
  CHECK(world2.temperature(position, 100, 10) == Approx(100.0));
  CHECK(world2.temperature(position, 500, 10) == Approx(100.0));
  CHECK(world2.temperature(position, 1000, 10) == Approx(100.0));
  CHECK(world2.temperature(position, 5000, 10) == Approx(100.0));
  CHECK(world2.temperature(position, std::sqrt(2) * 50e3/2, 10) == Approx(100.0));
  CHECK(world2.temperature(position, std::sqrt(2) * 50e3 - 1, 10) == Approx(100.0));
  CHECK(world2.temperature(position, std::sqrt(2) * 50e3 + 1, 10) == Approx(1631.9945206949));
  CHECK(world2.composition(position, 0, 0) == 0.0);
  CHECK(world2.composition(position, 0, 1) == 0.0);
  CHECK(world2.composition(position, 0, 2) == 0.0);
  CHECK(world2.composition(position, 0, 3) == 0.0);
  CHECK(world2.composition(position, 10, 0) == 1.0);
  CHECK(world2.composition(position, 10, 1) == 0.0);
  CHECK(world2.composition(position, 10, 2) == 0.0);
  CHECK(world2.composition(position, 10, 3) == 0.0);

  CHECK(world2.composition(position, std::sqrt(2) * 33e3 * 0.5 - 1, 0) == 1.0);
  CHECK(world2.composition(position, std::sqrt(2) * 33e3 * 0.5 - 1, 1) == 0.0);
  CHECK(world2.composition(position, std::sqrt(2) * 33e3 * 0.5 - 1, 2) == 0.0);
  CHECK(world2.composition(position, std::sqrt(2) * 33e3 * 0.5 + 1, 0) == 0.0);
  CHECK(world2.composition(position, std::sqrt(2) * 33e3 * 0.5 + 1, 1) == 1.0);
  CHECK(world2.composition(position, std::sqrt(2) * 33e3 * 0.5 + 1, 2) == 0.0);
  CHECK(world2.composition(position, std::sqrt(2) * 66e3 * 0.5 - 1, 1) == 1.0);
  CHECK(world2.composition(position, std::sqrt(2) * 66e3 * 0.5 + 1, 1) == 0.0);
  CHECK(world2.composition(position, std::sqrt(2) * 99e3 * 0.5 - 1, 2) == 0.25);
  CHECK(world2.composition(position, std::sqrt(2) * 99e3 * 0.5 + 1, 2) == 0.0);
  CHECK(world2.composition(position, std::sqrt(2) * 99e3 * 0.5 - 1, 3) == 0.75);
  CHECK(world2.composition(position, std::sqrt(2) * 99e3 * 0.5 + 1, 3) == 0.0);
  CHECK(world2.composition(position, std::sqrt(2) * 100e3 * 0.5 - 1, 3) == 0.0);
  CHECK(world2.composition(position, std::sqrt(2) * 100e3 * 0.5 + 1, 3) == 0.0);
  CHECK(world2.composition(position, 0, 4) == 0.0);
  CHECK(world2.composition(position, 0, 5) == 0.0);
  CHECK(world2.composition(position, 0, 6) == 0.0);

  position = {{250e3,250e3,800e3}};
  CHECK(world2.composition(position, 1, 0) == 1.0);
  CHECK(world2.composition(position, 1, 1) == 0.0);
  CHECK(world2.composition(position, 1, 2) == 0.0);
  position = {{250e3,250e3-std::sqrt(2) * 33e3 * 0.5 + 1, 800e3}};
  CHECK(world2.composition(position, 1, 0) == 1.0);
  CHECK(world2.composition(position, 1, 1) == 0.0);
  CHECK(world2.composition(position, 1, 2) == 0.0);
  position = {{250e3,250e3-std::sqrt(2) * 66e3 * 0.5 + 1, 800e3}};
  CHECK(world2.composition(position, 1, 0) == 0.0);
  CHECK(world2.composition(position, 1, 1) == 1.0);
  CHECK(world2.composition(position, 1, 2) == 0.0);
  position = {{250e3,250e3-std::sqrt(2) * 99e3 * 0.5 + 1, 800e3}};
  CHECK(world2.composition(position, 1, 0) == 0.0);
  CHECK(world2.composition(position, 1, 1) == 0.0);
  CHECK(world2.composition(position, 1, 2) == 0.25);
  CHECK(world2.composition(position, 1, 3) == 0.75);


  // Cartesian
  file_name = WorldBuilder::Data::WORLD_BUILDER_SOURCE_DIR + "/tests/data/fault_constant_angles_cartesian_2.wb";
  WorldBuilder::World world3(file_name);

  // Check fault directly (upper case should automatically turn into lower case).
  std::unique_ptr<Features::Interface> continental_plate = Features::Interface::create("Fault", &world3);

  // Check fault through the world
  position = {{0,0,800e3}};
  CHECK(world3.temperature(position, 0, 10) == Approx(1600.0));
  CHECK(world3.temperature(position, 240e3, 10) == Approx(1711.2149738521));
  CHECK(world3.temperature(position, 260e3, 10) == Approx(1720.8246597128));
  CHECK(world3.composition(position, 0, 0) == 0.0);
  CHECK(world3.composition(position, 0, 1) == 0.0);
  CHECK(world3.composition(position, 0, 2) == 0.0);
  CHECK(world3.composition(position, 0, 3) == 0.0);
  CHECK(world3.composition(position, 0, 4) == 0.0);
  CHECK(world3.composition(position, 0, 5) == 0.0);
  CHECK(world3.composition(position, 0, 6) == 0.0);

  position = {{250e3,500e3,800e3}};
  //adibatic temperature
  CHECK(world3.temperature(position, 0, 10) == Approx(1600));
  CHECK(world3.temperature(position, 1, 10) == Approx(1600.0004480001));
  CHECK(world3.temperature(position, 5000, 10) == Approx(1602.241568732));
  CHECK(world3.temperature(position, 10e3, 10) == Approx(1604.486277858));
  CHECK(world3.temperature(position, 25e3, 10) == Approx(1611.239291627));
  CHECK(world3.temperature(position, 50e3, 10) == Approx(1622.5575343016));
  CHECK(world3.temperature(position, 75e3, 10) == Approx(1633.95528262));
  CHECK(world3.temperature(position, 150e3, 10) == Approx(1668.6311660012));
  //CHECK(world3.temperature(position, std::sqrt(2) * 100e3 - 1, 10) == Approx(150.0));
  //CHECK(world3.temperature(position, std::sqrt(2) * 100e3 + 1, 10) == Approx(1664.6283561404));
  CHECK(world3.composition(position, 0, 0) == 0.0);
  CHECK(world3.composition(position, 0, 1) == 0.0);
  CHECK(world3.composition(position, 0, 2) == 0.0);
  CHECK(world3.composition(position, 0, 3) == 0.0);
  CHECK(world3.composition(position, 0, 4) == 0.0);
  CHECK(world3.composition(position, 10, 0) == 1.0);
  CHECK(world3.composition(position, 10, 1) == 0.0);
  CHECK(world3.composition(position, 10, 2) == 0.0);
  CHECK(world3.composition(position, 10, 3) == 0.0);
  CHECK(world3.composition(position, 10, 4) == 0.0);
  //todo: recheck these results
  CHECK(world3.composition(position, std::sqrt(2) * 16.5e3 - 1, 0) == 1.0);
  CHECK(world3.composition(position, std::sqrt(2) * 16.5e3 + 1, 0) == 1.0);
  CHECK(world3.composition(position, std::sqrt(2) * 16.5e3 - 1, 1) == 0.0);
  CHECK(world3.composition(position, std::sqrt(2) * 16.5e3 + 1, 1) == 0.0);
  CHECK(world3.composition(position, std::sqrt(2) * 33e3 - 1, 1) == 0.0);
  CHECK(world3.composition(position, std::sqrt(2) * 33e3 + 1, 1) == 1.0);
  CHECK(world3.composition(position, std::sqrt(2) * 33e3 - 1, 2) == 0.0);
  CHECK(world3.composition(position, std::sqrt(2) * 33e3 + 1, 2) == 0.0);
  CHECK(world3.composition(position, std::sqrt(2) * 33e3 + 1, 3) == 0.0);
  CHECK(world3.composition(position, std::sqrt(2) * 49.5e3 - 1, 2) == 0.0);
  CHECK(world3.composition(position, std::sqrt(2) * 49.5e3 - 1, 3) == 0.0);
  CHECK(world3.composition(position, std::sqrt(2) * 49.5e3 + 1, 2) == 0.0);
  // this comes form the first subducting plate
  CHECK(world3.composition(position, std::sqrt(2) * 49.5e3 + 1, 3) == 0.0);
  CHECK(world3.composition(position, std::sqrt(2) * 50e3 - 1, 3) == 0.0);
  CHECK(world3.composition(position, std::sqrt(2) * 50e3 + 1, 3) == 0.0);
  CHECK(world3.composition(position, 0, 4) == 0.0);
  CHECK(world3.composition(position, 0, 5) == 0.0);
  CHECK(world3.composition(position, 0, 6) == 0.0);


  position = {{250e3,600e3,800e3}};
  CHECK(world3.temperature(position, 0, 10) == Approx(1600));
  CHECK(world3.temperature(position, 10, 10) == Approx(1600.0044800063));
  CHECK(world3.temperature(position, 100e3, 10) == Approx(1.0));
  CHECK(world3.temperature(position, 100e3+1, 10) == Approx(1.0000424264));
  CHECK(world3.temperature(position, 101e3, 10) == Approx(1.0424264069));
  CHECK(world3.temperature(position, 110e3, 10) == Approx(1.4242640687));
  CHECK(world3.temperature(position, 150e3, 10) == Approx(3.1213203436));
  CHECK(world3.temperature(position, 200e3, 10) == Approx(1692.1562939786));
  CHECK(world3.composition(position, 0, 0) == 0.0);
  CHECK(world3.composition(position, 0, 1) == 0.0);
  CHECK(world3.composition(position, 0, 2) == 0.0);
  CHECK(world3.composition(position, 0, 3) == 0.0);
  CHECK(world3.composition(position, 10, 0) == 0.0);
  CHECK(world3.composition(position, 10, 1) == 0.0);
  CHECK(world3.composition(position, 10, 2) == 0.0);
  CHECK(world3.composition(position, 10, 3) == 0.0);
  CHECK(world3.composition(position, 100e3, 0) == Approx(0.0));
  CHECK(world3.composition(position, 100e3, 1) == Approx(0.0));
  CHECK(world3.composition(position, 100e3, 2) == Approx(0.0));
  CHECK(world3.composition(position, 100e3, 3) == Approx(1.0));
  CHECK(world3.composition(position, 100e3, 4) == Approx(0.0));
  CHECK(world3.composition(position, 100e3+1, 0) == Approx(0.0));
  CHECK(world3.composition(position, 100e3+1, 1) == Approx(0.0));
  CHECK(world3.composition(position, 100e3+1, 2) == Approx(0.0));
  CHECK(world3.composition(position, 100e3+1, 3) == Approx(1.0));
  CHECK(world3.composition(position, 100e3+1, 4) == Approx(0.0));
  CHECK(world3.composition(position, 101e3, 0) == Approx(0.0));
  CHECK(world3.composition(position, 101e3+1, 1) == Approx(0.0));
  CHECK(world3.composition(position, 101e3+1, 2) == Approx(0.0));
  CHECK(world3.composition(position, 101e3+1, 3) == Approx(1.0));
  CHECK(world3.composition(position, 101e3+1, 4) == Approx(0.0));
  CHECK(world3.composition(position, 150e3, 0) == Approx(0.0));
  CHECK(world3.composition(position, 150e3, 1) == Approx(0.0));
  CHECK(world3.composition(position, 150e3, 2) == Approx(0.0));
  CHECK(world3.composition(position, 150e3, 3) == Approx(1.0));
  CHECK(world3.composition(position, 150e3, 4) == Approx(0.0));
  CHECK(world3.composition(position, 200e3, 0) == Approx(0.0));
  CHECK(world3.composition(position, 200e3, 1) == Approx(0.0));
  CHECK(world3.composition(position, 200e3, 2) == Approx(0.0));
  CHECK(world3.composition(position, 200e3, 3) == Approx(0.0));
  CHECK(world3.composition(position, 200e3, 4) == Approx(0.0));
  CHECK(world3.composition(position, 0, 4) == 0.0);
  CHECK(world3.composition(position, 0, 5) == 0.0);
  CHECK(world3.composition(position, 0, 6) == 0.0);

  position = {{650e3,650e3,800e3}};
  CHECK(world3.temperature(position, 0, 10) == Approx(4.3590710784));
  CHECK(world3.temperature(position, 10, 10) == Approx(4.3590710784));
  CHECK(world3.temperature(position, 100e3, 10) == Approx(4.3590710784));
  CHECK(world3.temperature(position, 100e3+1, 10) == Approx(4.3590710784));
  CHECK(world3.temperature(position, 101e3, 10) == Approx(4.3590710784));
  CHECK(world3.temperature(position, 110e3, 10) == Approx(4.3590710784));
  CHECK(world3.temperature(position, 150e3, 10) == Approx(1668.6311660012));
  CHECK(world3.temperature(position, 200e3, 10) == Approx(1692.1562939786));
  CHECK(world3.composition(position, 0, 0) == 0.0);
  CHECK(world3.composition(position, 0, 1) == 0.0);
  CHECK(world3.composition(position, 0, 2) == 0.0);
  CHECK(world3.composition(position, 0, 3) == Approx(0.0897677696));
  CHECK(world3.composition(position, 10, 0) == 0.0);
  CHECK(world3.composition(position, 10, 1) == 0.0);
  CHECK(world3.composition(position, 10, 2) == 0.0);
  CHECK(world3.composition(position, 10, 3) == Approx(0.0897677696));
  CHECK(world3.composition(position, 100e3, 0) == Approx(0.0));
  CHECK(world3.composition(position, 100e3, 1) == Approx(0.0));
  CHECK(world3.composition(position, 100e3, 2) == Approx(0.0));
  CHECK(world3.composition(position, 100e3, 3) == Approx(0.0897677696));
  CHECK(world3.composition(position, 100e3, 4) == Approx(0.2693033088));
  CHECK(world3.composition(position, 100e3+1, 0) == Approx(0.0));
  CHECK(world3.composition(position, 100e3+1, 1) == Approx(0.0));
  CHECK(world3.composition(position, 100e3+1, 2) == Approx(0.0));
  CHECK(world3.composition(position, 100e3+1, 3) == Approx(0.0897677696));
  CHECK(world3.composition(position, 100e3+1, 4) == Approx(0.2693033088));
  CHECK(world3.composition(position, 101e3, 0) == Approx(0.0));
  CHECK(world3.composition(position, 101e3+1, 1) == Approx(0.0));
  CHECK(world3.composition(position, 101e3+1, 2) == Approx(0.0));
  CHECK(world3.composition(position, 101e3+1, 3) == Approx(0.0897677696));
  CHECK(world3.composition(position, 101e3+1, 4) == Approx(0.2693033088));
  CHECK(world3.composition(position, 150e3, 0) == Approx(0.0));
  CHECK(world3.composition(position, 150e3, 1) == Approx(0.0));
  CHECK(world3.composition(position, 150e3, 2) == Approx(0.0));
  CHECK(world3.composition(position, 150e3, 3) == Approx(0.0));
  CHECK(world3.composition(position, 150e3, 4) == Approx(0.0));
  CHECK(world3.composition(position, 200e3, 0) == Approx(0.0));
  CHECK(world3.composition(position, 200e3, 1) == Approx(0.0));
  CHECK(world3.composition(position, 200e3, 2) == Approx(0.0));
  CHECK(world3.composition(position, 200e3, 3) == Approx(0.0));
  CHECK(world3.composition(position, 200e3, 4) == Approx(0.0));
  CHECK(world3.composition(position, 0, 4) == Approx(0.2693033088));
  CHECK(world3.composition(position, 0, 5) == Approx(0.6409289216));
  CHECK(world3.composition(position, 0, 6) == 0.0);

  position = {{700e3,675e3,800e3}};
  CHECK(world3.temperature(position, 0, 10) == Approx(4.4592312705));
  CHECK(world3.temperature(position, 10, 10) == Approx(4.4592312705));
  CHECK(world3.temperature(position, 100e3, 10) == Approx(4.4592312705));
  CHECK(world3.temperature(position, 100e3+1, 10) == Approx(4.4592312705));
  CHECK(world3.temperature(position, 101e3, 10) == Approx(4.4592312705));
  CHECK(world3.temperature(position, 110e3, 10) == Approx(4.4592312705));
  CHECK(world3.temperature(position, 150e3, 10) == Approx(1668.6311660012));
  CHECK(world3.temperature(position, 200e3, 10) == Approx(1692.1562939786));
  CHECK(world3.composition(position, 0, 0) == 0.0);
  CHECK(world3.composition(position, 0, 1) == 0.0);
  CHECK(world3.composition(position, 0, 2) == 0.0);
  CHECK(world3.composition(position, 0, 3) == Approx(0.1148078176));
  CHECK(world3.composition(position, 10, 0) == 0.0);
  CHECK(world3.composition(position, 10, 1) == 0.0);
  CHECK(world3.composition(position, 10, 2) == 0.0);
  CHECK(world3.composition(position, 10, 3) == Approx(0.1148078176));
  CHECK(world3.composition(position, 100e3, 0) == Approx(0.0));
  CHECK(world3.composition(position, 100e3, 1) == Approx(0.0));
  CHECK(world3.composition(position, 100e3, 2) == Approx(0.0));
  CHECK(world3.composition(position, 100e3, 3) == Approx(0.1148078176));
  CHECK(world3.composition(position, 100e3, 4) == Approx(0.3444234529));
  CHECK(world3.composition(position, 100e3+1, 0) == Approx(0.0));
  CHECK(world3.composition(position, 100e3+1, 1) == Approx(0.0));
  CHECK(world3.composition(position, 100e3+1, 2) == Approx(0.0));
  CHECK(world3.composition(position, 100e3+1, 3) == Approx(0.1148078176));
  CHECK(world3.composition(position, 100e3+1, 4) == Approx(0.3444234529));
  CHECK(world3.composition(position, 101e3, 0) == Approx(0.0));
  CHECK(world3.composition(position, 101e3+1, 1) == Approx(0.0));
  CHECK(world3.composition(position, 101e3+1, 2) == Approx(0.0));
  CHECK(world3.composition(position, 101e3+1, 3) == Approx(0.1148078176));
  CHECK(world3.composition(position, 101e3+1, 4) == Approx(0.3444234529));
  CHECK(world3.composition(position, 150e3, 0) == Approx(0.0));
  CHECK(world3.composition(position, 150e3, 1) == Approx(0.0));
  CHECK(world3.composition(position, 150e3, 2) == Approx(0.0));
  CHECK(world3.composition(position, 150e3, 3) == Approx(0.0));
  CHECK(world3.composition(position, 150e3, 4) == Approx(0.0));
  CHECK(world3.composition(position, 200e3, 0) == Approx(0.0));
  CHECK(world3.composition(position, 200e3, 1) == Approx(0.0));
  CHECK(world3.composition(position, 200e3, 2) == Approx(0.0));
  CHECK(world3.composition(position, 200e3, 3) == Approx(0.0));
  CHECK(world3.composition(position, 200e3, 4) == Approx(0.0));
  CHECK(world3.composition(position, 0, 4) == Approx(0.3444234529));
  CHECK(world3.composition(position, 0, 5) == Approx(0.5407687295));
  CHECK(world3.composition(position, 0, 6) == 0.0);


  file_name = WorldBuilder::Data::WORLD_BUILDER_SOURCE_DIR + "/tests/data/fault_different_angles_cartesian.wb";
  WorldBuilder::World world4(file_name);

  position = {{250e3,501e3,800e3}};
  CHECK(world4.temperature(position, 0, 10) == Approx(-1));
  CHECK(world4.composition(position, 0, 0) == 1.0);
  CHECK(world4.temperature(position, 1, 10) == Approx(-1));
  CHECK(world4.composition(position, 1, 0) == 1.0);
  CHECK(world4.temperature(position, 1e3, 10) == Approx(-1));
  CHECK(world4.composition(position, 1e3, 0) == 1.0);
  CHECK(world4.temperature(position, 10e3, 10) == Approx(-1));
  CHECK(world4.composition(position, 10e3, 0) == 1.0);
  CHECK(world4.temperature(position, 20e3, 10) == Approx(573.2769020077));
  CHECK(world4.composition(position, 20e3, 0) == 0.0);
  CHECK(world4.temperature(position, 30e3, 10) == Approx(1374.2941781431));
  CHECK(world4.composition(position, 30e3, 0) == 0.0);
  CHECK(world4.temperature(position, 35e3, 10) == Approx(-3));
  CHECK(world4.composition(position, 35e3, 0) == 0.0);
  CHECK(world4.temperature(position, 40e3, 10) == Approx(-3));
  CHECK(world4.composition(position, 40e3, 0) == 0.0);
  CHECK(world4.temperature(position, 45e3, 10) == Approx(-3));
  CHECK(world4.composition(position, 45e3, 0) == 0.0);
  CHECK(world4.temperature(position, 50e3, 10) == Approx(1622.5575343016));
  CHECK(world4.composition(position, 50e3, 0) == 0.0);
  CHECK(world4.temperature(position, 60e3, 10) == Approx(1627.1070617637));
  CHECK(world4.composition(position, 60e3, 0) == 0.0);
  CHECK(world4.temperature(position, 80e3, 10) == Approx(1636.2444220394));
  CHECK(world4.composition(position, 80e3, 0) == 0.0);
  CHECK(world4.temperature(position, 100e3, 10) == Approx(1645.4330950743));
  CHECK(world4.composition(position, 100e3, 0) == 0.0);
}

TEST_CASE("WorldBuilder Features: coordinate interpolation")
{
  {
    std::string file_name = WorldBuilder::Data::WORLD_BUILDER_SOURCE_DIR + "/tests/data/interpolation_none_cartesian.wb";
    WorldBuilder::World world1(file_name);

    std::array<double,3> position = {{374e3,875e3,800e3}};
    CHECK(world1.temperature(position, 0, 10) == Approx(150));
    CHECK(world1.composition(position, 0, 0) == 1.0);
    position = {{376e3,875e3,800e3}};
    CHECK(world1.temperature(position, 0, 10) == Approx(1600));
    CHECK(world1.composition(position, 0, 0) == 0.0);

    position = {{375e3,874e3,800e3}};
    CHECK(world1.temperature(position, 0, 10) == Approx(1600));
    CHECK(world1.composition(position, 0, 0) == 0.0);
    position = {{375e3,876e3,800e3}};
    CHECK(world1.temperature(position, 0, 10) == Approx(150));
    CHECK(world1.composition(position, 0, 0) == 1.0);


    position = {{374e3,625e3,800e3}};
    CHECK(world1.temperature(position, 0, 10) == Approx(150));
    CHECK(world1.composition(position, 0, 0) == 1.0);
    position = {{376e3,625e3,800e3}};
    CHECK(world1.temperature(position, 0, 10) == Approx(1600));
    CHECK(world1.composition(position, 0, 0) == 0.0);

    position = {{375e3,624e3,800e3}};
    CHECK(world1.temperature(position, 0, 10) == Approx(150));
    CHECK(world1.composition(position, 0, 0) == 1.0);
    position = {{375e3,626e3,800e3}};
    CHECK(world1.temperature(position, 0, 10) == Approx(1600));
    CHECK(world1.composition(position, 0, 0) == 0.0);


    position = {{638e3,425e3,800e3}};
    CHECK(world1.temperature(position, 1e3, 10) == Approx(150));
    CHECK(world1.composition(position, 1e3, 0) == 1.0);
    position = {{637e3,425e3,800e3}};
    CHECK(world1.temperature(position, 10, 10) == Approx(1600));
    CHECK(world1.composition(position, 10, 0) == 0.0);

    position = {{625e3,200e3,800e3}};
    CHECK(world1.temperature(position, 10, 10) == Approx(150));
    CHECK(world1.composition(position, 10, 0) == 1.0);
    position = {{624e3,200e3,800e3}};
    CHECK(world1.temperature(position, 10, 10) == Approx(1600));
    CHECK(world1.composition(position, 10, 0) == 0.0);
  }

  {
    std::string file_name = WorldBuilder::Data::WORLD_BUILDER_SOURCE_DIR + "/tests/data/interpolation_linear_cartesian.wb";
    WorldBuilder::World world1(file_name);

    std::array<double,3> position = {{374e3,875e3,800e3}};
    CHECK(world1.temperature(position, 0, 10) == Approx(150));
    CHECK(world1.composition(position, 0, 0) == 1.0);
    position = {{376e3,875e3,800e3}};
    CHECK(world1.temperature(position, 0, 10) == Approx(1600));
    CHECK(world1.composition(position, 0, 0) == 0.0);

    position = {{375e3,874e3,800e3}};
    CHECK(world1.temperature(position, 0, 10) == Approx(1600));
    CHECK(world1.composition(position, 0, 0) == 0.0);
    position = {{375e3,876e3,800e3}};
    CHECK(world1.temperature(position, 0, 10) == Approx(150));
    CHECK(world1.composition(position, 0, 0) == 1.0);


    position = {{374e3,625e3,800e3}};
    CHECK(world1.temperature(position, 0, 10) == Approx(150));
    CHECK(world1.composition(position, 0, 0) == 1.0);
    position = {{376e3,625e3,800e3}};
    CHECK(world1.temperature(position, 0, 10) == Approx(1600));
    CHECK(world1.composition(position, 0, 0) == 0.0);

    position = {{375e3,624e3,800e3}};
    CHECK(world1.temperature(position, 0, 10) == Approx(150));
    CHECK(world1.composition(position, 0, 0) == 1.0);
    position = {{375e3,626e3,800e3}};
    CHECK(world1.temperature(position, 0, 10) == Approx(1600));
    CHECK(world1.composition(position, 0, 0) == 0.0);


    position = {{638e3,425e3,800e3}};
    CHECK(world1.temperature(position, 1e3, 10) == Approx(150));
    CHECK(world1.composition(position, 1e3, 0) == 1.0);
    position = {{637e3,425e3,800e3}};
    CHECK(world1.temperature(position, 10, 10) == Approx(1600));
    CHECK(world1.composition(position, 10, 0) == 0.0);

    position = {{625e3,200e3,800e3}};
    CHECK(world1.temperature(position, 10, 10) == Approx(150));
    CHECK(world1.composition(position, 10, 0) == 1.0);
    position = {{624e3,200e3,800e3}};
    CHECK(world1.temperature(position, 10, 10) == Approx(1600));
    CHECK(world1.composition(position, 10, 0) == 0.0);
  }

  {
    std::string file_name = WorldBuilder::Data::WORLD_BUILDER_SOURCE_DIR + "/tests/data/interpolation_monotone_spline_cartesian.wb";
    WorldBuilder::World world1(file_name);

    std::array<double,3> position = {{374e3,875e3,800e3}};
    CHECK(world1.temperature(position, 0, 10) == Approx(1600));
    CHECK(world1.composition(position, 0, 0) == 0.0);
    position = {{376e3,875e3,800e3}};
    CHECK(world1.temperature(position, 0, 10) == Approx(1600));
    CHECK(world1.composition(position, 0, 0) == 0.0);
    position = {{350e3,900e3,800e3}};
    CHECK(world1.temperature(position, 0, 10) == Approx(150));
    CHECK(world1.composition(position, 0, 0) == 1.0);

    position = {{375e3,874e3,800e3}};
    CHECK(world1.temperature(position, 0, 10) == Approx(1600));
    CHECK(world1.composition(position, 0, 0) == 0.0);
    position = {{375e3,876e3,800e3}};
    CHECK(world1.temperature(position, 0, 10) == Approx(1600));
    CHECK(world1.composition(position, 0, 0) == 0.0);


    position = {{374e3,625e3,800e3}};
    CHECK(world1.temperature(position, 0, 10) == Approx(1600));
    CHECK(world1.composition(position, 0, 0) == 0.0);
    position = {{376e3,625e3,800e3}};
    CHECK(world1.temperature(position, 0, 10) == Approx(1600));
    CHECK(world1.composition(position, 0, 0) == 0.0);
    position = {{350e3,600e3,800e3}};
    CHECK(world1.temperature(position, 0, 10) == Approx(150));
    CHECK(world1.composition(position, 0, 0) == 1.0);

    position = {{375e3,624e3,800e3}};
    CHECK(world1.temperature(position, 0, 10) == Approx(1600));
    CHECK(world1.composition(position, 0, 0) == 0.0);
    position = {{375e3,626e3,800e3}};
    CHECK(world1.temperature(position, 0, 10) == Approx(1600));
    CHECK(world1.composition(position, 0, 0) == 0.0);


    position = {{638e3,425e3,800e3}};
    CHECK(world1.temperature(position, 10, 10) == Approx(1600));
    CHECK(world1.composition(position, 10, 0) == 0.0);
    position = {{637e3,425e3,800e3}};
    CHECK(world1.temperature(position, 10, 10) == Approx(1600));
    CHECK(world1.composition(position, 10, 0) == 0.0);
    position = {{617.5e3,445e3,800e3}};
    CHECK(world1.temperature(position, 1e3, 10) == Approx(150));
    CHECK(world1.composition(position, 1e3, 0) == 1.0);

    position = {{625e3,200e3,800e3}};
    CHECK(world1.temperature(position, 10, 10) == Approx(1600));
    CHECK(world1.composition(position, 10, 0) == 0.0);
    position = {{624e3,200e3,800e3}};
    CHECK(world1.temperature(position, 10, 10) == Approx(1600));
    CHECK(world1.composition(position, 10, 0) == 0.0);
    position = {{607e3,180e3,800e3}};
    CHECK(world1.temperature(position, 3e3, 10) == Approx(150));
    CHECK(world1.composition(position, 3e3, 0) == 1.0);
  }
}

TEST_CASE("WorldBuilder Types: Double")
{
#define TYPE Double
  Types::TYPE type(1);
  CHECK(type.default_value == 1);
  CHECK(type.get_type() == Types::type::TYPE);

  Types::TYPE type_copy(type);
  CHECK(type_copy.default_value == 1);
  CHECK(type_copy.get_type() == Types::type::TYPE);

  Types::TYPE type_explicit(3);
  CHECK(type_explicit.default_value == 3);
  CHECK(type_explicit.get_type() == Types::type::TYPE);

  std::unique_ptr<Types::Interface> type_clone = type_explicit.clone();
  Types::TYPE *type_clone_natural = dynamic_cast<Types::TYPE *>(type_clone.get());
  CHECK(type_clone_natural->default_value == 3);
  CHECK(type_clone_natural->get_type() == Types::type::TYPE);
#undef TYPE
}

/*TEST_CASE("WorldBuilder Types: UnsignedInt")
{
#define TYPE UnsignedInt
  Types::TYPE type(1,"test");
  CHECK(type.value == 1);
  CHECK(type.default_value == 1);
  CHECK(type.description == "test");
  CHECK(type.get_type() == Types::type::TYPE);

  Types::TYPE type_copy(type);
  CHECK(type_copy.value == 1);
  CHECK(type_copy.default_value == 1);
  CHECK(type_copy.description == "test");
  CHECK(type_copy.get_type() == Types::type::TYPE);

  Types::TYPE type_explicit(2, 3, "test explicit");
  CHECK(type_explicit.value == 2);
  CHECK(type_explicit.default_value == 3);
  CHECK(type_explicit.description == "test explicit");
  CHECK(type_explicit.get_type() == Types::type::TYPE);

  std::unique_ptr<Types::Interface> type_clone = type_explicit.clone();
  Types::TYPE *type_clone_natural = dynamic_cast<Types::TYPE *>(type_clone.get());
  CHECK(type_clone_natural->value == 2);
  CHECK(type_clone_natural->default_value == 3);
  CHECK(type_clone_natural->description == "test explicit");
  CHECK(type_clone_natural->get_type() == Types::type::TYPE);
#undef TYPE
}*/

TEST_CASE("WorldBuilder Types: String")
{
#define TYPE String
  Types::TYPE type("1","test");
  CHECK(type.default_value == "1");
  CHECK(type.get_type() == Types::type::TYPE);

  Types::TYPE type_copy(type);
  CHECK(type_copy.default_value == "1");
  CHECK(type_copy.get_type() == Types::type::TYPE);

  Types::TYPE type_explicit("2", "3", "test explicit");
  CHECK(type_explicit.default_value == "3");
  CHECK(type_explicit.get_type() == Types::type::TYPE);

  std::unique_ptr<Types::Interface> type_clone = type_explicit.clone();
  Types::TYPE *type_clone_natural = dynamic_cast<Types::TYPE *>(type_clone.get());
  CHECK(type_clone_natural->default_value == "3");
  CHECK(type_clone_natural->get_type() == Types::type::TYPE);
#undef TYPE
}

TEST_CASE("WorldBuilder Types: Point 2d")
{
#define TYPE Point<2>
  Types::TYPE type(TYPE(1,2,cartesian),"test");
  CHECK(type.value[0] == TYPE(1,2,cartesian)[0]);
  CHECK(type.value[1] == TYPE(1,2,cartesian)[1]);
  CHECK(type.default_value[0] == TYPE(1,2,cartesian)[0]);
  CHECK(type.default_value[1] == TYPE(1,2,cartesian)[1]);
  CHECK(type.description == "test");
  CHECK(type.get_type() == Types::type::Point2D);

  Types::TYPE type_copy(type);
  CHECK(type_copy.value[0] == TYPE(1,2,cartesian)[0]);
  CHECK(type_copy.value[1] == TYPE(1,2,cartesian)[1]);
  CHECK(type.default_value[0] == TYPE(1,2,cartesian)[0]);
  CHECK(type.default_value[1] == TYPE(1,2,cartesian)[1]);
  CHECK(type_copy.description == "test");
  CHECK(type_copy.get_type() == Types::type::Point2D);

  Types::TYPE type_explicit(TYPE(3,4,cartesian), TYPE(5,6,cartesian), "test explicit");
  CHECK(type_explicit.value[0] == TYPE(3,4,cartesian)[0]);
  CHECK(type_explicit.value[1] == TYPE(3,4,cartesian)[1]);
  CHECK(type_explicit.default_value[0] == TYPE(5,6,cartesian)[0]);
  CHECK(type_explicit.default_value[1] == TYPE(5,6,cartesian)[1]);
  CHECK(type_explicit.description == "test explicit");
  CHECK(type_explicit.get_type() == Types::type::Point2D);

  std::unique_ptr<Types::Interface> type_clone = type_explicit.clone();
  Types::TYPE *type_clone_natural = dynamic_cast<Types::TYPE *>(type_clone.get());
  CHECK(type_clone_natural->value[0] == TYPE(3,4,cartesian)[0]);
  CHECK(type_clone_natural->value[1] == TYPE(3,4,cartesian)[1]);
  CHECK(type_clone_natural->default_value[0] == TYPE(5,6,cartesian)[0]);
  CHECK(type_clone_natural->default_value[1] == TYPE(5,6,cartesian)[1]);
  CHECK(type_clone_natural->description == "test explicit");
  CHECK(type_clone_natural->get_type() == Types::type::Point2D);



  // Test Point operators

  const TYPE point_array(std::array<double,2> {{1,2}},cartesian);
  const TYPE point_explicit(3,4,cartesian);

  Types::TYPE type_point_array(point_array, point_array, "test array");
  Types::TYPE type_point_explicit(point_explicit, point_explicit, "test array");

  CHECK(type_point_array.value.get_array() == std::array<double,2> {{1,2}});
  CHECK(type_point_explicit.value.get_array() == std::array<double,2> {{3,4}});

  // Test multiply operator
  TYPE point = 2 * type_point_array * 1.0;

  CHECK(point.get_array() == std::array<double,2> {{2,4}});

  // Test dot operator
  CHECK(type_point_array * type_point_explicit == 11);

  // Test add operator
  point = type_point_array + type_point_explicit;

  CHECK(point.get_array() == std::array<double,2> {{4,6}});

  // Test subtract operator
  point = type_point_explicit - type_point_array;

  CHECK(point.get_array() == std::array<double,2> {{2,2}});

  // test the access operator
  CHECK(type_point_array[0] == 1);

  type_point_array[0] = 2;
  CHECK(type_point_array[0] == 2);
#undef TYPE
}

TEST_CASE("WorldBuilder Types: Point 3d")
{
#define TYPE Point<3>
  Types::TYPE type(TYPE(1,2,3,cartesian),"test");
  CHECK(type.value[0] == 1);
  CHECK(type.value[1] == 2);
  CHECK(type.value[2] == 3);
  CHECK(type.default_value[0] == 1);
  CHECK(type.default_value[1] == 2);
  CHECK(type.default_value[2] == 3);
  CHECK(type.description == "test");
  CHECK(type.get_type() == Types::type::Point3D);

  Types::TYPE type_copy(type);
  CHECK(type_copy.value[0] == 1);
  CHECK(type_copy.value[1] == 2);
  CHECK(type_copy.value[2] == 3);
  CHECK(type_copy.default_value[0] == 1);
  CHECK(type_copy.default_value[1] == 2);
  CHECK(type_copy.default_value[2] == 3);
  CHECK(type_copy.description == "test");
  CHECK(type_copy.get_type() == Types::type::Point3D);

  Types::TYPE type_explicit(TYPE(4,5,6,cartesian), TYPE(7,8,9,cartesian), "test explicit");
  CHECK(type_explicit.value[0] == 4);
  CHECK(type_explicit.value[1] == 5);
  CHECK(type_explicit.value[2] == 6);
  CHECK(type_explicit.default_value[0] == 7);
  CHECK(type_explicit.default_value[1] == 8);
  CHECK(type_explicit.default_value[2] == 9);
  CHECK(type_explicit.description == "test explicit");
  CHECK(type_explicit.get_type() == Types::type::Point3D);

  std::unique_ptr<Types::Interface> type_clone = type_explicit.clone();
  Types::TYPE *type_clone_natural = dynamic_cast<Types::TYPE *>(type_clone.get());
  CHECK(type_clone_natural->value[0] == 4);
  CHECK(type_clone_natural->value[1] == 5);
  CHECK(type_clone_natural->value[2] == 6);
  CHECK(type_clone_natural->default_value[0] == 7);
  CHECK(type_clone_natural->default_value[1] == 8);
  CHECK(type_clone_natural->default_value[2] == 9);
  CHECK(type_clone_natural->description == "test explicit");
  CHECK(type_clone_natural->get_type() == Types::type::Point3D);


  // Test Point operators

  const TYPE point_array(std::array<double,3> {{1,2,3}},cartesian);
  const TYPE point_explicit(4,5,6,cartesian);

  Types::TYPE type_point_array(point_array, point_array, "test array");
  Types::TYPE type_point_explicit(point_explicit, point_explicit, "test array");

  CHECK(type_point_array.value.get_array() == std::array<double,3> {{1,2,3}});
  CHECK(type_point_explicit.value.get_array() == std::array<double,3> {{4,5,6}});

  // Test multiply operator
  TYPE point = 2 * type_point_array;

  CHECK(point.get_array() == std::array<double,3> {{2,4,6}});

  // Test multiply operator
  point = type_point_array * 2;

  CHECK(point.get_array() == std::array<double,3> {{2,4,6}});

  // Test dot operator
  CHECK(type_point_array * type_point_explicit == 32);

  // Test add operator
  point = type_point_array + type_point_explicit;

  CHECK(point.get_array() == std::array<double,3> {{5,7,9}});

  // Test subtract operator
  point = type_point_explicit - type_point_array;

  CHECK(point.get_array() == std::array<double,3> {{3,3,3}});

  // test the access operator
  CHECK(type_point_array[0] == 1);

  type_point_array[0] = 2;
  CHECK(type_point_array[0] == 2);

  // const test the access operator
  CHECK(point_array[0] == 1);

#undef TYPE
}
/*
TEST_CASE("WorldBuilder Types: Coordinate System")
{
#define TYPE CoordinateSystem
  Types::TYPE type("1","test");
  CHECK(type.value == nullptr);
  CHECK(type.default_value == "1");
  CHECK(type.description == "test");
  CHECK(type.get_type() == Types::type::TYPE);

  std::unique_ptr<Types::Interface> type_clone = type.clone();
  Types::TYPE *type_clone_natural = dynamic_cast<Types::TYPE *>(type_clone.get());
  CHECK(type_clone_natural->value == nullptr);
  CHECK(type_clone_natural->default_value == "1");
  CHECK(type_clone_natural->description == "test");
  CHECK(type_clone_natural->get_type() == Types::type::TYPE);


  // todo: test the set value function.

#undef TYPE
}*/

// not sure how to unit test this.
TEST_CASE("WorldBuilder Types: PluginSystem")
{
#define TYPE PluginSystem
  Types::TYPE type("test", Features::ContinentalPlate::declare_entries, std::vector<std::string> {{"test required"}}, false);
  CHECK(type.default_value == "test");
  CHECK(type.required_entries[0] == "test required");
  CHECK(type.allow_multiple == false);
  CHECK(type.get_type() == Types::type::TYPE);

  Types::TYPE type_copy(type);
  CHECK(type_copy.default_value == "test");
  CHECK(type_copy.required_entries[0] == "test required");
  CHECK(type_copy.allow_multiple == false);
  CHECK(type_copy.get_type() == Types::type::TYPE);

  std::unique_ptr<Types::Interface> type_clone = type_copy.clone();
  Types::TYPE *type_clone_natural = dynamic_cast<Types::TYPE *>(type_clone.get());
  CHECK(type_clone_natural->default_value == "test");
  CHECK(type_clone_natural->required_entries[0] == "test required");
  CHECK(type_clone_natural->allow_multiple == false);
  CHECK(type_clone_natural->get_type() == Types::type::TYPE);

#undef TYPE
}


// not sure how to unit test this.
TEST_CASE("WorldBuilder Types: Segment Object")
{
#define TYPE Segment
  WorldBuilder::Point<2> thickness(1,2,invalid);
  WorldBuilder::Point<2> top_trucation(3,4,invalid);
  WorldBuilder::Point<2> angle(5,6,invalid);
  Objects::TYPE<Features::FaultModels::Temperature::Interface, Features::FaultModels::Composition::Interface>
  type (1.0, thickness, top_trucation, angle,
        std::vector<std::shared_ptr<Features::FaultModels::Temperature::Interface> >(),
        std::vector<std::shared_ptr<Features::FaultModels::Composition::Interface> >());
  CHECK(type.value_length == 1);
  CHECK(type.value_thickness[0] == 1);
  CHECK(type.value_thickness[1] == 2);
  CHECK(type.value_top_truncation[0] == 3);
  CHECK(type.value_top_truncation[1] == 4);
  CHECK(type.value_angle[0] == 5);
  CHECK(type.value_angle[1] == 6);
  CHECK(type.get_type() == Types::type::TYPE);

  Objects::TYPE<Features::FaultModels::Temperature::Interface, Features::FaultModels::Composition::Interface>
  type_copy(type);
  CHECK(type_copy.value_length == 1);
  CHECK(type_copy.value_thickness[0] == 1);
  CHECK(type_copy.value_thickness[1] == 2);
  CHECK(type_copy.value_top_truncation[0] == 3);
  CHECK(type_copy.value_top_truncation[1] == 4);
  CHECK(type_copy.value_angle[0] == 5);
  CHECK(type_copy.value_angle[1] == 6);
  CHECK(type_copy.get_type() == Types::type::TYPE);

  std::unique_ptr<Types::Interface> type_clone = type_copy.clone();
  Objects::TYPE<Features::FaultModels::Temperature::Interface, Features::FaultModels::Composition::Interface>
  *type_clone_natural = dynamic_cast<Objects::TYPE<Features::FaultModels::Temperature::Interface,
   Features::FaultModels::Composition::Interface> *>(type_clone.get());
  CHECK(type_clone_natural->value_length == 1);
  CHECK(type_clone_natural->value_thickness[0] == 1);
  CHECK(type_clone_natural->value_thickness[1] == 2);
  CHECK(type_clone_natural->value_top_truncation[0] == 3);
  CHECK(type_clone_natural->value_top_truncation[1] == 4);
  CHECK(type_clone_natural->value_angle[0] == 5);
  CHECK(type_clone_natural->value_angle[1] == 6);
  CHECK(type_clone_natural->get_type() == Types::type::TYPE);

#undef TYPE
}

TEST_CASE("WorldBuilder Types: Array")
{
#define TYPE Array
  Types::TYPE type(Types::Double(0));
  CHECK(type.inner_type == Types::type::Double);
  CHECK(type.inner_type_ptr.get() != nullptr);
  CHECK(type.get_type() == Types::type::TYPE);

  /*Types::TYPE type_copy(type);
  CHECK(type_copy.inner_type == Types::type::Double);
  CHECK(type_copy.inner_type_ptr.get() == nullptr);
  CHECK(type_copy.get_type() == Types::type::TYPE);

  Types::TYPE type_explicit(std::vector<unsigned int> {1,2}, Types::type::Double, "array test explicit");
  CHECK(type_explicit.inner_type == Types::type::Double);
  CHECK(type_explicit.inner_type_ptr.get() == nullptr);
  CHECK(type_explicit.inner_type_index.size() == 2);
  CHECK(type_explicit.description == "array test explicit");
  CHECK(type_explicit.get_type() == Types::type::TYPE);*/

  CHECK_THROWS_WITH(type.clone(),Contains("Error: Cloning Arrays is currently not possible."));
  /*
  std::unique_ptr<Types::Interface> type_clone = type_explicit.clone()
  Types::TYPE *type_clone_natural = dynamic_cast<Types::TYPE *>(type_clone.get());
  CHECK(type_clone_natural->inner_type == Types::type::Double);
  CHECK(type_clone_natural->inner_type_ptr.get() == nullptr);
  CHECK(type_clone_natural->inner_type_index.size() == 2);
  CHECK(type_clone_natural->description == "array test explicit");
  CHECK(type_clone_natural->get_type() == Types::type::TYPE);

  Types::TYPE type_copy2(*type_clone_natural);
  CHECK(type_copy2.inner_type == Types::type::Double);
  CHECK(type_copy2.inner_type_ptr.get() == nullptr);
  CHECK(type_copy2.inner_type_index.size() == 2);
  CHECK(type_copy2.description == "array test explicit");
  CHECK(type_copy2.get_type() == Types::type::TYPE);*/

#undef TYPE
}

TEST_CASE("WorldBuilder Types: Object")
{
#define TYPE Object
  Types::TYPE type(std::vector<std::string> {"test1","test2"},true);
  CHECK(type.required.size() == 2);
  CHECK(type.required[0] == "test1");
  CHECK(type.required[1] == "test2");
  CHECK(type.additional_properties == true);
  CHECK(type.get_type() == Types::type::TYPE);

  Types::TYPE type_copy(type);
  CHECK(type_copy.required.size() == 2);
  CHECK(type_copy.required[0] == "test1");
  CHECK(type_copy.required[1] == "test2");
  CHECK(type_copy.additional_properties == true);
  CHECK(type_copy.get_type() == Types::type::TYPE);

  std::unique_ptr<Types::Interface> type_clone = type_copy.clone();
  Types::TYPE *type_clone_natural = dynamic_cast<Types::TYPE *>(type_clone.get());
  CHECK(type_clone_natural->required.size() == 2);
  CHECK(type_clone_natural->required[0] == "test1");
  CHECK(type_clone_natural->required[1] == "test2");
  CHECK(type_clone_natural->additional_properties == true);
  CHECK(type_clone_natural->get_type() == Types::type::TYPE);

  Types::TYPE type_copy2(*type_clone_natural);
  CHECK(type_copy2.required.size() == 2);
  CHECK(type_copy2.required[0] == "test1");
  CHECK(type_copy2.required[1] == "test2");
  CHECK(type_copy2.additional_properties == true);
  CHECK(type_copy2.get_type() == Types::type::TYPE);

#undef TYPE
}


TEST_CASE("WorldBuilder Types: Bool")
{
#define TYPE Bool
  Types::TYPE type(true);
  CHECK(type.default_value == true);
  CHECK(type.get_type() == Types::type::TYPE);


  Types::TYPE type_copy(type);
  CHECK(type_copy.default_value == true);
  CHECK(type_copy.get_type() == Types::type::TYPE);

  std::unique_ptr<Types::Interface> type_clone = type_copy.clone();
  Types::TYPE *type_clone_natural = dynamic_cast<Types::TYPE *>(type_clone.get());
  CHECK(type_clone_natural->default_value == true);
  CHECK(type_clone_natural->get_type() == Types::type::TYPE);

  Types::TYPE type_copy2(*type_clone_natural);
  CHECK(type_copy2.default_value == true);
  CHECK(type_copy2.get_type() == Types::type::TYPE);

#undef TYPE
}

/*
TEST_CASE("WorldBuilder Types: print_tree")
{
  std::string file_name = WorldBuilder::Data::WORLD_BUILDER_SOURCE_DIR + "/tests/data/simple_wb1.json";
  boost::property_tree::ptree tree;
  std::ifstream json_input_stream(file_name.c_str());

  boost::property_tree::json_parser::read_json (json_input_stream, tree);
  std::stringstream output;
  output <<
         "{\n"
         "  \"version\": \"0.1\",\n"
         "  \"coordinate system\": \n"
         "  {\n"
         "    \"model\": \"cartesian\"\n"
         "   },\n"
         "  \"cross section\": \n"
         "  {\n"
         "    \"\": \n"
         "    {\n"
         "      \"\": \"100e3\",\n"
         "      \"\": \"100e3\"\n"
         "     },\n"
         "    \"\": \n"
         "    {\n"
         "      \"\": \"400e3\",\n"
         "      \"\": \"500e3\"\n"
         "     }\n"
         "   },\n"
         "  \"maximum distance between coordinates\": \"5\",\n"
         "  \"features\": \n"
         "  {\n"
         "    \"\": \n"
         "    {\n"
         "      \"model\": \"continental plate\",\n"
         "      \"name\": \"Carribean\",\n"
         "      \"max depth\": \"300e3\",\n"
         "      \"coordinates\": \n"
         "      {\n"
         "        \"\": \n"
         "        {\n"
         "          \"\": \"-1e3\",\n"
         "          \"\": \"500e3\"\n"
         "         },\n"
         "        \"\": \n"
         "        {\n"
         "          \"\": \"500e3\",\n"
         "          \"\": \"500e3\"\n"
         "         },\n"
         "        \"\": \n"
         "        {\n"
         "          \"\": \"500e3\",\n"
         "          \"\": \"1000e3\"\n"
         "         },\n"
         "        \"\": \n"
         "        {\n"
         "          \"\": \"-1e3\",\n"
         "          \"\": \"1000e3\"\n"
         "         }\n"
         "       },\n"
         "      \"temperature models\": \n"
         "      {\n"
         "        \"\": \n"
         "        {\n"
         "          \"model\": \"uniform\",\n"
         "          \"min depth\": \"0\",\n"
         "          \"max depth\": \"250e3\",\n"
         "          \"temperature\": \"150\"\n"
         "         }\n"
         "       }\n"
         "     },\n"
         "    \"\": \n"
         "    {\n"
         "      \"model\": \"continental plate\",\n"
         "      \"name\": \"Rest\",\n"
         "      \"max depth\": \"300e3\",\n"
         "      \"coordinates\": \n"
         "      {\n"
         "        \"\": \n"
         "        {\n"
         "          \"\": \"2000e3\",\n"
         "          \"\": \"2000e3\"\n"
         "         },\n"
         "        \"\": \n"
         "        {\n"
         "          \"\": \"1000e3\",\n"
         "          \"\": \"2000e3\"\n"
         "         },\n"
         "        \"\": \n"
         "        {\n"
         "          \"\": \"1000e3\",\n"
         "          \"\": \"1000e3\"\n"
         "         },\n"
         "        \"\": \n"
         "        {\n"
         "          \"\": \"2000e3\",\n"
         "          \"\": \"1000e3\"\n"
         "         }\n"
         "       },\n"
         "      \"temperature models\": \n"
         "      {\n"
         "        \"\": \n"
         "        {\n"
         "          \"model\": \"uniform\",\n"
         "          \"max depth\": \"250e3\",\n"
         "          \"temperature\": \"20\"\n"
         "         }\n"
         "       },\n"
         "      \"composition models\": \n"
         "      {\n"
         "        \"\": \n"
         "        {\n"
         "          \"model\": \"uniform\",\n"
         "          \"max depth\": \"250e3\",\n"
         "          \"compositions\": \n"
         "          {\n"
         "            \"\": \"2\"\n"
         "           }\n"
         "         }\n"
         "       }\n"
         "     },\n"
         "    \"\": \n"
         "    {\n"
         "      \"model\": \"continental plate\",\n"
         "      \"name\": \"Carribean2\",\n"
         "      \"max depth\": \"300e3\",\n"
         "      \"coordinates\": \n"
         "      {\n"
         "        \"\": \n"
         "        {\n"
         "          \"\": \"-1e3\",\n"
         "          \"\": \"500e3\"\n"
         "         },\n"
         "        \"\": \n"
         "        {\n"
         "          \"\": \"500e3\",\n"
         "          \"\": \"500e3\"\n"
         "         },\n"
         "        \"\": \n"
         "        {\n"
         "          \"\": \"500e3\",\n"
         "          \"\": \"1000e3\"\n"
         "         },\n"
         "        \"\": \n"
         "        {\n"
         "          \"\": \"-1e3\",\n"
         "          \"\": \"1000e3\"\n"
         "         }\n"
         "       },\n"
         "      \"temperature models\": \n"
         "      {\n"
         "        \"\": \n"
         "        {\n"
         "          \"model\": \"uniform\",\n"
         "          \"min depth\": \"0\",\n"
         "          \"max depth\": \"250e3\",\n"
         "          \"temperature\": \"150\"\n"
         "         }\n"
         "       },\n"
         "      \"composition models\": \n"
         "      {\n"
         "        \"\": \n"
         "        {\n"
         "          \"model\": \"uniform\",\n"
         "          \"min depth\": \"0\",\n"
         "          \"max depth\": \"250e3\",\n"
         "          \"compositions\": \n"
         "          {\n"
         "            \"\": \"3\"\n"
         "           }\n"
         "         }\n"
         "       }\n"
         "     }\n"
         "   }\n"
         " }";
  CHECK(Utilities::print_tree(tree, 0) == output.str());
}*/

TEST_CASE("WorldBuilder Parameters")
{
  // First test a world builder file with a cross section defined
  std::string file = WorldBuilder::Data::WORLD_BUILDER_SOURCE_DIR + "/tests/data/type_data.json";
  std::string file_name = WorldBuilder::Data::WORLD_BUILDER_SOURCE_DIR + "/tests/data/subducting_plate_different_angles_spherical.wb";
  WorldBuilder::World world(file_name);

  Parameters prm(world);
  prm.initialize(file_name);

  //prm.load_entry("Coordinate system", false, Types::CoordinateSystem("cartesian","This determines the coordinate system"));

  // Test the UnsignedInt functions
  /*CHECK_THROWS_WITH(prm.load_entry("non existent unsigned int", true, Types::UnsignedInt(1,"description")),
                    Contains("Entry undeclared: non existent unsigned int"));

  #ifndef NDEBUG
  CHECK_THROWS_WITH(prm.get_unsigned_int("non existent unsigned int"),
                    Contains("Could not find entry 'non existent unsigned int' not found. Make sure it is loaded or set"));
  #endif
  CHECK(prm.load_entry("non existent unsigned int", false, Types::UnsignedInt(1,"description")) == false);
  CHECK(prm.get_unsigned_int("non existent unsigned int") == 1);

  prm.set_entry("new unsigned int", Types::UnsignedInt(2,"description"));
  CHECK(prm.get_unsigned_int("new unsigned int") == 2);

  prm.load_entry("unsigned int", true, Types::UnsignedInt(3,"description"));
  CHECK(prm.get_unsigned_int("unsigned int") == 4);*/
// Todo: redo this with the new type system.
  /*
    // Test the Double functions
    CHECK_THROWS_WITH(prm.load_entry("non existent double", true, Types::Double(1,"description")),
                      Contains("Entry undeclared: non existent"));

  #ifndef NDEBUG
    CHECK_THROWS_WITH(prm.get_double("non existent double"),
                      Contains("Could not find entry 'non existent double' not found. Make sure it is loaded or set"));
  #endif
    CHECK(prm.load_entry("non existent double", false, Types::Double(1,"description")) == false);
    CHECK(prm.get_double("non existent double") == 1);

    prm.set_entry("new double", Types::Double(2,"description"));
    CHECK(prm.get_double("new double") == 2);

    prm.load_entry("double", true, Types::Double(3,"description"));
    CHECK(prm.get_double("double") == 1.23456e2);


    // Test the String functions
    CHECK_THROWS_WITH(prm.load_entry("non existent string", true, Types::String("1","description")),
                      Contains("Entry undeclared: non existent string"));

  #ifndef NDEBUG
    CHECK_THROWS_WITH(prm.get_string("non existent string"),
                      Contains("Could not find entry 'non existent string' not found. Make sure it is loaded or set"));
  #endif
    CHECK(prm.load_entry("non exitent string", false, Types::String("1","description")) == false);
    CHECK(prm.get_string("non exitent string") == "1");

    prm.set_entry("new string", Types::String("2","description"));
    CHECK(prm.get_string("new string") == "2");

    prm.load_entry("string", true, Types::String("3","description"));
    CHECK(prm.get_string("string") == "mystring 0");

    // Test the Point functions
    CHECK_THROWS_WITH(prm.load_entry("non existent 2d Point", true, Types::Point<2>(Point<2>(1,2,cartesian),"description")),
                      Contains("Could not find .non existent 2d Point, while it is set as required."));
    CHECK_THROWS_WITH(prm.load_entry("non existent 3d Point", true, Types::Point<3>(Point<3>(1,2,3,cartesian),"description")),
                      Contains("Could not find .non existent 3d Point, while it is set as required."));

  #ifndef NDEBUG
    CHECK_THROWS_WITH(prm.get_point<2>("non existent 2d Point"),
                      Contains("Could not find entry 'non existent 2d Point' not found. Make sure it is loaded or set"));
    CHECK_THROWS_WITH(prm.get_point<3>("non existent 3d Point"),
                      Contains("Could not find entry 'non existent 3d Point' not found. Make sure it is loaded or set"));
  #endif

    CHECK(prm.load_entry("non existent 2d Point", false, Types::Point<2>(Point<2>(1,2,cartesian),"description")) == false);
    CHECK(prm.load_entry("non existent 3d Point", false, Types::Point<3>(Point<3>(1,2,3,cartesian),"description")) == false);

    CHECK(prm.get_point<2>("non existent 2d Point").get_array() == std::array<double,2> {1,2});
    CHECK(prm.get_point<3>("non existent 3d Point").get_array() == std::array<double,3> {1,2,3});

    prm.set_entry("new Point 2d", Types::Point<2>(Point<2>(3,4,cartesian),"description"));
    prm.set_entry("new Point 3d", Types::Point<3>(Point<3>(5,6,7,cartesian),"description"));

    CHECK(prm.get_point<2>("new Point 2d").get_array() == std::array<double,2> {3,4});
    CHECK(prm.get_point<3>("new Point 3d").get_array() == std::array<double,3> {5,6,7});


    prm.load_entry("2d point", true, Types::Point<2>(Point<2>(1,2,cartesian),"description"));
    prm.load_entry("3d point", true, Types::Point<3>(Point<3>(3,4,5,cartesian),"description"));

    CHECK(prm.get_point<2>("2d point").get_array() == std::array<double,2> {10,11});
    CHECK(prm.get_point<3>("3d point").get_array() == std::array<double,3> {12,13,14});

    // Test the Array<Types::Double> functions
    CHECK_THROWS_WITH(prm.load_entry("non existent double array", true, Types::Array(Types::Double(1,"description"),"description")),
                      Contains("Could not find .non existent double array, while it is set as required."));

  #ifndef NDEBUG
    CHECK_THROWS_WITH(prm.get_array("non existent double array"),
                      Contains("Could not find entry 'non existent double array' not found. Make sure it is loaded or set"));
  #endif

    CHECK(prm.load_entry("non exitent double array", false, Types::Array(Types::Double(2,"description"),"description")) == false);

  #ifndef NDEBUG
    CHECK(prm.get_array<Types::Double>("non exitent double array").size() == 1);
    CHECK(prm.get_array<Types::Double>("non exitent double array")[0].value == 2.0);
    CHECK(prm.get_array<Types::Double>("non exitent double array")[0].description == "description");
  #endif
    // This is not desired behavior, but it is not implemented yet.

    prm.set_entry("new double array", Types::Array(Types::Double(3,"description"),"description"));
    std::vector<Types::Double> set_typed_double =  prm.get_array<Types::Double >("new double array");
    CHECK(set_typed_double.size() == 0);
    // This is not desired behavior, but it is not implemented yet.

    prm.load_entry("double array", true, Types::Array(Types::Double(4,"description"),"description"));
    std::vector<Types::Double> true_loaded_typed_double =  prm.get_array<Types::Double >("double array");
    CHECK(true_loaded_typed_double.size() == 3);
    CHECK(true_loaded_typed_double[0].value == 25);
    CHECK(true_loaded_typed_double[1].value == 26);
    CHECK(true_loaded_typed_double[2].value == 27);


    // Test the Array<Types::Point<2> > functions
    CHECK_THROWS_WITH(prm.load_entry("non existent point<2> array", true, Types::Array(Types::Point<2>(Point<2>(1,2,cartesian),"description"),"description")),
                      Contains("Could not find .non existent point<2> array, while it is set as required."));

  #ifndef NDEBUG
    CHECK_THROWS_WITH(prm.get_array("non existent point<2> array"),
                      Contains("Could not find entry 'non existent point<2> array' not found. Make sure it is loaded or set"));
  #endif
    CHECK(prm.load_entry("non exitent double array", false, Types::Array(Types::Point<2>(Point<2>(3,4,cartesian),"description"),"description")) == false);

  #ifndef NDEBUG
    CHECK_THROWS_WITH(prm.get_array<Types::Point<2> >("non existent point<2> array"),
                      Contains("Could not find entry 'non existent point<2> array' not found. Make sure it is loaded or set."));
    // This is not desired behavior, but it is not implemented yet.
  #endif

    prm.set_entry("new point<2> array", Types::Array(Types::Point<2>(Point<2>(5,6,cartesian),"description"),"description"));
    std::vector<Types::Point<2> > set_typed_point_2d = prm.get_array<Types::Point<2> >("new point<2> array");
    CHECK(set_typed_point_2d.size() == 0);
    // This is not desired behavior, but it is not implemented yet.

    prm.load_entry("point<2> array", true, Types::Array(Types::Point<2>(Point<2>(7,8,cartesian),"description"),"description"));
    std::vector<Types::Point<2> > true_loaded_typed_point_2d =  prm.get_array<Types::Point<2> >("point<2> array");
    CHECK(true_loaded_typed_point_2d.size() == 3);
    CHECK(true_loaded_typed_point_2d[0].value.get_array() == std::array<double,2> {10,11});
    CHECK(true_loaded_typed_point_2d[1].value.get_array() == std::array<double,2> {12,13});
    CHECK(true_loaded_typed_point_2d[2].value.get_array() == std::array<double,2> {14,15});


    // Test the Array<Types::Point<3> > functions
    CHECK_THROWS_WITH(prm.load_entry("non existent point<3> array", true, Types::Array(Types::Point<3>(Point<3>(1,2,3,cartesian),"description"),"description")),
                      Contains("Could not find .non existent point<3> array, while it is set as required."));

  #ifndef NDEBUG
    CHECK_THROWS_WITH(prm.get_array("non existent point<3> array"),
                      Contains("Could not find entry 'non existent point<3> array' not found. Make sure it is loaded or set"));
  #endif
    CHECK(prm.load_entry("non exitent double array", false, Types::Array(Types::Point<3>(Point<3>(4,5,6,cartesian),"description"),"description")) == false);

  #ifndef NDEBUG
    CHECK_THROWS_WITH(prm.get_array<Types::Point<3> >("non existent point<3> array"),
                      Contains("Could not find entry 'non existent point<3> array' not found. Make sure it is loaded or set."));
    // This is not desired behavior, but it is not implemented yet.
  #endif

    prm.set_entry("new point<3> array", Types::Array(Types::Point<3>(Point<3>(7,8,9,cartesian),"description"),"description"));
    std::vector<Types::Point<3> > set_typed_point_3d = prm.get_array<Types::Point<3> >("new point<3> array");
    CHECK(set_typed_point_3d.size() == 0);
    // This is not desired behavior, but it is not implemented yet.

    prm.load_entry("point<3> array", true, Types::Array(Types::Point<3>(Point<3>(10,11,12,cartesian),"description"),"description"));
    std::vector<Types::Point<3> > true_loaded_typed_point_3d =  prm.get_array<Types::Point<3> >("point<3> array");
    CHECK(true_loaded_typed_point_3d.size() == 3);
    CHECK(true_loaded_typed_point_3d[0].value.get_array() == std::array<double,3> {20,21,22});
    CHECK(true_loaded_typed_point_3d[1].value.get_array() == std::array<double,3> {23,24,25});
    CHECK(true_loaded_typed_point_3d[2].value.get_array() == std::array<double,3> {26,27,28});

  #ifndef NDEBUG
    CHECK_THROWS_WITH(prm.get_array<Types::Double >("point<2> array"),
                      Contains("Could not get point<2> array, because it is not a 2d Point."));
    CHECK_THROWS_WITH(prm.get_array<Types::Double >("point<3> array"),
                      Contains("Could not get point<3> array, because it is not a 3d Point."));

    CHECK_THROWS_WITH(prm.get_array<Types::Point<2> >("point<3> array"),
                      Contains("Could not get point<3> array, because it is not a 3d Point."));
    CHECK_THROWS_WITH(prm.get_array<Types::Point<2> >("double array"),
                      Contains("Could not get double array, because it is not a Double."));

    CHECK_THROWS_WITH(prm.get_array<Types::Point<3> >("point<2> array"),
                      Contains("Could not get point<2> array, because it is not a 2d Point."));
    CHECK_THROWS_WITH(prm.get_array<Types::Point<3> >("double array"),
                      Contains("Could not get double array, because it is not a Double."));
  #endif

    // Test the enter_subsection and leave_subsection functions
    prm.enter_subsection("subsection 1");
    {
      // Test the UnsignedInt functions
      / *CHECK_THROWS_WITH(prm.load_entry("non existent unsigned int", true, Types::UnsignedInt(1,"description")),
                        Contains("Entry undeclared: subsection 1.non existent unsigned int"));

      #ifndef NDEBUG
      CHECK_THROWS_WITH(prm.get_unsigned_int("non existent unsigned int"),
                        Contains("Could not find entry 'non existent unsigned int' not found. Make sure it is loaded or set"));
      #endif

      CHECK(prm.load_entry("non existent unsigned int", false, Types::UnsignedInt(1,"description")) == false);
      CHECK(prm.get_unsigned_int("non existent unsigned int") == 1);

      prm.set_entry("new unsigned int", Types::UnsignedInt(2,"description"));
      CHECK(prm.get_unsigned_int("new unsigned int") == 2);

      prm.load_entry("unsigned int", true, Types::UnsignedInt(3,"description"));
      CHECK(prm.get_unsigned_int("unsigned int") == 5);* /


      // Test the Double functions
      CHECK_THROWS_WITH(prm.load_entry("non existent double", true, Types::Double(1,"description")),
                        Contains("Entry undeclared: subsection 1.non existent"));

  #ifndef NDEBUG
      CHECK_THROWS_WITH(prm.get_double("non existent double"),
                        Contains("Could not find entry 'non existent double' not found. Make sure it is loaded or set"));
  #endif

      CHECK(prm.load_entry("non existent double", false, Types::Double(2,"description")) == false);
      CHECK(prm.get_double("non existent double") == 2);

      prm.set_entry("new double", Types::Double(3,"description"));
      CHECK(prm.get_double("new double") == 3);

      prm.load_entry("double", true, Types::Double(4,"description"));
      CHECK(prm.get_double("double") == 2.23456e2);


      // Test the String functions
      CHECK_THROWS_WITH(prm.load_entry("non existent string", true, Types::String("2","description")),
                        Contains("Entry undeclared: subsection 1.non existent string"));

  #ifndef NDEBUG
      CHECK_THROWS_WITH(prm.get_string("non existent string"),
                        Contains("Could not find entry 'non existent string' not found. Make sure it is loaded or set"));
  #endif

      CHECK(prm.load_entry("non exitent string", false, Types::String("3","description")) == false);
      CHECK(prm.get_string("non exitent string") == "3");

      prm.set_entry("new string", Types::String("4","description"));
      CHECK(prm.get_string("new string") == "4");

      prm.load_entry("string", true, Types::String("5","description"));
      CHECK(prm.get_string("string") == "mystring 1");

      // Test the Point functions
      CHECK_THROWS_WITH(prm.load_entry("non existent 2d Point", true, Types::Point<2>(Point<2>(3,4,cartesian),"description")),
                        Contains("Could not find subsection 1.non existent 2d Point, while it is set as required."));
      CHECK_THROWS_WITH(prm.load_entry("non existent 3d Point", true, Types::Point<3>(Point<3>(4,5,6,cartesian),"description")),
                        Contains("Could not find subsection 1.non existent 3d Point, while it is set as required."));

  #ifndef NDEBUG
      CHECK_THROWS_WITH(prm.get_point<2>("non existent 2d Point"),
                        Contains("Could not find entry 'non existent 2d Point' not found. Make sure it is loaded or set"));
      CHECK_THROWS_WITH(prm.get_point<3>("non existent 3d Point"),
                        Contains("Could not find entry 'non existent 3d Point' not found. Make sure it is loaded or set"));
  #endif

      CHECK(prm.load_entry("non existent 2d Point", false, Types::Point<2>(Point<2>(3,4,cartesian),"description")) == false);
      CHECK(prm.load_entry("non existent 3d Point", false, Types::Point<3>(Point<3>(4,5,6,cartesian),"description")) == false);

      CHECK(prm.get_point<2>("non existent 2d Point").get_array() == std::array<double,2> {3,4});
      CHECK(prm.get_point<3>("non existent 3d Point").get_array() == std::array<double,3> {4,5,6});

      prm.set_entry("new Point 2d", Types::Point<2>(Point<2>(5,6,cartesian),"description"));
      prm.set_entry("new Point 3d", Types::Point<3>(Point<3>(7,8,9,cartesian),"description"));

      CHECK(prm.get_point<2>("new Point 2d").get_array() == std::array<double,2> {5,6});
      CHECK(prm.get_point<3>("new Point 3d").get_array() == std::array<double,3> {7,8,9});


      prm.load_entry("2d point", true, Types::Point<2>(Point<2>(1,2,cartesian),"description"));
      prm.load_entry("3d point", true, Types::Point<3>(Point<3>(3,4,5,cartesian),"description"));

      CHECK(prm.get_point<2>("2d point").get_array() == std::array<double,2> {15,16});
      CHECK(prm.get_point<3>("3d point").get_array() == std::array<double,3> {17,18,19});


      // Test the Array functions
      CHECK_THROWS_WITH(prm.load_entry("non existent double array", true, Types::Array(Types::Double(1,"description"),"description")),
                        Contains("Could not find subsection 1.non existent double array, while it is set as required."));

  #ifndef NDEBUG
      CHECK_THROWS_WITH(prm.get_array("non existent double array"),
                        Contains("Could not find entry 'non existent double array' not found. Make sure it is loaded or set"));

      CHECK(prm.load_entry("non exitent double array", false, Types::Array(Types::Double(2,"description"),"description")) == false);
      CHECK(prm.get_array<Types::Double>("non exitent double array").size() == 1);
      CHECK(prm.get_array<Types::Double>("non exitent double array")[0].value == 2.0);
      CHECK(prm.get_array<Types::Double>("non exitent double array")[0].description == "description");
      // This is not desired behavior, but it is not implemented yet.
  #endif

      prm.set_entry("new double array", Types::Array(Types::Double(3,"description"),"description"));
      std::vector<Types::Double > set_typed_double =  prm.get_array<Types::Double >("new double array");
      CHECK(set_typed_double.size() == 0);
      // This is not desired behavior, but it is not implemented yet.

      prm.load_entry("double array", true, Types::Array(Types::Double(4,"description"),"description"));
      std::vector<Types::Double > true_loaded_typed_double =  prm.get_array<Types::Double >("double array");
      CHECK(true_loaded_typed_double.size() == 3);
      CHECK(true_loaded_typed_double[0].value == 35);
      CHECK(true_loaded_typed_double[1].value == 36);
      CHECK(true_loaded_typed_double[2].value == 37);

      // Test the Array<Types::Point<2> > functions
      CHECK_THROWS_WITH(prm.load_entry("non existent point<2> array", true, Types::Array(Types::Point<2>(Point<2>(1,2,cartesian),"description"),"description")),
                        Contains("Could not find subsection 1.non existent point<2> array, while it is set as required."));

  #ifndef NDEBUG
      CHECK_THROWS_WITH(prm.get_array("non existent point<2> array"),
                        Contains("Could not find entry 'non existent point<2> array' not found. Make sure it is loaded or set"));

      CHECK(prm.load_entry("non exitent double array", false, Types::Array(Types::Point<2>(Point<2>(3,4,cartesian),"description"),"description")) == false);
      CHECK_THROWS_WITH(prm.get_array<Types::Point<2> >("non existent point<2> array"),
                        Contains("Could not find entry 'non existent point<2> array' not found. Make sure it is loaded or set."));
      // This is not desired behavior, but it is not implemented yet.
  #endif

      prm.set_entry("new point<2> array", Types::Array(Types::Point<2>(Point<2>(5,6,cartesian),"description"),"description"));
      std::vector<Types::Point<2> > set_typed_point_2d = prm.get_array<Types::Point<2> >("new point<2> array");
      CHECK(set_typed_point_2d.size() == 0);
      // This is not desired behavior, but it is not implemented yet.

      prm.load_entry("point<2> array", true, Types::Array(Types::Point<2>(Point<2>(7,8,cartesian),"description"),"description"));
      std::vector<Types::Point<2> > true_loaded_typed_point_2d =  prm.get_array<Types::Point<2> >("point<2> array");
      CHECK(true_loaded_typed_point_2d.size() == 3);
      CHECK(true_loaded_typed_point_2d[0].value.get_array() == std::array<double,2> {20,21});
      CHECK(true_loaded_typed_point_2d[1].value.get_array() == std::array<double,2> {22,23});
      CHECK(true_loaded_typed_point_2d[2].value.get_array() == std::array<double,2> {24,25});


      // Test the Array<Types::Point<3> > functions
      CHECK_THROWS_WITH(prm.load_entry("non existent point<3> array", true, Types::Array(Types::Point<3>(Point<3>(1,2,3,cartesian),"description"),"description")),
                        Contains("Could not find subsection 1.non existent point<3> array, while it is set as required."));

  #ifndef NDEBUG
      CHECK_THROWS_WITH(prm.get_array("non existent point<3> array"),
                        Contains("Could not find entry 'non existent point<3> array' not found. Make sure it is loaded or set"));

      CHECK(prm.load_entry("non exitent double array", false, Types::Array(Types::Point<3>(Point<3>(4,5,6,cartesian),"description"),"description")) == false);
      CHECK_THROWS_WITH(prm.get_array<Types::Point<3> >("non existent point<3> array"),
                        Contains("Could not find entry 'non existent point<3> array' not found. Make sure it is loaded or set."));
      // This is not desired behavior, but it is not implemented yet.
  #endif

      prm.set_entry("new point<3> array", Types::Array(Types::Point<3>(Point<3>(7,8,9,cartesian),"description"),"description"));
      std::vector<Types::Point<3> > set_typed_point_3d = prm.get_array<Types::Point<3> >("new point<3> array");
      CHECK(set_typed_point_3d.size() == 0);
      // This is not desired behavior, but it is not implemented yet.

      prm.load_entry("point<3> array", true, Types::Array(Types::Point<3>(Point<3>(10,11,12,cartesian),"description"),"description"));
      std::vector<Types::Point<3> > true_loaded_typed_point_3d =  prm.get_array<Types::Point<3> >("point<3> array");
      CHECK(true_loaded_typed_point_3d.size() == 3);
      CHECK(true_loaded_typed_point_3d[0].value.get_array() == std::array<double,3> {30,31,32});
      CHECK(true_loaded_typed_point_3d[1].value.get_array() == std::array<double,3> {33,34,35});
      CHECK(true_loaded_typed_point_3d[2].value.get_array() == std::array<double,3> {36,37,38});

  #ifndef NDEBUG
      CHECK_THROWS_WITH(prm.get_array<Types::Double >("point<2> array"),
                        Contains("Could not get subsection 1.point<2> array, because it is not a 2d Point."));
      CHECK_THROWS_WITH(prm.get_array<Types::Double >("point<3> array"),
                        Contains("Could not get subsection 1.point<3> array, because it is not a 3d Point."));

      CHECK_THROWS_WITH(prm.get_array<Types::Point<2> >("point<3> array"),
                        Contains("Could not get subsection 1.point<3> array, because it is not a 3d Point."));
      CHECK_THROWS_WITH(prm.get_array<Types::Point<2> >("double array"),
                        Contains("Could not get subsection 1.double array, because it is not a Double."));

      CHECK_THROWS_WITH(prm.get_array<Types::Point<3> >("point<2> array"),
                        Contains("Could not get subsection 1.point<2> array, because it is not a 2d Point."));
      CHECK_THROWS_WITH(prm.get_array<Types::Point<3> >("double array"),
                        Contains("Could not get subsection 1.double array, because it is not a Double."));
  #endif

      prm.enter_subsection("subsection 2");
      {
        // Test the UnsignedInt functions
        / *CHECK_THROWS_WITH(prm.load_entry("non existent unsigned int", true, Types::UnsignedInt(1,"description")),
                          Contains("Entry undeclared: subsection 1.subsection 2.non existent unsigned int"));

        #ifndef NDEBUG
        CHECK_THROWS_WITH(prm.get_unsigned_int("non existent unsigned int"),
                          Contains("Could not find entry 'non existent unsigned int' not found. Make sure it is loaded or set"));
        #endif

        CHECK(prm.load_entry("non existent unsigned int", false, Types::UnsignedInt(1,"description")) == false);
        CHECK(prm.get_unsigned_int("non existent unsigned int") == 1);

        prm.set_entry("new unsigned int", Types::UnsignedInt(2,"description"));
        CHECK(prm.get_unsigned_int("new unsigned int") == 2);

        prm.load_entry("unsigned int", true, Types::UnsignedInt(3,"description"));
        CHECK(prm.get_unsigned_int("unsigned int") == 6);* /


        // Test the Double functions
        CHECK_THROWS_WITH(prm.load_entry("non existent double", true, Types::Double(3,"description")),
                          Contains("Entry undeclared: subsection 1.subsection 2.non existent"));

  #ifndef NDEBUG
        CHECK_THROWS_WITH(prm.get_double("non existent double"),
                          Contains("Could not find entry 'non existent double' not found. Make sure it is loaded or set"));
  #endif

        CHECK(prm.load_entry("non existent double", false, Types::Double(4,"description")) == false);
        CHECK(prm.get_double("non existent double") == 4);

        prm.set_entry("new double", Types::Double(5,"description"));
        CHECK(prm.get_double("new double") == 5);

        prm.load_entry("double", true, Types::Double(6,"description"));
        CHECK(prm.get_double("double") == 3.23456e2);


        // Test the String functions
        CHECK_THROWS_WITH(prm.load_entry("non existent string", true, Types::String("3","description")),
                          Contains("Entry undeclared: subsection 1.subsection 2.non existent string"));

  #ifndef NDEBUG
        CHECK_THROWS_WITH(prm.get_string("non existent string"),
                          Contains("Could not find entry 'non existent string' not found. Make sure it is loaded or set"));
  #endif

        CHECK(prm.load_entry("non exitent string", false, Types::String("4","description")) == false);
        CHECK(prm.get_string("non exitent string") == "4");

        prm.set_entry("new string", Types::String("5","description"));
        CHECK(prm.get_string("new string") == "5");

        prm.load_entry("string", true, Types::String("6","description"));
        CHECK(prm.get_string("string") == "mystring 2");

        // Test the Point functions

        CHECK_THROWS_WITH(prm.load_entry("non existent 2d Point", true, Types::Point<2>(Point<2>(1,2,cartesian),"description")),
                          Contains("Could not find subsection 1.subsection 2.non existent 2d Point, while it is set as required."));
        CHECK_THROWS_WITH(prm.load_entry("non existent 3d Point", true, Types::Point<3>(Point<3>(1,2,3,cartesian),"description")),
                          Contains("Could not find subsection 1.subsection 2.non existent 3d Point, while it is set as required."));

  #ifndef NDEBUG
        CHECK_THROWS_WITH(prm.get_point<2>("non existent 2d Point"),
                          Contains("Could not find entry 'non existent 2d Point' not found. Make sure it is loaded or set"));
        CHECK_THROWS_WITH(prm.get_point<3>("non existent 3d Point"),
                          Contains("Could not find entry 'non existent 3d Point' not found. Make sure it is loaded or set"));
  #endif

        CHECK(prm.load_entry("non existent 2d Point", false, Types::Point<2>(Point<2>(1,2,cartesian),"description")) == false);
        CHECK(prm.load_entry("non existent 3d Point", false, Types::Point<3>(Point<3>(1,2,3,cartesian),"description")) == false);

        CHECK(prm.get_point<2>("non existent 2d Point").get_array() == std::array<double,2> {1,2});
        CHECK(prm.get_point<3>("non existent 3d Point").get_array() == std::array<double,3> {1,2,3});

        prm.set_entry("new Point 2d", Types::Point<2>(Point<2>(3,4,cartesian),"description"));
        prm.set_entry("new Point 3d", Types::Point<3>(Point<3>(5,6,7,cartesian),"description"));

        CHECK(prm.get_point<2>("new Point 2d").get_array() == std::array<double,2> {3,4});
        CHECK(prm.get_point<3>("new Point 3d").get_array() == std::array<double,3> {5,6,7});


        prm.load_entry("2d point", true, Types::Point<2>(Point<2>(1,2,cartesian),"description"));
        prm.load_entry("3d point", true, Types::Point<3>(Point<3>(3,4,5,cartesian),"description"));

        CHECK(prm.get_point<2>("2d point").get_array() == std::array<double,2> {20,21});
        CHECK(prm.get_point<3>("3d point").get_array() == std::array<double,3> {22,23,24});

        // Test the Array functions

        CHECK_THROWS_WITH(prm.load_entry("non existent double array", true, Types::Array(Types::Double(1,"description"),"description")),
                          Contains("Could not find subsection 1.subsection 2.non existent double array, while it is set as required."));
  #ifndef NDEBUG
        CHECK_THROWS_WITH(prm.get_array("non existent double array"),
                          Contains("Could not find entry 'non existent double array' not found. Make sure it is loaded or set"));
  #endif

        CHECK(prm.load_entry("non exitent double array", false, Types::Array(Types::Double(2,"description"),"description")) == false);

        CHECK(prm.get_array<Types::Double >("non exitent double array").size() == 1);
        CHECK(prm.get_array<Types::Double >("non exitent double array")[0].value == 2.0);
        CHECK(prm.get_array<Types::Double >("non exitent double array")[0].description == "description");

        prm.set_entry("new double array", Types::Array(Types::Double(3,"description"),"description"));
        std::vector<Types::Double > set_typed_double =  prm.get_array<Types::Double >("new double array");
        CHECK(set_typed_double.size() == 0);
        // This is not desired behavior, but it is not implemented yet.

        prm.load_entry("double array", true, Types::Array(Types::Double(4,"description"),"description"));
        std::vector<Types::Double > true_loaded_typed_double =  prm.get_array<Types::Double >("double array");
        CHECK(true_loaded_typed_double.size() == 3);
        CHECK(true_loaded_typed_double[0].value == 45);
        CHECK(true_loaded_typed_double[1].value == 46);
        CHECK(true_loaded_typed_double[2].value == 47);


        // Test the Array<Types::Point<2> > functions
  #ifndef NDEBUG
        CHECK_THROWS_WITH(prm.load_entry("non existent point<2> array", true, Types::Array(Types::Point<2>(Point<2>(1,2,cartesian),"description"),"description")),
                          Contains("Could not find subsection 1.subsection 2.non existent point<2> array, while it is set as required."));

        CHECK_THROWS_WITH(prm.get_array("non existent point<2> array"),
                          Contains("Could not find entry 'non existent point<2> array' not found. Make sure it is loaded or set"));
  #endif

        CHECK(prm.load_entry("non exitent double array", false, Types::Array(Types::Point<2>(Point<2>(3,4,cartesian),"description"),"description")) == false);

  #ifndef NDEBUG
        CHECK_THROWS_WITH(prm.get_array<Types::Point<2> >("non existent point<2> array"),
                          Contains("Could not find entry 'non existent point<2> array' not found. Make sure it is loaded or set."));
        // This is not desired behavior, but it is not implemented yet.
  #endif

        prm.set_entry("new point<2> array", Types::Array(Types::Point<2>(Point<2>(5,6,cartesian),"description"),"description"));
        std::vector<Types::Point<2> > set_typed_point_2d = prm.get_array<Types::Point<2> >("new point<2> array");
        CHECK(set_typed_point_2d.size() == 0);
        // This is not desired behavior, but it is not implemented yet.

        prm.load_entry("point<2> array", true, Types::Array(Types::Point<2>(Point<2>(7,8,cartesian),"description"),"description"));
        std::vector<Types::Point<2> > true_loaded_typed_point_2d =  prm.get_array<Types::Point<2> >("point<2> array");
        CHECK(true_loaded_typed_point_2d.size() == 3);
        CHECK(true_loaded_typed_point_2d[0].value.get_array() == std::array<double,2> {40,41});
        CHECK(true_loaded_typed_point_2d[1].value.get_array() == std::array<double,2> {42,43});
        CHECK(true_loaded_typed_point_2d[2].value.get_array() == std::array<double,2> {44,45});


        // Test the Array<Types::Point<3> > functions
  #ifndef NDEBUG
        CHECK_THROWS_WITH(prm.load_entry("non existent point<3> array", true, Types::Array(Types::Point<3>(Point<3>(1,2,3,cartesian),"description"),"description")),
                          Contains("Could not find subsection 1.subsection 2.non existent point<3> array, while it is set as required."));

        CHECK_THROWS_WITH(prm.get_array("non existent point<3> array"),
                          Contains("Could not find entry 'non existent point<3> array' not found. Make sure it is loaded or set"));
  #endif

        CHECK(prm.load_entry("non exitent double array", false, Types::Array(Types::Point<3>(Point<3>(4,5,6,cartesian),"description"),"description")) == false);

  #ifndef NDEBUG
        CHECK_THROWS_WITH(prm.get_array<Types::Point<3> >("non existent point<3> array"),
                          Contains("Could not find entry 'non existent point<3> array' not found. Make sure it is loaded or set."));
        // This is not desired behavior, but it is not implemented yet.
  #endif

        prm.set_entry("new point<3> array", Types::Array(Types::Point<3>(Point<3>(7,8,9,cartesian),"description"),"description"));
        std::vector<Types::Point<3> > set_typed_point_3d = prm.get_array<Types::Point<3> >("new point<3> array");
        CHECK(set_typed_point_3d.size() == 0);
        // This is not desired behavior, but it is not implemented yet.

        prm.load_entry("point<3> array", true, Types::Array(Types::Point<3>(Point<3>(10,11,12,cartesian),"description"),"description"));
        std::vector<Types::Point<3> > true_loaded_typed_point_3d =  prm.get_array<Types::Point<3> >("point<3> array");
        CHECK(true_loaded_typed_point_3d.size() == 3);
        CHECK(true_loaded_typed_point_3d[0].value.get_array() == std::array<double,3> {40,41,42});
        CHECK(true_loaded_typed_point_3d[1].value.get_array() == std::array<double,3> {43,44,45});
        CHECK(true_loaded_typed_point_3d[2].value.get_array() == std::array<double,3> {46,47,48});

  #ifndef NDEBUG
        CHECK_THROWS_WITH(prm.get_array<Types::Double >("point<2> array"),
                          Contains("Could not get subsection 1.subsection 2.point<2> array, because it is not a 2d Point."));
        CHECK_THROWS_WITH(prm.get_array<Types::Double >("point<3> array"),
                          Contains("Could not get subsection 1.subsection 2.point<3> array, because it is not a 3d Point."));

        CHECK_THROWS_WITH(prm.get_array<Types::Point<2> >("point<3> array"),
                          Contains("Could not get subsection 1.subsection 2.point<3> array, because it is not a 3d Point."));
        CHECK_THROWS_WITH(prm.get_array<Types::Point<2> >("double array"),
                          Contains("Could not get subsection 1.subsection 2.double array, because it is not a Double."));

        CHECK_THROWS_WITH(prm.get_array<Types::Point<3> >("point<2> array"),
                          Contains("Could not get subsection 1.subsection 2.point<2> array, because it is not a 2d Point."));
        CHECK_THROWS_WITH(prm.get_array<Types::Point<3> >("double array"),
                          Contains("Could not get subsection 1.subsection 2.double array, because it is not a Double."));
  #endif
      }
      prm.leave_subsection();
    }
    prm.leave_subsection();
  */

  // Todo: add tests for list,feature and coordinate system.




}

TEST_CASE("WorldBuilder Utilities function: distance_point_from_curved_planes cartesian part 1")
{
  std::unique_ptr<CoordinateSystems::Interface> cartesian_system = CoordinateSystems::Interface::create("cartesian", NULL);;

  //Todo:fix
  //cartesian_system->declare_entries();

  Point<3> position(10,0,0,cartesian);
  Point<2> reference_point(0,0,cartesian);

  std::vector<Point<2> > coordinates;
  coordinates.push_back(Point<2>(0,10,cartesian));
  coordinates.push_back(Point<2>(20,10,cartesian));

  std::vector<std::vector<double> > slab_segment_lengths(2);
  slab_segment_lengths[0].push_back(std::sqrt(10*10+10*10));
  slab_segment_lengths[0].push_back(200);
  slab_segment_lengths[1].push_back(std::sqrt(10*10+10*10));
  slab_segment_lengths[1].push_back(200);

  double dtr = Utilities::const_pi/180;
  std::vector<std::vector<Point<2> > > slab_segment_angles(2);
  slab_segment_angles[0].push_back(Point<2>(45 * dtr,45 * dtr,cartesian));
  slab_segment_angles[0].push_back(Point<2>(45 * dtr,45 * dtr,cartesian));
  slab_segment_angles[1].push_back(Point<2>(45 * dtr,45 * dtr,cartesian));
  slab_segment_angles[1].push_back(Point<2>(45 * dtr,45 * dtr,cartesian));

  double starting_radius = 10;

  std::map<std::string,double> distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(distance_from_planes["distanceFromPlane"] == Approx(-3.97205e-15)); // practically zero
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(std::sqrt(10*10+10*10)));
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 0);
  CHECK(distance_from_planes["segmentFraction"] == Approx(1.0));


  // center square test 2
  reference_point[1] = 20;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(distance_from_planes["distanceFromPlane"] == Approx(std::sqrt(10*10+10*10)));
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(7.10543e-16)); // practically zero
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 0);
  CHECK(distance_from_planes["segmentFraction"] == Approx(5.0243e-17)); // practically zero

  // center square test 3
  position[1] = 20;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(distance_from_planes["distanceFromPlane"] == Approx(-3.97205e-15)); // practically zero
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(std::sqrt(10*10+10*10)));
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 0);
  CHECK(distance_from_planes["segmentFraction"] == Approx(1.0));

  // center square test 4
  reference_point[1] = 0;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(distance_from_planes["distanceFromPlane"] == Approx(std::sqrt(10*10+10*10)));
  CHECK(std::fabs(distance_from_planes["distanceAlongPlane"]) < 1e-14); // practically zero
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 0);
  CHECK(std::fabs(distance_from_planes["segmentFraction"]) < 1e-14); // practically zero

  // center square test 5
  position[1] = -10;
  position[2] = -10;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(std::fabs(distance_from_planes["distanceFromPlane"]) < 1e-14);
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(sqrt(20*20+20*20))); // practically zero
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 1);
  CHECK(distance_from_planes["segmentFraction"] == Approx(0.0707106781)); // practically zero

  // begin section square test 6
  position[0] = 0;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(std::fabs(distance_from_planes["distanceFromPlane"]) < 1e-14);
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(sqrt(20*20+20*20))); // practically zero
  CHECK(std::fabs(distance_from_planes["sectionFraction"]) < 1e-14);
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 1);
  CHECK(distance_from_planes["segmentFraction"] == Approx(0.0707106781)); // practically zero


  // end section square test 7
  position[0] = 20;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(std::fabs(distance_from_planes["distanceFromPlane"]) < 1e-14);
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(sqrt(20*20+20*20))); // practically zero
  CHECK(distance_from_planes["sectionFraction"] == Approx(1.0));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 1);
  CHECK(distance_from_planes["segmentFraction"] == Approx(0.0707106781)); // practically zero

  // before begin section square test 8
  position[0] = -10;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(distance_from_planes["distanceFromPlane"] == INFINITY);
  CHECK(distance_from_planes["distanceAlongPlane"] == INFINITY); // practically zero
  CHECK(distance_from_planes["sectionFraction"] == 0.0);
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 0);
  CHECK(distance_from_planes["segmentFraction"] == 0.0); // practically zero

  // beyond end section square test 9
  position[0] = 25;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(distance_from_planes["distanceFromPlane"] == INFINITY);
  CHECK(distance_from_planes["distanceAlongPlane"] == INFINITY); // practically zero
  CHECK(distance_from_planes["sectionFraction"] == 0.0);
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 0);
  CHECK(distance_from_planes["segmentFraction"] == 0.0); // practically zero


  // beyond end section square test 10
  position[0] = 10;
  position[1] = 0;
  position[2] = 5;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(distance_from_planes["distanceFromPlane"] == Approx(-3.5355339059));
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(10.6066017178));
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 0);
  CHECK(distance_from_planes["segmentFraction"] == Approx(0.75));

  // beyond end section square test 10 (only positive version)
  position[0] = 10;
  position[1] = 0;
  position[2] = 5;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 true);

  CHECK(distance_from_planes["distanceFromPlane"] == Approx(3.5355339059));
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(10.6066017178));
  CHECK(distance_from_planes["sectionFraction"] ==  Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 0);
  CHECK(distance_from_planes["segmentFraction"] == Approx(0.75));


  // beyond end section square test 11
  position[0] = 10;
  position[1] = 0;
  position[2] = -5;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(distance_from_planes["distanceFromPlane"] == Approx(3.5355339059));
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(17.6776695297)); // practically zero
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 1);
  CHECK(distance_from_planes["segmentFraction"] == Approx(0.0176776695)); // practically zero


  // beyond end section square test 11 (only positve version)
  position[0] = 10;
  position[1] = 0;
  position[2] = -5;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 true);

  CHECK(distance_from_planes["distanceFromPlane"] == Approx(3.5355339059));
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(17.6776695297)); // practically zero
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 1);
  CHECK(distance_from_planes["segmentFraction"] == Approx(0.0176776695)); // practically zero

  // add coordinate
  position[0] = 25;
  position[1] = 0;
  position[2] = 0;

  coordinates.push_back(Point<2>(30,10,cartesian));

  slab_segment_lengths.resize(3);
  slab_segment_lengths[2].push_back(std::sqrt(10*10+10*10));
  slab_segment_lengths[2].push_back(200);

  slab_segment_angles.resize(3);
  slab_segment_angles[2].push_back(Point<2>(45 * dtr,45 * dtr,cartesian));
  slab_segment_angles[2].push_back(Point<2>(45 * dtr,45 * dtr,cartesian));

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(distance_from_planes["distanceFromPlane"] == Approx(-3.97205e-15)); // practically zero
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(std::sqrt(10*10+10*10)));
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 1);
  CHECK(distance_from_planes["segment"] == 0);
  CHECK(distance_from_planes["segmentFraction"] == Approx(1.0));

  // different angle
  slab_segment_angles[0][0][0] = 22.5 * dtr;
  slab_segment_angles[0][0][1] = 22.5 * dtr;
  slab_segment_angles[0][1][0] = 22.5 * dtr;
  slab_segment_angles[0][1][1] = 22.5 * dtr;
  slab_segment_angles[1][0][0] = 22.5 * dtr;
  slab_segment_angles[1][0][1] = 22.5 * dtr;
  slab_segment_angles[1][1][0] = 22.5 * dtr;
  slab_segment_angles[1][1][1] = 22.5 * dtr;

  position[0] = 10;
  position[1] = 0;
  position[2] = 10-10*tan(22.5*dtr);

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(std::fabs(distance_from_planes["distanceFromPlane"]) < 1e-14); // practically zero
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(10.8239219938));
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 0);
  CHECK(distance_from_planes["segmentFraction"] == Approx(0.7653668647));

  // check interpolation 1 (in the middle of a segment with 22.5 degree and a segement with 45)
  position[0] = 25;
  position[1] = 0;
  position[2] = 10-10*tan((22.5*1.5)*dtr);

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(std::fabs(distance_from_planes["distanceFromPlane"]) < 1e-14);
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(12.0268977387)); // practically zero
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 1);
  CHECK(distance_from_planes["segment"] == 0);
  CHECK(distance_from_planes["segmentFraction"] == Approx(0.8504300948));

  // check interpolation 2 (at the end of the segment at 45 degree)
  position[0] = 30;
  position[1] = 0;
  position[2] = 10-10*tan(45*dtr);

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(std::fabs(distance_from_planes["distanceFromPlane"]) < 1e-14);
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(14.1421356237)); // practically zero
  CHECK(distance_from_planes["sectionFraction"] == Approx(1.0));
  CHECK(distance_from_planes["section"] == 1);
  CHECK(distance_from_planes["segment"] == 0);
  CHECK(distance_from_planes["segmentFraction"] == Approx(1.0));

  // check length interpolation with 90 degree angles for simplicity
  // check length interpolation first segment center 1
  slab_segment_angles[0][0][0] = 90 * dtr;
  slab_segment_angles[0][0][1] = 90 * dtr;
  slab_segment_angles[0][1][0] = 90 * dtr;
  slab_segment_angles[0][1][1] = 90 * dtr;
  slab_segment_angles[1][0][0] = 90 * dtr;
  slab_segment_angles[1][0][1] = 90 * dtr;
  slab_segment_angles[1][1][0] = 90 * dtr;
  slab_segment_angles[1][1][1] = 90 * dtr;
  slab_segment_angles[2][0][0] = 90 * dtr;
  slab_segment_angles[2][0][1] = 90 * dtr;
  slab_segment_angles[2][1][0] = 90 * dtr;
  slab_segment_angles[2][1][1] = 90 * dtr;

  slab_segment_lengths[0][0] = 100;
  slab_segment_lengths[0][1] = 100;
  slab_segment_lengths[1][0] = 100;
  slab_segment_lengths[1][1] = 100;
  slab_segment_lengths[2][0] = 50;
  slab_segment_lengths[2][1] = 50;

  position[0] = 10;
  position[1] = 10;
  position[2] = 10-100;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(std::fabs(distance_from_planes["distanceFromPlane"]) < 1e-14);
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(100.0)); // practically zero
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 0);
  CHECK(distance_from_planes["segmentFraction"] == Approx(1.0));

  // check length interpolation first segment center 2
  position[0] = 10;
  position[1] = 10;
  position[2] = 10-101;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(std::fabs(distance_from_planes["distanceFromPlane"]) < 1e-14);
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(101.0)); // practically zero
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 1);
  CHECK(distance_from_planes["segmentFraction"] == Approx(0.01));

  // check length interpolation first segment center 3
  position[0] = 10;
  position[1] = 10;
  position[2] = 10-200;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(std::fabs(distance_from_planes["distanceFromPlane"]) < 1e-14);
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(200.0));
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 1);
  CHECK(distance_from_planes["segmentFraction"] == Approx(1.0));



  // check length interpolation first segment center 4
  position[0] = 10;
  position[1] = 10;
  position[2] = 10-201;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(distance_from_planes["distanceFromPlane"] == INFINITY);
  CHECK(distance_from_planes["distanceAlongPlane"] == INFINITY);
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.0));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 0);
  CHECK(distance_from_planes["segmentFraction"] == Approx(0.0));


  // Now check the center of the second segment, each segment should have a length of 75.
  // check length interpolation second segment center 1
  position[0] = 25;
  position[1] = 10;
  position[2] = 10-75;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(std::fabs(distance_from_planes["distanceFromPlane"]) < 1e-14);
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(75.0)); // practically zero
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 1);
  CHECK(distance_from_planes["segment"] == 0);
  CHECK(distance_from_planes["segmentFraction"] == Approx(1.0));

  // check length interpolation second segment center 2
  position[0] = 25;
  position[1] = 10;
  position[2] = 10-76;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(std::fabs(distance_from_planes["distanceFromPlane"]) < 1e-14);
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(76.0)); // practically zero
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 1);
  CHECK(distance_from_planes["segment"] == 1);
  CHECK(distance_from_planes["segmentFraction"] == Approx(0.01333333333333));

  // check length interpolation second segment center 3
  position[0] = 25;
  position[1] = 10;
  position[2] = 10-150;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(std::fabs(distance_from_planes["distanceFromPlane"]) < 1e-14);
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(150.0));
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 1);
  CHECK(distance_from_planes["segment"] == 1);
  CHECK(distance_from_planes["segmentFraction"] == Approx(1.0));



  // check length interpolation second segment center 4
  position[0] = 25;
  position[1] = 10;
  position[2] = 10-151;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(distance_from_planes["distanceFromPlane"] == INFINITY);
  CHECK(distance_from_planes["distanceAlongPlane"] == INFINITY);
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.0));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 0);
  CHECK(distance_from_planes["segmentFraction"] == Approx(0.0));

  // Now check the end of the second segment, each segment should have a length of 50.
  // check length interpolation second segment center 1
  position[0] = 30;
  position[1] = 10;
  position[2] = 10-50;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(std::fabs(distance_from_planes["distanceFromPlane"]) < 1e-14);
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(50.0)); // practically zero
  CHECK(distance_from_planes["sectionFraction"] == Approx(1.0));
  CHECK(distance_from_planes["section"] == 1);
  CHECK(distance_from_planes["segment"] == 0);
  CHECK(distance_from_planes["segmentFraction"] == Approx(1.0));

  // check length interpolation second segment center 2
  position[0] = 30;
  position[1] = 10;
  position[2] = 10-51;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(std::fabs(distance_from_planes["distanceFromPlane"]) < 1e-14);
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(51.0)); // practically zero
  CHECK(distance_from_planes["sectionFraction"] == Approx(1.0));
  CHECK(distance_from_planes["section"] == 1);
  CHECK(distance_from_planes["segment"] == 1);
  CHECK(distance_from_planes["segmentFraction"] == Approx(0.02));

  // check length interpolation second segment center 3
  position[0] = 30;
  position[1] = 10;
  position[2] = 10-100;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(std::fabs(distance_from_planes["distanceFromPlane"]) < 1e-14);
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(100.0));
  CHECK(distance_from_planes["sectionFraction"] == Approx(1.0));
  CHECK(distance_from_planes["section"] == 1);
  CHECK(distance_from_planes["segment"] == 1);
  CHECK(distance_from_planes["segmentFraction"] == Approx(1.0));



  // check length interpolation second segment center 4
  position[0] = 30;
  position[1] = 10;
  position[2] = 10-101;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(distance_from_planes["distanceFromPlane"] == INFINITY);
  CHECK(distance_from_planes["distanceAlongPlane"] == INFINITY);
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.0));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 0);
  CHECK(distance_from_planes["segmentFraction"] == Approx(0.0));
}

TEST_CASE("WorldBuilder Utilities function: distance_point_from_curved_planes cartesian part 2")
{
  std::unique_ptr<CoordinateSystems::Interface> cartesian_system = CoordinateSystems::Interface::create("cartesian", NULL);;

  //todo: fix
  //cartesian_system->declare_entries();

  Point<3> position(10,0,0,cartesian);
  Point<2> reference_point(0,0,cartesian);

  std::vector<Point<2> > coordinates;
  coordinates.push_back(Point<2>(0,10,cartesian));
  coordinates.push_back(Point<2>(20,10,cartesian));
  coordinates.push_back(Point<2>(30,10,cartesian));

  std::vector<std::vector<double> > slab_segment_lengths(3);
  slab_segment_lengths[0].push_back(std::sqrt(10*10+10*10));
  slab_segment_lengths[0].push_back(200);
  slab_segment_lengths[1].push_back(std::sqrt(10*10+10*10));
  slab_segment_lengths[1].push_back(200);
  slab_segment_lengths[2].push_back(std::sqrt(10*10+10*10));
  slab_segment_lengths[2].push_back(200);

  double dtr = Utilities::const_pi/180;
  std::vector<std::vector<Point<2> > > slab_segment_angles(3);
  slab_segment_angles[0].push_back(Point<2>(45 * dtr,45 * dtr,cartesian));
  slab_segment_angles[0].push_back(Point<2>(45 * dtr,45 * dtr,cartesian));
  slab_segment_angles[1].push_back(Point<2>(45 * dtr,45 * dtr,cartesian));
  slab_segment_angles[1].push_back(Point<2>(45 * dtr,45 * dtr,cartesian));
  slab_segment_angles[2].push_back(Point<2>(45 * dtr,45 * dtr,cartesian));
  slab_segment_angles[2].push_back(Point<2>(45 * dtr,45 * dtr,cartesian));

  double starting_radius = 10;
  // Now test the curves into the depth
  // curve test 1

  slab_segment_angles[0][0][0] = 0.0 * dtr;
  slab_segment_angles[0][0][1] = 45.0 * dtr;
  slab_segment_angles[0][1][0] = 45.0 * dtr;
  slab_segment_angles[0][1][1] = 90.0 * dtr;
  slab_segment_angles[1][0][0] = 0.0 * dtr;
  slab_segment_angles[1][0][1] = 45.0 * dtr;
  slab_segment_angles[1][1][0] = 45.0 * dtr;
  slab_segment_angles[1][1][1] = 90.0 * dtr;
  slab_segment_angles[2][0][0] = 90 * dtr;
  slab_segment_angles[2][0][1] = 90 * dtr;
  slab_segment_angles[2][1][0] = 90 * dtr;
  slab_segment_angles[2][1][1] = 90 * dtr;

  slab_segment_lengths[0][0] = 10 * 45 * dtr;
  slab_segment_lengths[0][1] = 10 * 45 * dtr;
  slab_segment_lengths[1][0] = 10 * 45 * dtr;
  slab_segment_lengths[1][1] = 10 * 45 * dtr;
  slab_segment_lengths[2][0] = 5;
  slab_segment_lengths[2][1] = 5;


  position[0] = 10;
  position[1] = 0;
  position[2] = 0;

  std::map<std::string,double> distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(std::fabs(distance_from_planes["distanceFromPlane"]) < 1e-14); // checked that it should be about 0 this with a drawing
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(90.0 * Utilities::const_pi/180 * 10));
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 1);
  CHECK(distance_from_planes["segmentFraction"] == Approx(1.0));

  // curve test 2
  position[0] = 10;
  position[1] = 5;
  position[2] = 0;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(distance_from_planes["distanceFromPlane"] == Approx(5.0)); // checked that it should be about 5 this with a drawing
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(90.0 * Utilities::const_pi/180 * 10));
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 1);
  CHECK(distance_from_planes["segmentFraction"] == Approx(1.0));

  // curve test 3
  position[0] = 10;
  position[1] = -5;
  position[2] = 0;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(distance_from_planes["distanceFromPlane"] == Approx(-5.0)); // checked that it should be about -5 this with a drawing
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(90.0 * Utilities::const_pi/180 * 10));
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 1);
  CHECK(distance_from_planes["segmentFraction"] == Approx(1.0));


  // curve test 4
  position[0] = 10;
  position[1] = 10 - 10 * sqrt(2)/2;
  position[2] = 10 * sqrt(2)/2;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(std::fabs(distance_from_planes["distanceFromPlane"]) < 1e-14); // checked that it should be about 0 this with a drawing
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(45.0 * Utilities::const_pi/180 * 10));
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 0);
  CHECK(distance_from_planes["segmentFraction"] == Approx(1.0));

  // curve test 5
  position[0] = 10;
  position[1] = 10 - 10 * sqrt(2);
  position[2] = 10 * sqrt(2);

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(distance_from_planes["distanceFromPlane"] == Approx(-10.0)); // checked that it should be about -10 this with a drawing
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(45.0 * Utilities::const_pi/180 * 10));
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 0);
  CHECK(distance_from_planes["segmentFraction"] == Approx(1.0));

  // curve test 6
  position[0] = 10;
  position[1] = 10;
  position[2] = 0;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(distance_from_planes["distanceFromPlane"] == Approx(10.0)); // checked that it should be about 10 this with a drawing
  // This is a special case where the point coincides with the center of the circle.
  // Because all the points on the circle are equally close, we have chosen in the
  // code to define this case as that this point belongs to the top of the top segment
  // where the check point has angle 0. This means that the distanceAlongPlate is zero.
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(0.0));
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 0);
  CHECK(distance_from_planes["segmentFraction"] == Approx(0.0));


  // curve test 7
  position[0] = 10;
  position[1] = -5;
  position[2] = -1;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(distance_from_planes["distanceFromPlane"] == INFINITY);
  CHECK(distance_from_planes["distanceAlongPlane"] == INFINITY);
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.0));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 0);
  CHECK(distance_from_planes["segmentFraction"] == Approx(0.0));

  // curve test 8
  slab_segment_lengths[0][0] = 5 * 45 * dtr;
  slab_segment_lengths[0][1] = 5 * 45 * dtr;
  slab_segment_lengths[1][0] = 5 * 45 * dtr;
  slab_segment_lengths[1][1] = 5 * 45 * dtr;

  position[0] = 10;
  position[1] = 5;
  position[2] = 5;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(std::fabs(distance_from_planes["distanceFromPlane"]) < 1e-14); // checked that it should be about 0 this with a drawing
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(90.0 * Utilities::const_pi/180 * 5));
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 1);
  CHECK(distance_from_planes["segmentFraction"] == Approx(1.0));

  // curve test 9
  position[0] = 10;
  position[1] = 10 - 5 * sqrt(2)/2;
  position[2] = 5 + 5 * sqrt(2)/2;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(std::fabs(distance_from_planes["distanceFromPlane"]) < 1e-14); // checked that it should be about 0 this with a drawing
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(45.0 * Utilities::const_pi/180 * 5));
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 0);
  CHECK(distance_from_planes["segmentFraction"] == Approx(1.0));


  // curve test 10
  slab_segment_angles[0][0][0] = 0.0 * dtr;
  slab_segment_angles[0][0][1] = 90.0 * dtr;
  slab_segment_angles[0][1][0] = 90.0 * dtr;
  slab_segment_angles[0][1][1] = 180.0 * dtr;
  slab_segment_angles[1][0][0] = 0.0 * dtr;
  slab_segment_angles[1][0][1] = 90.0 * dtr;
  slab_segment_angles[1][1][0] = 90.0 * dtr;
  slab_segment_angles[1][1][1] = 180.0 * dtr;

  slab_segment_lengths[0][0] = 10 * 90 * dtr;
  slab_segment_lengths[0][1] = 10 * 90 * dtr;
  slab_segment_lengths[1][0] = 10 * 90 * dtr;
  slab_segment_lengths[1][1] = 10 * 90 * dtr;

  position[0] = 10;
  position[1] = 0;
  position[2] = 0;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(std::fabs(distance_from_planes["distanceFromPlane"]) < 1e-14); // checked that it should be about 0 this with a drawing
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(90.0 * Utilities::const_pi/180 * 10));
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 0);
  CHECK(distance_from_planes["segmentFraction"] == Approx(1.0));

  // curve test 11
  position[0] = 10;
  position[1] = 10 - 10 * sqrt(2)/2;
  position[2] = 10 * sqrt(2)/2;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(std::fabs(distance_from_planes["distanceFromPlane"]) < 1e-14); // checked that it should be about 0 this with a drawing
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(45.0 * Utilities::const_pi/180 * 10));
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 0);
  CHECK(distance_from_planes["segmentFraction"] == Approx(0.5));

  // curve test 12
  position[0] = 10;
  position[1] = 10 - 10 * sqrt(2)/2;
  position[2] = -10 * sqrt(2)/2;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(std::fabs(distance_from_planes["distanceFromPlane"]) < 1e-14); // checked that it should be about 0 this with a drawing
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(135.0 * Utilities::const_pi/180 * 10));
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 1);
  CHECK(distance_from_planes["segmentFraction"] == Approx(0.5));


  // curve test 13
  position[0] = 10;
  position[1] = 10;
  position[2] = -10;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(std::fabs(distance_from_planes["distanceFromPlane"]) < 1e-14); // checked that it should be about 0 this with a drawing
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(180.0 * Utilities::const_pi/180 * 10));
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 1);
  CHECK(distance_from_planes["segmentFraction"] == Approx(1.0));

  // curve test 14
  slab_segment_angles[0][0][0] = 0.0 * dtr;
  slab_segment_angles[0][0][1] = 180.0 * dtr;
  slab_segment_angles[0][1][0] = 180.0 * dtr;
  slab_segment_angles[0][1][1] = 270.0 * dtr;
  slab_segment_angles[1][0][0] = 0.0 * dtr;
  slab_segment_angles[1][0][1] = 180.0 * dtr;
  slab_segment_angles[1][1][0] = 180.0 * dtr;
  slab_segment_angles[1][1][1] = 270.0 * dtr;

  slab_segment_lengths[0][0] = 10 * 180 * dtr;
  slab_segment_lengths[0][1] = 10 * 90 * dtr;
  slab_segment_lengths[1][0] = 10 * 180 * dtr;
  slab_segment_lengths[1][1] = 10 * 90 * dtr;

  position[0] = 10;
  position[1] = 0;
  position[2] = 0;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(std::fabs(distance_from_planes["distanceFromPlane"]) < 1e-14); // checked that it should be about 0 this with a drawing
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(90.0 * Utilities::const_pi/180 * 10));
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 0);
  CHECK(distance_from_planes["segmentFraction"] == Approx(0.5));

  // curve test 15
  position[0] = 10;
  position[1] = 10;
  position[2] = -10;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(std::fabs(distance_from_planes["distanceFromPlane"]) < 1e-14); // checked that it should be about 0 this with a drawing
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(180.0 * Utilities::const_pi/180 * 10));
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 0);
  CHECK(distance_from_planes["segmentFraction"] == Approx(1.0));

  // curve test 16
  position[0] = 10;
  position[1] = 10;
  position[2] = -11;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(distance_from_planes["distanceFromPlane"] == Approx(-1.0)); // checked that it should be about -1 this with a drawing
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(180.0 * Utilities::const_pi/180 * 10));
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 0);
  CHECK(distance_from_planes["segmentFraction"] == Approx(1.0));

  // curve test 16
  position[0] = 10;
  position[1] = 10;
  position[2] = -9;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(distance_from_planes["distanceFromPlane"] == Approx(1.0)); // checked that it should be about -1 this with a drawing
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(180.0 * Utilities::const_pi/180 * 10));
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 0);
  CHECK(distance_from_planes["segmentFraction"] == Approx(1.0));

  // curve test 17
  position[0] = 10;
  position[1] = 20;
  position[2] = 0;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(std::fabs(distance_from_planes["distanceFromPlane"]) < 1e-14); // checked that it should be about 0 this with a drawing
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(270.0 * Utilities::const_pi/180 * 10));
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 1);
  CHECK(distance_from_planes["segmentFraction"] == Approx(1.0));


  // curve test 18
  position[0] = 10;
  position[1] = 21;
  position[2] = 0;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(distance_from_planes["distanceFromPlane"] == Approx(-1.0)); // checked that it should be about 1 this with a drawing
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(270.0 * Utilities::const_pi/180 * 10));
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 1);
  CHECK(distance_from_planes["segmentFraction"] == Approx(1.0));

  // curve test 19
  position[0] = 10;
  position[1] = 19;
  position[2] = 0;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(distance_from_planes["distanceFromPlane"] == Approx(1.0)); // checked that it should be about 1 this with a drawing
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(270.0 * Utilities::const_pi/180 * 10));
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 1);
  CHECK(distance_from_planes["segmentFraction"] == Approx(1.0));


  // curve test 20
  slab_segment_angles[0][0][0] = 0.0 * dtr;
  slab_segment_angles[0][0][1] = 270.0 * dtr;
  slab_segment_angles[0][1][0] = 270.0 * dtr;
  slab_segment_angles[0][1][1] = 315.0 * dtr;
  slab_segment_angles[1][0][0] = 0.0 * dtr;
  slab_segment_angles[1][0][1] = 270.0 * dtr;
  slab_segment_angles[1][1][0] = 270.0 * dtr;
  slab_segment_angles[1][1][1] = 315.0 * dtr;

  slab_segment_lengths[0][0] = 10 * 270 * dtr;
  slab_segment_lengths[0][1] = 10 * 45 * dtr;
  slab_segment_lengths[1][0] = 10 * 270 * dtr;
  slab_segment_lengths[1][1] = 10 * 45 * dtr;

  position[0] = 10;
  position[1] = 0;
  position[2] = 0;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(std::fabs(distance_from_planes["distanceFromPlane"]) < 1e-14); // checked that it should be about 0 this with a drawing
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(90.0 * Utilities::const_pi/180 * 10));
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 0);
  CHECK(distance_from_planes["segmentFraction"] == Approx(1.0/3.0));

  // curve test 21
  position[0] = 10;
  position[1] = 10;
  position[2] = -10;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(std::fabs(distance_from_planes["distanceFromPlane"]) < 1e-14); // checked that it should be about 0 this with a drawing
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(180.0 * Utilities::const_pi/180 * 10));
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 0);
  CHECK(distance_from_planes["segmentFraction"] == Approx(2.0/3.0));

  // curve test 21
  position[0] = 10;
  position[1] = 20;
  position[2] = 0;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(std::fabs(distance_from_planes["distanceFromPlane"]) < 1e-14); // checked that it should be about 0 this with a drawing
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(270.0 * Utilities::const_pi/180 * 10));
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 0);
  CHECK(distance_from_planes["segmentFraction"] == Approx(1.0));

  // curve test 22
  position[0] = 10;
  position[1] = 10 + 1e-14 + 10 * sqrt(2)/2; // somehow it doesn't get the exact value here, so adding an epsiolon of 1e-14.
  position[2] = 10 * sqrt(2)/2;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(std::fabs(distance_from_planes["distanceFromPlane"]) < 1e-14); // checked that it should be about 0 this with a drawing
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(315.0 * Utilities::const_pi/180 * 10));
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 1);
  CHECK(distance_from_planes["segmentFraction"] == Approx(1.0));

  // curve test start 45 degree 1
  slab_segment_angles[0][0][0] = 45.0 * dtr;
  slab_segment_angles[0][0][1] = 90.0 * dtr;
  slab_segment_angles[0][1][0] = 90.0 * dtr;
  slab_segment_angles[0][1][1] = 135.0 * dtr;
  slab_segment_angles[1][0][0] = 45.0 * dtr;
  slab_segment_angles[1][0][1] = 90.0 * dtr;
  slab_segment_angles[1][1][0] = 90.0 * dtr;
  slab_segment_angles[1][1][1] = 135.0 * dtr;

  slab_segment_lengths[0][0] = 10 * 45 * dtr;
  slab_segment_lengths[0][1] = 10 * 45 * dtr;
  slab_segment_lengths[1][0] = 10 * 45 * dtr;
  slab_segment_lengths[1][1] = 10 * 45 * dtr;
  slab_segment_lengths[2][0] = 5;
  slab_segment_lengths[2][1] = 5;


  position[0] = 10;
  position[1] = 0;
  position[2] = 0;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(distance_from_planes["distanceFromPlane"] == Approx(-7.3205080757)); // checked that it should be about -7.3 this with a drawing
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(9.5531661812));
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 1);
  CHECK(distance_from_planes["segmentFraction"] == Approx(0.2163468959));

  // curve test change reference point 1
  reference_point[0] = 50;
  reference_point[1] = 50;

  position[0] = 10;
  position[1] = 0;
  position[2] = 0;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  // checked that distanceFromPlane should be infinity (it is on the other side of the circle this with a drawing
  CHECK(distance_from_planes["distanceFromPlane"] == INFINITY);
  CHECK(distance_from_planes["distanceAlongPlane"] == INFINITY);
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.0));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 0);
  CHECK(distance_from_planes["segmentFraction"] == Approx(0.0));

  // curve test change reference point 2
  position[0] = 10;
  position[1] = 10;
  position[2] = 0;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(distance_from_planes["distanceFromPlane"] == Approx(2.3463313527)); // checked that it should be about 2.3 this with a drawing
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(11.780972451));
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 1);
  CHECK(distance_from_planes["segmentFraction"] == Approx(0.5));

  // curve test angle interpolation 1
  reference_point[0] = 0;
  reference_point[1] = 0;

  slab_segment_angles[0][0][0] = 0.0 * dtr;
  slab_segment_angles[0][0][1] = 180.0 * dtr;
  slab_segment_angles[0][1][0] = 180.0 * dtr;
  slab_segment_angles[0][1][1] = 270.0 * dtr;
  slab_segment_angles[1][0][0] = 0.0 * dtr;
  slab_segment_angles[1][0][1] = 90.0 * dtr;
  slab_segment_angles[1][1][0] = 90.0 * dtr;
  slab_segment_angles[1][1][1] = 135.0 * dtr;

  slab_segment_lengths[0][0] = 10 * 135 * dtr;
  slab_segment_lengths[0][1] = 10 * 67.5 * dtr;
  slab_segment_lengths[1][0] = 10 * 135 * dtr;
  slab_segment_lengths[1][1] = 10 * 67.5 * dtr;


  position[0] = 10;
  position[1] = 0;
  position[2] = 0;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(std::fabs(distance_from_planes["distanceFromPlane"]) < 1e-14); // checked that it should be about 0 this with a drawing
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(90.0 * Utilities::const_pi/180 * 10));
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 0);
  CHECK(distance_from_planes["segmentFraction"] == Approx(2.0/3.0));

  // curve test reverse angle 1
  reference_point[0] = 0;
  reference_point[1] = 0;

  slab_segment_angles[0][0][0] = 0.0 * dtr;
  slab_segment_angles[0][0][1] = 90.0 * dtr;
  slab_segment_angles[0][1][0] = 90.0 * dtr;
  slab_segment_angles[0][1][1] = 0.0 * dtr;
  slab_segment_angles[1][0][0] = 0.0 * dtr;
  slab_segment_angles[1][0][1] = 90.0 * dtr;
  slab_segment_angles[1][1][0] = 90.0 * dtr;
  slab_segment_angles[1][1][1] = 0.0 * dtr;

  slab_segment_lengths[0][0] = 10 * 90 * dtr;
  slab_segment_lengths[0][1] = 10 * 90 * dtr;
  slab_segment_lengths[1][0] = 10 * 90 * dtr;
  slab_segment_lengths[1][1] = 10 * 90 * dtr;

  position[0] = 10;
  position[1] = 0;
  position[2] = 0;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(std::fabs(distance_from_planes["distanceFromPlane"]) < 1e-14); // checked that it should be about 0 this with a drawing
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(90.0 * Utilities::const_pi/180 * 10));
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 0);
  CHECK(distance_from_planes["segmentFraction"] == Approx(1.0));

  // curve test reverse angle 2
  position[0] = 10;
  position[1] = -10;
  position[2] = -10;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(std::fabs(distance_from_planes["distanceFromPlane"]) < 1e-14); // checked that it should be about 0 this with a drawing
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(180.0 * Utilities::const_pi/180 * 10));
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 1);
  CHECK(distance_from_planes["segmentFraction"] == Approx(1.0));

  // curve test reverse angle 3
  position[0] = 10;
  position[1] = 10 - (20 - 10 * sqrt(2)/2);
  position[2] = -10 * sqrt(2)/2;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(std::fabs(distance_from_planes["distanceFromPlane"]) < 1e-14); // checked that it should be about 0 this with a drawing
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(135.0 * Utilities::const_pi/180 * 10));
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 1);
  CHECK(distance_from_planes["segmentFraction"] == Approx(0.5));

  // curve test reverse angle 4
  position[0] = 10;

  double angle = 180+0.1;
  position[1] = 10 - (20 * std::cos(0 * Utilities::const_pi/180) + 10 * std::cos((angle) * Utilities::const_pi/180));
  position[2] = 0 * std::cos(0 * Utilities::const_pi/180) + 10 * std::sin((angle) * Utilities::const_pi/180);

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(std::fabs(distance_from_planes["distanceFromPlane"]) < 1e-14); // checked that it should be about 0 this with a drawing
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(90.1 * Utilities::const_pi/180 * 10));
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 1);
  CHECK(distance_from_planes["segmentFraction"] == Approx(0.0011111111));


  // curve test reverse angle 5
  position[0] = 10;
  position[1] = 10 - (20 - 10 * std::cos(0.001 * Utilities::const_pi/180));
  position[2] = - 10 * std::sin(0.001 * Utilities::const_pi/180);

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(std::fabs(distance_from_planes["distanceFromPlane"]) < 1e-14); // checked that it should be about 0 this with a drawing
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(90.001 * Utilities::const_pi/180 * 10));
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 1);
  CHECK(distance_from_planes["segmentFraction"] == Approx(0.000011111111));

  // curve test reverse angle 6
  slab_segment_angles[0][0][0] = 0.0 * dtr;
  slab_segment_angles[0][0][1] = 45.0 * dtr;
  slab_segment_angles[0][1][0] = 45.0 * dtr;
  slab_segment_angles[0][1][1] = 0.0 * dtr;
  slab_segment_angles[1][0][0] = 0.0 * dtr;
  slab_segment_angles[1][0][1] = 45.0 * dtr;
  slab_segment_angles[1][1][0] = 45.0 * dtr;
  slab_segment_angles[1][1][1] = 0.0 * dtr;

  slab_segment_lengths[0][0] = 10 * 45 * dtr;
  slab_segment_lengths[0][1] = 10 * 45 * dtr;
  slab_segment_lengths[1][0] = 10 * 45 * dtr;
  slab_segment_lengths[1][1] = 10 * 45 * dtr;

  position[0] = 10;
  position[1] = 10 - 10 * std::cos(45.000 * Utilities::const_pi/180);
  position[2] = 10 * std::sin(45.000 * Utilities::const_pi/180);

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(std::fabs(distance_from_planes["distanceFromPlane"]) < 1e-14); // checked that it should be about 0 this with a drawing
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(45 * Utilities::const_pi/180 * 10));
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 0);
  CHECK(distance_from_planes["segmentFraction"] == Approx(1.0));

  // curve test reverse angle 6
  position[0] = 10;
  angle = 45;
  position[1] = 10 - (10 * std::cos((angle) * Utilities::const_pi/180));
  position[2] = 10 * std::sin((angle) * Utilities::const_pi/180);

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(std::fabs(distance_from_planes["distanceFromPlane"]) < 1e-14); // checked that it should be about 0 this with a drawing
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(45 * Utilities::const_pi/180 * 10));
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 0.0);
  CHECK(distance_from_planes["segmentFraction"] == Approx(1.0));

  // curve test reverse angle 6
  position[0] = 10;
  angle = 180+45;
  position[1] = 10 - (20 * std::cos(45 * Utilities::const_pi/180) + 10 * std::cos((angle) * Utilities::const_pi/180));
  position[2] = 20 * std::cos(45 * Utilities::const_pi/180) + 10 * std::sin((angle) * Utilities::const_pi/180);

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(std::fabs(distance_from_planes["distanceFromPlane"]) < 1e-14); // checked that it should be about 0 this with a drawing
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(45 * Utilities::const_pi/180 * 10));
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 0);
  CHECK(distance_from_planes["segmentFraction"] == Approx(1.0));


  // curve test reverse angle 7
  position[0] = 10;
  angle = 180+46;
  position[1] = 10 - (20 * std::cos(45 * Utilities::const_pi/180) + 10 * std::cos((angle) * Utilities::const_pi/180));
  position[2] = 20 * std::cos(45 * Utilities::const_pi/180) + 10 * std::sin((angle) * Utilities::const_pi/180);

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(std::fabs(distance_from_planes["distanceFromPlane"]) < 1e-14); // checked that it should be about 0 this with a drawing
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(46 * Utilities::const_pi/180 * 10));
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 1);
  CHECK(distance_from_planes["segmentFraction"] == Approx(0.0222222222));



  // curve test reverse angle 8
  position[0] = 10;
  angle = 180+46;
  position[1] = 10 - (20 * std::cos(45 * Utilities::const_pi/180) + 10 * std::cos((angle) * Utilities::const_pi/180))+0.1;
  position[2] = 20 * std::cos(45 * Utilities::const_pi/180) + 10 * std::sin((angle) * Utilities::const_pi/180);

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(distance_from_planes["distanceFromPlane"] == Approx(0.0697227738)); // checked that it should be small positive this with a drawing
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx((90 - 44.4093) * Utilities::const_pi/180 * 10));
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 1);
  CHECK(distance_from_planes["segmentFraction"] == Approx(0.0131266424));

  // curve test reverse angle 9
  position[0] = 10;
  angle = 180+46;
  position[1] = 10 - (20 * std::cos(45 * Utilities::const_pi/180) + 10 * std::cos((angle) * Utilities::const_pi/180))-0.1;
  position[2] = 20 * std::cos(45 * Utilities::const_pi/180) + 10 * std::sin((angle) * Utilities::const_pi/180);

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(distance_from_planes["distanceFromPlane"] == Approx(-0.0692053058)); // checked that it should be small negative this with a drawing
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx((90 - 43.585) * Utilities::const_pi/180 * 10));
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 1);
  CHECK(distance_from_planes["segmentFraction"] == Approx(0.031445048));

  // curve test reverse angle 10
  position[0] = 10;
  angle = 180+90;
  position[1] = 10 - (20 * std::cos(45 * Utilities::const_pi/180) + 10 * std::cos((angle) * Utilities::const_pi/180));
  position[2] = 20 * std::cos(45 * Utilities::const_pi/180) + 10 * std::sin((angle) * Utilities::const_pi/180);

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false);

  CHECK(std::fabs(distance_from_planes["distanceFromPlane"]) < 1e-14); // checked that it should be about 0 this with a drawing
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(90 * Utilities::const_pi/180 * 10));
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 1);
  CHECK(distance_from_planes["segmentFraction"] == Approx(1.0));


  // global_x_list test 1
  // reminder, coordinates are {0,10},{20,10},{30,10}
  position[0] = 10;
  position[1] = 10;
  position[2] = 10;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false,
  {0,1,2});

  CHECK(std::fabs(distance_from_planes["distanceFromPlane"]) < 1e-14); // checked that it should be about 0 this with a drawing
  CHECK(std::fabs(distance_from_planes["distanceAlongPlane"]) < 1e-14);
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 0);
  CHECK(distance_from_planes["segmentFraction"] == Approx(0.0));

  // global_x_list test 2
  position[0] = 10;
  position[1] = 10;
  position[2] = 10;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false,
  {0,0.5,1});

  CHECK(std::fabs(distance_from_planes["distanceFromPlane"]) < 1e-14); // checked that it should be about 0 this with a drawing
  CHECK(std::fabs(distance_from_planes["distanceAlongPlane"]) < 1e-14);
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.25));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 0);
  CHECK(distance_from_planes["segmentFraction"] == Approx(0.0));

  // global_x_list test 3
  position[0] = 15;
  position[1] = 10;
  position[2] = 10;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false,
  {0,0.5,1});

  CHECK(std::fabs(distance_from_planes["distanceFromPlane"]) < 1e-14); // checked that it should be about 0 this with a drawing
  CHECK(std::fabs(distance_from_planes["distanceAlongPlane"]) < 1e-14);
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.375));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 0);
  CHECK(distance_from_planes["segmentFraction"] == Approx(0.0));

  // global_x_list test 4
  position[0] = 20;
  position[1] = 10;
  position[2] = 10;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false,
  {0,0.5,1});

  CHECK(std::fabs(distance_from_planes["distanceFromPlane"]) < 1e-14); // checked that it should be about 0 this with a drawing
  CHECK(std::fabs(distance_from_planes["distanceAlongPlane"]) < 1e-14);
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 0);
  CHECK(distance_from_planes["segmentFraction"] == Approx(0.0));

  // global_x_list test 5
  position[0] = 25;
  position[1] = 10;
  position[2] = 10;
  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false,
  {0,0.5,1});

  CHECK(std::fabs(distance_from_planes["distanceFromPlane"]) < 1e-14); // checked that it should be about 0 this with a drawing
  CHECK(std::fabs(distance_from_planes["distanceAlongPlane"]) < 1e-12);
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.75));
  CHECK(distance_from_planes["section"] == 1);
  CHECK(distance_from_planes["segment"] == 0);
  CHECK(std::fabs(distance_from_planes["segmentFraction"]) < 1e-12);



  // global_x_list test 6
  position[0] = 30;
  position[1] = 10;
  position[2] = 10;

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 cartesian_system,
                                                 false,
  {0,0.5,1});

  CHECK(std::fabs(distance_from_planes["distanceFromPlane"]) < 1e-14); // checked that it should be about 0 this with a drawing
  CHECK(std::fabs(distance_from_planes["distanceAlongPlane"]) < 1e-14);
  CHECK(distance_from_planes["sectionFraction"] == Approx(1.0));
  CHECK(distance_from_planes["section"] == 1);
  CHECK(distance_from_planes["segment"] == 0);
  CHECK(distance_from_planes["segmentFraction"] == Approx(0.0));

}


TEST_CASE("WorldBuilder Utilities function: distance_point_from_curved_planes spherical")
{
  // Because most functionallity is already tested by the cartesian version
  // of this test case, the scope of this test case is only to test whether
  // the code which is different for the spherical case is correct.

  // spherical test 1
  std::string file_name = WorldBuilder::Data::WORLD_BUILDER_SOURCE_DIR + "/tests/data/subducting_plate_different_angles_spherical.wb";
  WorldBuilder::World world(file_name);

  const double dtr = Utilities::const_pi/180.0;
  Point<3> position(10,0 * dtr,10 * dtr,spherical);
  position = Point<3>(world.parameters.coordinate_system->natural_to_cartesian_coordinates(position.get_array()),cartesian);

  Point<2> reference_point(0,0,spherical);

  std::vector<Point<2> > coordinates;
  coordinates.push_back(Point<2>(0 * dtr,10 * dtr,spherical));
  coordinates.push_back(Point<2>(10 * dtr,10 * dtr,spherical));

  std::vector<std::vector<double> > slab_segment_lengths(2);
  slab_segment_lengths[0].push_back(std::sqrt(10*10+10*10));
  slab_segment_lengths[0].push_back(200);
  slab_segment_lengths[1].push_back(std::sqrt(10*10+10*10));
  slab_segment_lengths[1].push_back(200);

  //double dtr = Utilities::const_pi/180;
  std::vector<std::vector<Point<2> > > slab_segment_angles(2);
  slab_segment_angles[0].push_back(Point<2>(45 * dtr,45 * dtr,cartesian));
  slab_segment_angles[0].push_back(Point<2>(45 * dtr,45 * dtr,cartesian));
  slab_segment_angles[1].push_back(Point<2>(45 * dtr,45 * dtr,cartesian));
  slab_segment_angles[1].push_back(Point<2>(45 * dtr,45 * dtr,cartesian));

  double starting_radius = 10;

  std::map<std::string,double> distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 world.parameters.coordinate_system,
                                                 false);

  CHECK(std::fabs(distance_from_planes["distanceFromPlane"]) < 1e-14); // practically zero
  CHECK(std::fabs(distance_from_planes["distanceAlongPlane"]) < 1e-14);
  CHECK(std::fabs(distance_from_planes["sectionFraction"]) < 1e-14);
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 0);
  CHECK(std::fabs(distance_from_planes["segmentFraction"]) < 1e-14);


  // spherical test 2
  position = Point<3>(10,10 * dtr,10 * dtr,spherical);
  position = Point<3>(world.parameters.coordinate_system->natural_to_cartesian_coordinates(position.get_array()),cartesian);

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 world.parameters.coordinate_system,
                                                 false);

  CHECK(std::fabs(distance_from_planes["distanceFromPlane"]) < 1e-14); // practically zero
  CHECK(std::fabs(distance_from_planes["distanceAlongPlane"]) < 1e-14);
  CHECK(distance_from_planes["sectionFraction"] == Approx(1.0));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 0);
  CHECK(std::fabs(distance_from_planes["segmentFraction"]) < 1e-14);


  // spherical test 2
  coordinates[0][0] = -10 * dtr;
  coordinates[0][1] = 45 * dtr;
  coordinates[1][0] = 10 * dtr;
  coordinates[1][1] = 45 * dtr;
  position = Point<3>(10,0 * dtr,45 * dtr,spherical);
  position = Point<3>(world.parameters.coordinate_system->natural_to_cartesian_coordinates(position.get_array()),cartesian);

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 world.parameters.coordinate_system,
                                                 false);

  CHECK(std::fabs(distance_from_planes["distanceFromPlane"]) < 1e-14); // practically zero
  CHECK(std::fabs(distance_from_planes["distanceAlongPlane"]) < 1e-14);
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 0);
  CHECK(std::fabs(distance_from_planes["segmentFraction"]) < 1e-14);


// spherical test 3
  position = Point<3>(5,0 * dtr,45 * dtr,spherical);
  position = Point<3>(world.parameters.coordinate_system->natural_to_cartesian_coordinates(position.get_array()),cartesian);

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 world.parameters.coordinate_system,
                                                 false);

  CHECK(distance_from_planes["distanceFromPlane"] == Approx(10*sqrt(2)/4)); // checked it with a geometric drawing
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(10*sqrt(2)/4)); // checked it with a geometric drawing
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 0);
  CHECK(distance_from_planes["segmentFraction"] == Approx(0.25));


// spherical test 4
  position = Point<3>(10*sqrt(2)/2,0 * dtr,90 * dtr,spherical);
  position = Point<3>(world.parameters.coordinate_system->natural_to_cartesian_coordinates(position.get_array()),cartesian);

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 world.parameters.coordinate_system,
                                                 false);

  CHECK(distance_from_planes["distanceFromPlane"] == Approx(10*sqrt(2)/2)); // checked it with a geometric drawing
  CHECK(std::fabs(distance_from_planes["distanceAlongPlane"]) < 1e-14); // checked it with a geometric drawing
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 0);
  CHECK(std::fabs(distance_from_planes["segmentFraction"]) < 1e-14);


// spherical test 5
  position = Point<3>(10*sqrt(2)/2,0 * dtr,0 * dtr,spherical);
  position = Point<3>(world.parameters.coordinate_system->natural_to_cartesian_coordinates(position.get_array()),cartesian);

  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 world.parameters.coordinate_system,
                                                 false);

  CHECK(std::fabs(distance_from_planes["distanceFromPlane"]) < 1e-14);  // checked it with a geometric drawing
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(10*sqrt(2)/2)); // checked it with a geometric drawing
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 0);
  CHECK(distance_from_planes["segmentFraction"] == Approx(0.5));

// spherical curve test 1
// This test has not been checked analytically or with a drawing, but
// since the non-curved version works, and the visuals look oke, this
// test is used to see if this changes. Todo: Construct analytical
// solutions to test against.
  slab_segment_angles[0][0][0] = 0.0 * dtr;
  slab_segment_angles[0][0][1] = 45.0 * dtr;
  slab_segment_angles[0][1][0] = 45.0 * dtr;
  slab_segment_angles[0][1][1] = 0.0 * dtr;
  slab_segment_angles[1][0][0] = 0.0 * dtr;
  slab_segment_angles[1][0][1] = 45.0 * dtr;
  slab_segment_angles[1][1][0] = 45.0 * dtr;
  slab_segment_angles[1][1][1] = 0.0 * dtr;

  position = Point<3>(10*sqrt(2)/2,0 * dtr,0 * dtr,spherical);
  position = Point<3>(world.parameters.coordinate_system->natural_to_cartesian_coordinates(position.get_array()),cartesian);
  distance_from_planes =
    Utilities::distance_point_from_curved_planes(position,
                                                 reference_point,
                                                 coordinates,
                                                 slab_segment_lengths,
                                                 slab_segment_angles,
                                                 starting_radius,
                                                 world.parameters.coordinate_system,
                                                 false);

  CHECK(distance_from_planes["distanceFromPlane"] == Approx(4.072033215));  // see comment at the top of the test
  CHECK(distance_from_planes["distanceAlongPlane"] == Approx(6.6085171895)); // see comment at the top of the test
  CHECK(distance_from_planes["sectionFraction"] == Approx(0.5));
  CHECK(distance_from_planes["section"] == 0);
  CHECK(distance_from_planes["segment"] == 0);
  CHECK(distance_from_planes["segmentFraction"] == Approx(0.4672927318));
}

TEST_CASE("WorldBuilder Utilities function: distance_point_from_curved_planes spherical depth methods")
{

  {
    // starting point
    std::string file_name = WorldBuilder::Data::WORLD_BUILDER_SOURCE_DIR + "/tests/data/spherical_depth_method_starting_point.wb";
    WorldBuilder::World world(file_name);

    const double dtr = Utilities::const_pi/180.0;
    // slab goes down and up again
    // origin
    std::array<double,3> position = {{6371000 - 0, 0 * dtr, 0 * dtr}};
    position = world.parameters.coordinate_system->natural_to_cartesian_coordinates(position);
    CHECK(world.temperature(position, 0, 10) == Approx(0));
    CHECK(world.temperature(position, 1, 10) == 0);
    CHECK(world.temperature(position, 200e3, 10) == 0);
    CHECK(world.temperature(position, 210e3, 10) == Approx(1696.9009710498));

    // ~330 km
    position = {{6371000 - 0, 0 * dtr, -3 * dtr}};
    position = world.parameters.coordinate_system->natural_to_cartesian_coordinates(position);
    CHECK(world.temperature(position, 0, 10) == Approx(1600.0));
    CHECK(world.temperature(position, 50e3, 10) == Approx(1622.5575343016));
    CHECK(world.temperature(position, 75e3, 10) == 0);
    CHECK(world.temperature(position, 250e3, 10) == 0);
    CHECK(world.temperature(position, 275e3, 10) == Approx(1728.0673222282));


    // ~1100 km
    position = {{6371000 - 0, 0 * dtr, -10 * dtr}};
    position = world.parameters.coordinate_system->natural_to_cartesian_coordinates(position);
    CHECK(world.temperature(position, 0, 10) == Approx(1600.0));
    CHECK(world.temperature(position, 95e3, 10) == Approx(1643.1311005134));
    CHECK(world.temperature(position, 100e3, 10) == 0);
    CHECK(world.temperature(position, 300e3, 10) == 0);
    CHECK(world.temperature(position, 305e3, 10) == Approx(1742.6442250145));


    // ~2200 km
    position = {{6371000 - 0, 0 * dtr, -20 * dtr}};
    position = world.parameters.coordinate_system->natural_to_cartesian_coordinates(position);
    CHECK(world.temperature(position, 0, 10) == Approx(1600.0));
    CHECK(world.temperature(position, 1, 10) == 0.0);
    CHECK(world.temperature(position, 200e3, 10) == 0);
    CHECK(world.temperature(position, 205e3, 10) == Approx(1694.5269718775));
    CHECK(world.temperature(position, 570e3, 10) == Approx(1876.8664968188));
  }

  {
    // begin segment depth method
    std::string file_name = WorldBuilder::Data::WORLD_BUILDER_SOURCE_DIR + "/tests/data/spherical_depth_method_begin_segment.wb";
    WorldBuilder::World world(file_name);

    const double dtr = Utilities::const_pi/180.0;
    // origin
    std::array<double,3> position = {{6371000 - 0, 0 * dtr, 0 * dtr}};
    position = world.parameters.coordinate_system->natural_to_cartesian_coordinates(position);
    CHECK(world.temperature(position, 0, 10) == Approx(0.0));
    CHECK(world.temperature(position, 1, 10) == 0);
    CHECK(world.temperature(position, 200e3, 10) == 0);
    CHECK(world.temperature(position, 210e3, 10) == Approx(1696.9009710498));

    // ~330 km
    position = {{6371000 - 0, 0 * dtr, -3 * dtr}};
    position = world.parameters.coordinate_system->natural_to_cartesian_coordinates(position);
    CHECK(world.temperature(position, 0, 10) == Approx(1600.0));
    CHECK(world.temperature(position, 50e3, 10) == Approx(1622.5575343016));
    CHECK(world.temperature(position, 75e3, 10) == 0);
    CHECK(world.temperature(position, 250e3, 10) == 0);
    CHECK(world.temperature(position, 275e3, 10) == Approx(1728.0673222282));


    // ~1100 km
    position = {{6371000 - 0, 0 * dtr, -10 * dtr}};
    position = world.parameters.coordinate_system->natural_to_cartesian_coordinates(position);
    CHECK(world.temperature(position, 0, 10) == Approx(1600.0));
    CHECK(world.temperature(position, 150e3, 10) == Approx(1668.6311660012));
    CHECK(world.temperature(position, 175e3, 10) == 0);
    CHECK(world.temperature(position, 380e3, 10) == 0);
    CHECK(world.temperature(position, 385e3, 10) == Approx(1782.119932987));


    // ~1100 km
    position = {{6371000 - 0, 0 * dtr, -20 * dtr}};
    position = world.parameters.coordinate_system->natural_to_cartesian_coordinates(position);
    CHECK(world.temperature(position, 0, 10) == Approx(1600.0));
    CHECK(world.temperature(position, 350e3, 10) == Approx(1764.7404561736));
    CHECK(world.temperature(position, 355e3, 10) == 0);
    CHECK(world.temperature(position, 565e3, 10) == 0);
    CHECK(world.temperature(position, 570e3, 10) == Approx(1876.8664968188));
  }

}

TEST_CASE("WorldBuilder parameters: invalid 1")
{

  std::string file_name = WorldBuilder::Data::WORLD_BUILDER_SOURCE_DIR + "/tests/data/invalid_1.wb";
  CHECK_THROWS_WITH(WorldBuilder::World(file_name), Contains("Invalid keyword: additionalPropertiesInvalid schema: #/test"));

}

