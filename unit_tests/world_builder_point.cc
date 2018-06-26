/*
  Copyright (C) 2018 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/

#include "common.h"
#include <world_builder/point.h>

using namespace WorldBuilder;
using Catch::Matchers::Contains;

TEST_CASE("WorldBuilder Point: Initialize point to zero")
{

  const Point<2> p2;
  const Point<3> p3;

  REQUIRE(p2.get_array() == std::array<double,2> {0,0});
  REQUIRE(p3.get_array() == std::array<double,3> {0,0,0});
}

TEST_CASE("WorldBuilder Point: Initialize point from sources")
{
  const Point<2> p2_array(std::array<double,2> {1,2});
  const Point<3> p3_array(std::array<double,3> {1,2,3});

  REQUIRE(p2_array.get_array() == std::array<double,2> {1,2});
  REQUIRE(p3_array.get_array() == std::array<double,3> {1,2,3});

  const Point<2> p2_point(p2_array);
  const Point<3> p3_point(p3_array);

  REQUIRE(p2_point.get_array() == std::array<double,2> {1,2});
  REQUIRE(p3_point.get_array() == std::array<double,3> {1,2,3});

  const Point<2> p2_explicit(1,2);
  const Point<3> p3_explicit(1,2,3);

  REQUIRE(p2_explicit.get_array() == std::array<double,2> {1,2});
  REQUIRE(p3_explicit.get_array() == std::array<double,3> {1,2,3});

  REQUIRE_THROWS_WITH(Point<2>(1,2,3),Contains("Can't use the 3d constructor in 2d."));
  REQUIRE_THROWS_WITH(Point<3>(1,2),Contains("Can't use the 2d constructor in 3d."));

}
