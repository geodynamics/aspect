/*
  Copyright (C) 2021 by the authors of the World Builder code.

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

#include "world_builder/bounding_box.h"

using namespace WorldBuilder;
using doctest::Approx;


TEST_CASE("bounding box 2D")
{
  // Check default constructor
  const BoundingBox<2> bb1;
  CHECK(bb1.center().distance(Point<2>(0,0,CoordinateSystem::cartesian)) < 1e-12);
  CHECK(BoundingBox<2>().center().distance(Point<2>(0,0,CoordinateSystem::cartesian)) < 1e-12);

  CHECK(bb1.point_inside(Point<2>(2,2.5,CoordinateSystem::cartesian)) == true);
  const Point<2> p1 ({{1,1}}, CoordinateSystem::cartesian);
  const Point<2> p2 ({{2,3}}, CoordinateSystem::cartesian);


  CHECK(BoundingBox<2>().center().distance(Point<2>(0,0,CoordinateSystem::cartesian)) < 1e-12);
  const BoundingBox<2> bb2({p1, p2});
  CHECK(bb2.center().distance(Point<2>(1.5,2.,CoordinateSystem::cartesian)) < 1e-12);

  // Check that the function point_inside works as expected
  CHECK(bb2.point_inside(Point<2>(1.5,2,CoordinateSystem::cartesian)) == true);
  CHECK(bb2.point_inside(Point<2>(-1.5,2,CoordinateSystem::cartesian)) == false);
  CHECK(bb2.point_inside(Point<2>(1.5,-2,CoordinateSystem::cartesian)) == false);
  CHECK(bb2.point_inside(Point<2>(-1.5,-2,CoordinateSystem::cartesian)) == false);
}


TEST_CASE("bounding box 3D")
{
  // Check default constructor
  const BoundingBox<3> bb1;
  CHECK(bb1.center().distance(Point<2>(0,0,CoordinateSystem::cartesian)) < 1e-12);
  CHECK(BoundingBox<3>().center().distance(Point<2>(0,0,CoordinateSystem::cartesian)) < 1e-12);

  CHECK(bb1.point_inside(Point<3>(2,2.5,2.5,CoordinateSystem::cartesian)) == true);

  // Check constructor with provided points
  const Point<3> p1 (
  {
    {
      1,1,1
    }
  }, CoordinateSystem::cartesian);
  const Point<3> p2 (
  {
    {
      2,3,4
    }
  }, CoordinateSystem::cartesian);
  const BoundingBox<3> bb2({p1, p2});

  // Check the center function
  CHECK(bb2.center().distance(Point<2>(1.5,2,CoordinateSystem::cartesian)) < 1e-12);

  // Check that the function point_inside works as expected
  CHECK(bb2.point_inside(Point<3>(1.5,2,2.5,CoordinateSystem::cartesian)) == true);
  CHECK(bb2.point_inside(Point<3>(-1.5,2,2.5,CoordinateSystem::cartesian)) == false);
  CHECK(bb2.point_inside(Point<3>(1.5,-2,2.5,CoordinateSystem::cartesian)) == false);
  CHECK(bb2.point_inside(Point<3>(1.5,2,-2.5,CoordinateSystem::cartesian)) == false);

}
