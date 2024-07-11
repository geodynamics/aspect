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

#include "world_builder/config.h"
#include "world_builder/coordinate_systems/invalid.h"
#include "world_builder/gravity_model/uniform.h"
#include "world_builder/parameters.h"
#include "world_builder/point.h"
#include "world_builder/world.h"

using namespace WorldBuilder;
using doctest::Approx;
using doctest::Contains;

TEST_CASE("Gravity uniform")
{

  const std::string file_name = Data::WORLD_BUILDER_SOURCE_DIR + "/tests/data/subducting_plate_different_angles_cartesian.wb";
  World world(file_name);

  CHECK(world.parameters.gravity_model->gravity_norm(Point<3>(10,10,10,CoordinateSystem::cartesian)) == Approx(10.));

  auto vector_cartesian = world.parameters.gravity_model->gravity_vector(Point<3>(10,10,10,CoordinateSystem::cartesian));

  CHECK(vector_cartesian[0] == Approx(0));
  CHECK(vector_cartesian[1] == Approx(0));
  CHECK(vector_cartesian[2] == Approx(-10.));


  const std::string file_name_spherical = Data::WORLD_BUILDER_SOURCE_DIR + "/tests/data/subducting_plate_different_angles_spherical.wb";
  World world_spherical(file_name_spherical);

  auto vector_spherical = world_spherical.parameters.gravity_model->gravity_vector(Point<3>(1,2,4,CoordinateSystem::spherical));

  CHECK(vector_spherical[0] == Approx(-2.1821789024));
  CHECK(vector_spherical[1] == Approx(-4.3643578047));
  CHECK(vector_spherical[2] == Approx(-8.7287156094));
  CHECK(std::sqrt(vector_spherical[0]*vector_spherical[0]+
                  vector_spherical[1]*vector_spherical[1]+
                  vector_spherical[2]*vector_spherical[2]) == Approx(10.));

  world_spherical.parameters.coordinate_system =
    std::unique_ptr<CoordinateSystems::Interface>(new CoordinateSystems::Invalid(nullptr));

  CHECK_THROWS_WITH(world_spherical.parameters.gravity_model->gravity_vector(Point<3>(1,2,4,CoordinateSystem::invalid)),
                    Contains("Invalid coordinate system when using the gravity vector function."));

}