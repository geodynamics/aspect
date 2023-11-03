/*
  Copyright (C) 2019 by the authors of the ASPECT code.

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
#include <aspect/geometry_model/chunk.h>
#include <aspect/geometry_model/initial_topography_model/zero_topography.h>
#include <deal.II/grid/manifold_lib.h>

using namespace dealii;
using namespace aspect;

TEST_CASE("check chunk geometry derivatives")
{
  GeometryModel::internal::ChunkGeometry<2> manifold;

  const auto topo = new InitialTopographyModel::ZeroTopography<2>();
  manifold.initialize(topo);

  manifold.set_max_depth(7500e3);
  manifold.set_min_longitude(0.0);
  manifold.set_min_radius(5620.e3);

  SphericalManifold<2> spherical_manifold;

  Point<2> point1(3300.e3,-2200.e3);

  Point<2> spoint1 = manifold.pull_back(point1);
  std::array<double,2> spoint1analytical = aspect::Utilities::Coordinates::cartesian_to_spherical_coordinates(point1);

  std::array<double,2> spoint2(spoint1analytical);
  spoint2[0] += 1000000.0;
  spoint2[1] += numbers::PI/4;
  Point<2> point2 = aspect::Utilities::Coordinates::spherical_to_cartesian_coordinates<2>(spoint2);

  std::cout << spoint1 << std::endl;
  std::cout << spoint1analytical[0] << " " << spoint1analytical[1] << std::endl;

  Tensor<1,2> chunk_tangent = manifold.get_tangent_vector(point1,point2);
  chunk_tangent /= chunk_tangent.norm();

  Tensor<1,2> spherical_tangent = spherical_manifold.get_tangent_vector(point1,point2);
  spherical_tangent /= spherical_tangent.norm();

  std::cout << std::endl << chunk_tangent << std::endl;
  std::cout << spherical_tangent << std::endl << std::endl;

  std::cout << manifold.get_intermediate_point(point1,point2,0.3) << std::endl;
  std::cout << spherical_manifold.get_intermediate_point(point1,point2,0.3) << std::endl;

  Point<2> point1_back = manifold.push_forward(spoint1);

  std::cout << std::endl << point1 << std::endl;
  std::cout << point1_back << std::endl;



//  REQUIRE(time_temperature.front().second == 1.0);
//  REQUIRE(time_temperature.back().second == 3.0);
//
//  // Remove first entry because it is more than 1.0 away from last entry
//  aspect::TerminationCriteria::internal::trim_time_temperature_list(1.0,time_temperature);
//
//  REQUIRE(time_temperature.front().second == 2.0);
//  REQUIRE(time_temperature.back().second == 3.0);
//
//  // Try to remove first entry, should not change list, since it should keep at least 2 entries
//  aspect::TerminationCriteria::internal::trim_time_temperature_list(0.3,time_temperature);
//
//  REQUIRE(time_temperature.front().second == 2.0);
//  REQUIRE(time_temperature.back().second == 3.0);
}
