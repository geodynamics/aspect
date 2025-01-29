/*
  Copyright (C) 2022 - 2024 by the authors of the ASPECT code.

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

#include <aspect/geometry_model/ellipsoidal_chunk.h>
#include <aspect/geometry_model/initial_topography_model/zero_topography.h>

#include <deal.II/base/exceptions.h>

#include <iostream>

using namespace aspect;
bool test_point(const GeometryModel::internal::EllipsoidalChunkGeometry<3> ellipsoidal_manifold,
                const Point<3> &test_point)
{
  const Point<3> converted_point = ellipsoidal_manifold.pull_back(test_point);
  const Point<3> twice_converted_point = ellipsoidal_manifold.push_forward(converted_point);

  const Tensor<1,3> difference = twice_converted_point - test_point;
  const double error = difference.norm() / test_point.norm();

  std::cout << "Point: " << test_point
            << " becomes: " << converted_point
            << ". Back transformed: " << twice_converted_point
            << ". Error: " << error << std::endl;

  return true;
}

/*
 * Launch the following function when this plugin is created.
 */
int f()
{

  const unsigned int dim=3;

  std::vector<Point<2>> corners(2,Point<2>(-15.0,-15.0));
  corners[1] *= -1.0;

  std::vector<Point<3>> test_points;
  test_points.push_back(Point<3> (6371000.0,0,0));
  test_points.push_back(Point<3> (6171000.0,0,0));
  test_points.push_back(Point<3> (3000000.0,3000000.0,0));
  test_points.push_back(Point<3> (3000000.0,3000000.0,3000000.0));
  test_points.push_back(Point<3> (3000000.0,-3000000.0,-3000000.0));
  test_points.push_back(Point<3> (25000.0,25000.0,3481000.0));
  test_points.push_back(Point<3> (25000.0,25000.0,1000000.0));
  test_points.push_back(Point<3> (25000.0,25000.0,500000.0));
  test_points.push_back(Point<3> (25000.0,25000.0,50000.0));


  InitialTopographyModel::ZeroTopography<dim> topography;
  {
    std::cout << "Simple sphere test" << std::endl;
    GeometryModel::internal::EllipsoidalChunkGeometry<dim> ellipsoidal_manifold(topography,
                                                                                6371000.0,
                                                                                0.0,
                                                                                6371000.0,
                                                                                2890000.0,
                                                                                corners);

    for (unsigned int i=0; i<test_points.size(); ++i)
      test_point(ellipsoidal_manifold,test_points[i]);
  }

  const double semi_major_axis_a = 6378137.0;
  const double eccentricity = 8.1819190842622e-2;
  const double semi_minor = std::sqrt((1 - dealii::Utilities::fixed_power<2>(eccentricity)) * dealii::Utilities::fixed_power<2>(semi_major_axis_a));
  {
    std::cout << "WGS84 test" << std::endl;

    GeometryModel::internal::EllipsoidalChunkGeometry<dim> ellipsoidal_manifold(topography,
                                                                                semi_major_axis_a,
                                                                                eccentricity,
                                                                                semi_minor,
                                                                                2890000.0,
                                                                                corners);

    for (unsigned int i=0; i<test_points.size(); ++i)
      test_point(ellipsoidal_manifold,test_points[i]);
  }

  {
    std::cout << "Test for points outside of the provided depth range" << std::endl;
    GeometryModel::internal::EllipsoidalChunkGeometry<dim> ellipsoidal_manifold(topography,
                                                                                semi_major_axis_a,
                                                                                eccentricity,
                                                                                semi_minor,
                                                                                500000.0,
                                                                                corners);

    for (unsigned int i=0; i<test_points.size(); ++i)
      test_point(ellipsoidal_manifold,test_points[i]);
  }

  exit (0);
  return 42;
}


// run this function by initializing a global variable by it
int i = f();
