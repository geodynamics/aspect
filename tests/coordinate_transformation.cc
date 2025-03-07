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

// include utility functions
#include <aspect/utilities.h>

#include <iostream>

namespace aspect
{
  using namespace aspect::Utilities;
// Check various conversions between cartesian and spherical coordinates

  template <typename T, unsigned int dim>
  void check_point(T point1, T point2)
  {
    std::cout << std::endl << "Point 1: ";
    for (unsigned int i = 0; i < dim; ++i)
      std::cout << point1[i] << ' ';

    std::cout << std::endl << "Point 2: ";
    for (unsigned int i = 0; i < dim; ++i)
      std::cout << point2[i] << ' ';

    std::cout << std::endl;
  }

  int f()
  {
    const dealii::Point<2> origin2(0,0);
    const dealii::Point<3> origin3(0,0,0);

    const std::array<double,2> sorigin2 = {{0,0}};
    const std::array<double,3> sorigin3 = {{0,0,0}};

    const dealii::Point<2> one2(1,1);
    const dealii::Point<3> one3(1,1,1);

    const std::array<double,2> sone2 = {{std::sqrt(2),numbers::PI/4}};
    const std::array<double,3> sone3 = {{std::sqrt(3),numbers::PI/4,std::acos(1/std::sqrt(3))}};

    const dealii::Point<3> x(1,0,0);
    const dealii::Point<3> y(0,1,0);
    const dealii::Point<3> z(0,0,1);

    const std::array<double,3> sx = {{1,0,numbers::PI/2}};
    const std::array<double,3> sy = {{1,numbers::PI/2,numbers::PI/2}};
    const std::array<double,3> sz = {{1,0,0}};

    check_point<const std::array<double,2>,2>(Coordinates::cartesian_to_spherical_coordinates(origin2),sorigin2);
    check_point<const std::array<double,3>,3>(Coordinates::cartesian_to_spherical_coordinates(origin3),sorigin3);
    check_point<const Point<2>,2>(origin2, Coordinates::spherical_to_cartesian_coordinates<2>(sorigin2));
    check_point<const Point<3>,3>(origin3, Coordinates::spherical_to_cartesian_coordinates<3>(sorigin3));

    check_point<const std::array<double,2>,2>(Coordinates::cartesian_to_spherical_coordinates(one2),sone2);
    check_point<const std::array<double,3>,3>(Coordinates::cartesian_to_spherical_coordinates(one3),sone3);
    check_point<const Point<2>,2>(one2, Coordinates::spherical_to_cartesian_coordinates<2>(sone2));
    check_point<const Point<3>,3>(one3, Coordinates::spherical_to_cartesian_coordinates<3>(sone3));

    check_point<const Point<3>,3>(x, Coordinates::spherical_to_cartesian_coordinates<3>(sx));
    check_point<const Point<3>,3>(y, Coordinates::spherical_to_cartesian_coordinates<3>(sy));
    check_point<const Point<3>,3>(z, Coordinates::spherical_to_cartesian_coordinates<3>(sz));

    check_point<const std::array<double,3>,3>(Coordinates::cartesian_to_spherical_coordinates(x),sx);
    check_point<const std::array<double,3>,3>(Coordinates::cartesian_to_spherical_coordinates(y),sy);
    check_point<const std::array<double,3>,3>(Coordinates::cartesian_to_spherical_coordinates(z),sz);

    const dealii::Point<3> dateline(0,-1,0);
    const std::array<double,3> sdateline = {{1,3*numbers::PI/2,numbers::PI/2}};

    check_point<const Point<3>,3>(dateline, Coordinates::spherical_to_cartesian_coordinates<3>(sdateline));
    check_point<const std::array<double,3>,3>(Coordinates::cartesian_to_spherical_coordinates(dateline),sdateline);


    return 42;
  }

  int i = f();
}
