// include utility functions
#include <aspect/utilities.h>

#include <iostream>

using namespace dealii;
using namespace aspect::Utilities;
// Check various conversions between cartesian and spherical coordinates

template <typename T, unsigned int dim>
void check_point(T point1, T point2)
{
  std::cout << std::endl << "Point 1: ";
  for (unsigned int i = 0; i < dim; ++i)
    std::cout << point1[i] << " ";

  std::cout << std::endl << "Point 2: ";
  for (unsigned int i = 0; i < dim; ++i)
    std::cout << point2[i] << " ";

  std::cout << std::endl;
}

int f()
{
  const dealii::Point<2> origin2(0,0);
  const dealii::Point<3> origin3(0,0,0);

  const std_cxx1x::array<double,2> sorigin2 = {{0,0}};
  const std_cxx1x::array<double,3> sorigin3 = {{0,0,0}};

  const dealii::Point<2> one2(1,1);
  const dealii::Point<3> one3(1,1,1);

  const std_cxx1x::array<double,2> sone2 = {{sqrt(2),numbers::PI/4}};
  const std_cxx1x::array<double,3> sone3 = {{sqrt(3),numbers::PI/4,std::acos(1/sqrt(3))}};

  const dealii::Point<3> x(1,0,0);
  const dealii::Point<3> y(0,1,0);
  const dealii::Point<3> z(0,0,1);

  const std_cxx1x::array<double,3> sx = {{1,0,numbers::PI/2}};
  const std_cxx1x::array<double,3> sy = {{1,numbers::PI/2,numbers::PI/2}};
  const std_cxx1x::array<double,3> sz = {{1,0,0}};

  check_point<const std_cxx1x::array<double,2>,2>(spherical_coordinates(origin2),sorigin2);
  check_point<const std_cxx1x::array<double,3>,3>(spherical_coordinates(origin3),sorigin3);
  check_point<const Point<2>,2>(origin2, cartesian_coordinates<2>(sorigin2));
  check_point<const Point<3>,3>(origin3, cartesian_coordinates<3>(sorigin3));

  check_point<const std_cxx1x::array<double,2>,2>(spherical_coordinates(one2),sone2);
  check_point<const std_cxx1x::array<double,3>,3>(spherical_coordinates(one3),sone3);
  check_point<const Point<2>,2>(one2, cartesian_coordinates<2>(sone2));
  check_point<const Point<3>,3>(one3, cartesian_coordinates<3>(sone3));

  check_point<const Point<3>,3>(x, cartesian_coordinates<3>(sx));
  check_point<const Point<3>,3>(y, cartesian_coordinates<3>(sy));
  check_point<const Point<3>,3>(z, cartesian_coordinates<3>(sz));

  check_point<const std_cxx1x::array<double,3>,3>(spherical_coordinates(x),sx);
  check_point<const std_cxx1x::array<double,3>,3>(spherical_coordinates(y),sy);
  check_point<const std_cxx1x::array<double,3>,3>(spherical_coordinates(z),sz);

  const dealii::Point<3> dateline(0,-1,0);
  const std_cxx1x::array<double,3> sdateline = {{1,3*numbers::PI/2,numbers::PI/2}};

  check_point<const Point<3>,3>(dateline, cartesian_coordinates<3>(sdateline));
  check_point<const std_cxx1x::array<double,3>,3>(spherical_coordinates(dateline),sdateline);


  return 42;
}

int i = f();
