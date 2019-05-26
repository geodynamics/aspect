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

#include <limits>
#include <iostream>

#include <world_builder/point.h>
#include <world_builder/assert.h>

namespace WorldBuilder
{
  template<>
  Point<3>::Point(const CoordinateSystem coordinate_system_)
    :
    point({{0,0,0}}),
  coordinate_system(coordinate_system_)
  {}

  template<>
  Point<2>::Point(const CoordinateSystem coordinate_system_)
    :
    point({{0,0}}),
  coordinate_system(coordinate_system_)
  {}

  template<int dim>
  Point<dim>::Point(const std::array<double,dim> &location, const CoordinateSystem coordinate_system_)
    :
    point(location),
    coordinate_system(coordinate_system_)
  {}

  template<int dim>
  Point<dim>::Point(const Point<dim> &point_, const CoordinateSystem coordinate_system_)
    :
    point(point_.get_array()),
    coordinate_system(coordinate_system_)
  {}

  template<int dim>
  Point<dim>::Point(const Point<dim> &point_)
    :
    point(point_.get_array()),
    coordinate_system(point_.get_coordinate_system())
  {}


  template<>
  Point<2>::Point(const double x, const double y, const CoordinateSystem coordinate_system_)
    :
    point({{x,y}}),
  coordinate_system(coordinate_system_)
  {}

  template<>
  Point<3>::Point(const double /*x*/, const double /*y*/, CoordinateSystem coordinate_system_)
    :
    point({{std::numeric_limits<double>::signaling_NaN(),std::numeric_limits<double>::signaling_NaN(), std::numeric_limits<double>::signaling_NaN()}}),
  coordinate_system(coordinate_system_)
  {
    WBAssertThrow(false,"Can't use the 2d constructor in 3d.");
  }


  template<>
  Point<2>::Point(const double /*x*/, const double /*y*/, const double /*z*/, CoordinateSystem coordinate_system_)
    :
    point({{std::numeric_limits<double>::signaling_NaN(),std::numeric_limits<double>::signaling_NaN()}}),
  coordinate_system(coordinate_system_)
  {
    WBAssertThrow(false,"Can't use the 3d constructor in 2d.");
  }


  template<>
  Point<3>::Point(const double x, const double y, const double z, CoordinateSystem coordinate_system_)
    :
    point({{x,y,z}}),
  coordinate_system(coordinate_system_)
  {}


  template<int dim>
  Point<dim>::~Point()
  {}

  template<int dim>
  Point<dim> &Point<dim>::operator=(const Point<dim> &point_)
  {
    point = point_.point;
    coordinate_system = point_.coordinate_system;
    return *this;
  }

  template<int dim>
  double Point<dim>::operator*(const Point<dim> &point_) const
  {
    const std::array<double,dim> array = point_.get_array();
    double dot_product = 0;
    for (unsigned int i = 0; i < dim; ++i)
      dot_product += point[i] * array[i];
    return dot_product;
  }

  template<int dim>
  Point<dim> Point<dim>::operator*(const double scalar) const
  {
    // initialize the array to zero.
    std::array<double,dim> array = Point<dim>(coordinate_system).get_array();
    for (unsigned int i = 0; i < dim; ++i)
      array[i] += point[i] * scalar;
    return Point<dim>(array,coordinate_system);
  }

  template<int dim>
  Point<dim> Point<dim>::operator/(const double scalar) const
  {
    // initialize the array to zero.
    std::array<double,dim> array = Point<dim>(coordinate_system).get_array();
    const double one_over_scalar = 1/scalar;
    for (unsigned int i = 0; i < dim; ++i)
      array[i] += point[i] * one_over_scalar;
    return Point<dim>(array,coordinate_system);
  }

  template<int dim>
  Point<dim>
  Point<dim>::operator+(const Point<dim> &point_) const
  {
    WBAssert(coordinate_system == point_.get_coordinate_system(),
             "Cannot add two points which represent different coordinate systems.");
    Point<dim> point_tmp(point,coordinate_system);
    point_tmp += point_;

    return point_tmp;
  }

  template<int dim>
  Point<dim>
  Point<dim>::operator-(const Point<dim> &point_) const
  {
    WBAssert(coordinate_system == point_.get_coordinate_system(),
             "Cannot substract two points which represent different coordinate systems. Internal has type " << (int)coordinate_system << ", other point has type " << (int)point_.get_coordinate_system());
    Point<dim> point_tmp(point,coordinate_system);
    point_tmp -= point_;

    return point_tmp;
  }


  template<int dim>
  Point<dim> &
  Point<dim>::operator*=(const double scalar)
  {
    for (unsigned int i = 0; i < dim; ++i)
      point[i] *= scalar;
    return *this;
  }


  template<int dim>
  Point<dim> &
  Point<dim>::operator/=(const double scalar)
  {
    for (unsigned int i = 0; i < dim; ++i)
      point[i] /= scalar;
    return *this;
  }


  template<int dim>
  Point<dim> &
  Point<dim>::operator+=(const Point<dim> &point_)
  {
    for (unsigned int i = 0; i < dim; ++i)
      point[i] += point_[i];
    return *this;
  }


  template<int dim>
  Point<dim> &
  Point<dim>::operator-=(const Point<dim> &point_)
  {
    for (unsigned int i = 0; i < dim; ++i)
      point[i] -= point_[i];
    return *this;
  }


  /**
   * access index
   */
  template<int dim>
  const double &
  Point<dim>::operator[](const unsigned int index) const
  {
    WBAssertThrow(index <= dim, "Can't ask for element " << index << " in an point with dimension " << dim << ".");
    return point[index];
  }


  /**
   * access index
   */
  template<int dim>
  double &
  Point<dim>::operator[](const unsigned int index)
  {
    WBAssertThrow(index <= dim, "Can't ask for element " << index << " in an point with dimension " << dim << ".");
    return point[index];
  }


  template<int dim>
  const std::array<double,dim> &
  Point<dim>::get_array() const
  {
    return point;
  }


  template<int dim>
  CoordinateSystem
  Point<dim>::get_coordinate_system() const
  {
    return coordinate_system;
  }


  template<int dim>
  double
  Point<dim>::norm() const
  {
    return std::sqrt(this->norm_square());
  }


  template<>
  double
  Point<2>::norm_square() const
  {
    return (point[0] * point[0]) + (point[1] * point[1]);
  }

  template<>
  double
  Point<3>::norm_square() const
  {
    return (point[0] * point[0]) + (point[1] * point[1]) + (point[2] * point[2]);
  }


  /**
   * Multiplies a point with a scalar.
   */
  template<int dim>
  Point<dim>
  operator*(const double scalar, const Point<dim> &point)
  {
    return point*scalar;
  }

  /**
   * Divides a scalar by a point: output_vector[i] = scalar / point[i].
   */
  template<int dim>
  Point<dim>
  operator/(const double scalar, const Point<dim> &point)
  {
    // initialize the array to zero.
    std::array<double,dim> array = Point<dim>(point.coordinate_system).get_array();
    for (unsigned int i = 0; i < dim; ++i)
      array[i] = scalar / point[i];
    return Point<dim>(array,point.coordinate_system);
  }

  /**
   * The 2d version of the point class.
   */
  template class Point<2>;


  /**
   * The 3d version of the point class.
   */
  template class Point<3>;

  /**
   * Multiplies a 2d point with a scalar.
   */
  template Point<2> operator*(const double scalar, const Point<2> &point);

  /**
   * Multiplies a 3d point with a scalar.
   */
  template Point<3> operator*(const double scalar, const Point<3> &point);
}
