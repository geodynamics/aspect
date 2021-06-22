/*
  Copyright (C) 2018 - 2021 by the authors of the World Builder code.

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
    = default;

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
             "Cannot substract two points which represent different coordinate systems. Internal has type " << static_cast<int>(coordinate_system)
             << ", other point has type " << static_cast<int>(point_.get_coordinate_system()));
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


  template<int dim>
  double
  Point<dim>::distance(const Point<dim> &two) const
  {
    if (this->coordinate_system == spherical)
      {
        // spherical
        const double d_longitude = two[0] - this->point[0];
        const double d_lattitude = two[1] - this->point[1];
        const double sin_d_lat = std::sin(d_lattitude * 0.5);
        const double sin_d_long = std::sin(d_longitude * 0.5);
        return 2.0 * asin(sqrt((sin_d_lat * sin_d_lat) + (sin_d_long*sin_d_long) * std::cos(this->point[1]) * std::cos(two[1])));
      }

    // cartesian
    const double x_distance_to_reference_point = point[0]-two[0];
    const double y_distance_to_reference_point = point[1]-two[1];
    return sqrt((x_distance_to_reference_point*x_distance_to_reference_point) + (y_distance_to_reference_point*y_distance_to_reference_point));

  }



  template<int dim>
  double
  Point<dim>::cheap_relative_distance(const Point<dim> &two) const
  {
    if (this->coordinate_system == spherical)
      {
        // spherical
        const double d_longitude = two[0] - this->point[0];
        const double d_lattitude = two[1] - this->point[1];
        const double sin_d_lat = FT::sin(d_lattitude * 0.5);
        const double sin_d_long = FT::sin(d_longitude * 0.5);
        return (sin_d_lat * sin_d_lat) + (sin_d_long*sin_d_long) * FT::cos(this->point[1]) * FT::cos(two[1]);
      }

    // cartesian
    const double x_distance_to_reference_point = point[0]-two[0];
    const double y_distance_to_reference_point = point[1]-two[1];
    return (x_distance_to_reference_point*x_distance_to_reference_point) + (y_distance_to_reference_point*y_distance_to_reference_point);

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
