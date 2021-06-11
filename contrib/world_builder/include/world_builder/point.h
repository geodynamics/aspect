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

#ifndef _world_builder_point_h
#define _world_builder_point_h
#define _USE_MATH_DEFINES
#include <cmath>
#include <array>

#include <world_builder/coordinate_system.h>

namespace WorldBuilder
{
  /**
   * A class which stores 2d and 3d arrays of doubles (depending on the dimension),
   * and the coordinate system which the coordinates can be used for. It also
   * implements several operations such as the computation of the l2 norm and the
   * dot product.
   */
  template<int dim>
  class Point
  {
    public:
      /**
       * Constructor. Constructs a Point at (0,0) in 2d or (0,0,0) in 3d
       * with a Cartesian coordinate system.
       */
      Point(CoordinateSystem coordinate_system);

      /**
       * Constructor. Constructs a Point from a std::array<double,dim> and
       * a coordinate system.
       */
      Point(const std::array<double,dim> &location, CoordinateSystem coordinate_system);

      /**
       * Constructor. Constructs a Point from an other Point and
       * a coordinate system.
       */
      Point(const Point<dim> &point, CoordinateSystem coordinate_system);

      /**
       * Constructor. Constructs a Point from an other Point.
       */
      Point(const Point<dim> &point);

      /**
       * Constructor. Constructs a 2d Point from two doubles and
       * a coordinate system.
       */
      Point(const double x, const double y, CoordinateSystem coordinate_system);

      /**
       * Constructor. Constructs a 3d Point from three doubles and
       * a coordinate system.
       */
      Point(const double x, const double y, const double z, CoordinateSystem coordinate_system);

      /**
       * Destructor
       */
      ~Point();

      Point<dim> &operator=(const Point<dim> &point);

      /**
       * dot product
       */
      double operator*(const Point<dim> &point) const;


      /**
       * Multiply the vector with a scalar
       */
      Point<dim> operator*(const double scalar) const;

      /**
       * Divide the vector through a scalar
       */
      Point<dim> operator/(const double scalar) const;

      /**
       * add two points
       */
      Point<dim> operator+(const Point<dim> &point) const;


      /**
       * Substract two points
       */
      Point<dim> operator-(const Point<dim> &point) const;



      /**
       * Multiply the vector with a scalar
       */
      Point<dim> &operator*=(const double scalar);

      /**
       * Divide the vector through a scalar
       */
      Point<dim> &operator/=(const double scalar);

      /**
       * add two points
       */
      Point<dim> &operator+=(const Point<dim> &point);


      /**
       * substract two points
       */
      Point<dim> &operator-=(const Point<dim> &point);

      /**
       * access index (const)
       */
      const double &operator[](const unsigned int index) const;


      /**
       * access index
       */
      double &operator[](const unsigned int index);


      /**
       * return the internal array which stores the point data.
       */
      const std::array<double,dim> &get_array() const;


      /**
       * returns the coordinate system associated with the data.
       */
      CoordinateSystem get_coordinate_system() const;


      /**
      * Computes the L2 norm: sqrt(x_i * x_i + y_i * y_i + z_i * z_i) in 3d.
      */
      double norm() const;


      /**
      * Computes the square of the norm, which is the sum of the absolute squares
      * x_i * x_i + y_i * y_i + z_i * z_i in 3d.
      */
      double norm_square() const;


    private:
      std::array<double,dim> point;
      CoordinateSystem coordinate_system;

  };

  template<int dim>
  Point<dim> operator*(const double scalar, const Point<dim> &point);

  template<int dim>
  Point<dim> operator/(const double scalar, const Point<dim> &point);
}
#endif
