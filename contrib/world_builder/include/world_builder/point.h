/*
  Copyright (C) 2018 by the authors of the World Builder code.

  This file is part of the World Builder.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published
   by the Free Software Foundation, either version 3 of the License, or
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

#include <array>
#include <cmath>

#include <world_builder/coordinate_system.h>

namespace WorldBuilder
{
  template<int dim>
  class Point
  {
    public:
      /**
       * Constructor. Constructs a Point at (0,0) in 2d or (0,0,0) in 3d
       * with a Cartesian coordinate system.
       */
      Point();

      /**
       * todo
       */
      Point(const std::array<double,dim> &array, CoordinateSystem coordinate_system = CoordinateSystem::cartesian);

      /**
       * todo
       */
      Point(const Point<dim> &point, CoordinateSystem coordinate_system = CoordinateSystem::cartesian);

      /**
       * todo
       */
      Point(const double x, const double y, CoordinateSystem coordinate_system = CoordinateSystem::cartesian);

      /**
       * todo
       */
      Point(const double x, const double y, const double z, CoordinateSystem coordinate_system = CoordinateSystem::cartesian);

      /**
       * Destructor
       */
      ~Point();

      Point<dim> operator=(const Point<dim> &point);

      /**
       * dot product
       */
      double operator*(const Point<dim> &point) const;


      /**
       * Multiply the vector with a scalar
       */
      Point<dim> operator*(const double scalar) const;

      /**
       * add two points
       */
      Point<dim> operator+(const Point<dim> &point) const;


      /**
       * substract two points
       */
      Point<dim> operator-(const Point<dim> &point) const;



      /**
       * Multiply the vector with a scalar
       */
      Point<dim> &operator*=(const double scalar);

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
       * access index (const)
       */
      const double &operator()(const unsigned int index) const;

      /**
       * access index
       */
      double &operator()(const unsigned int index);


      /**
       * return the internal array which stores the point data.
       */
      std::array<double,dim> get_array() const;


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
}
#endif
