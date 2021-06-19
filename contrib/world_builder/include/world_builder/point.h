/*
  Copyright (C) 2018-2021 by the authors of the World Builder code.

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
#include <iostream>

#include <world_builder/coordinate_system.h>
#include <world_builder/assert.h>

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
      inline
      const double &operator[](const unsigned int index) const
      {
        WBAssert(index <= dim, "Can't ask for element " << index << " in an point with dimension " << dim << ".");
        return point[index];
      }


      /**
       * access index
       */
      inline double &operator[](const unsigned int index)
      {
        WBAssert(index <= dim, "Can't ask for element " << index << " in an point with dimension " << dim << ".");
        return point[index];
      }

      /**
       * Computes the distance between this and a given point.
       * In spherical coordinates it returns the central angle in radians.
       */
      double
      distance(const Point<dim> &two) const;

      /**
       * Computes the cheapest relative distance between this and a given point.
       * The return value itself is only guartenteed to have the property that a
       * larger value is further away.
       * In the current implementation that means for the cartasian case the squared
       * value is returned and for the spherical value the result of the havearsine
       * function without asin and sqrt is returned.
       */
      double
      cheap_relative_distance(const Point<dim> &two) const;

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

      /**
       * Outputs the values of the point to std cout separated by spaces. This does not
       * output the coordinate system.
       */
      friend std::ostream &operator<<( std::ostream &output, const Point<dim> &point )
      {
        for (unsigned int i = 0; i < dim-1; i++)
          {
            output <<  point[i] << " ";
          }
        output << point[dim-1];

        return output;
      }


    private:
      std::array<double,dim> point;
      CoordinateSystem coordinate_system;

  };

  /**
   * This namespace contains some faster but less accurate version of the
   * trigonomic functions and a faster version of the fmod function.
   */
  namespace FT
  {
    constexpr double const_pi = 3.141592653589793238462643383279502884;

    /**
     * Fast version of the fmod function.
     */
    inline double fmod(const double x, const double y)
    {
      const double x_div_y = x/y;
      return (x_div_y-(int)x_div_y)*y;
    }

    /**
     * Fast sin function, accurate for values between 0 and pi. The implemenation is
     * based on discussion at https://stackoverflow.com/a/6104692.
     *
     * The accuracy seem good enough for most purposes. The unit test tests in steps
     * of 0.01 from -4 pi to 4 pi and compares against the std sin function and the difference
     * is always smaller than 1.2e-5. If the test is run with intervals of 0.001 then there
     * are 12 entries which are (very slightly) above that (<3e-8) at angles of about
     * -174, -6, 6  and 174.
     *
     */
    inline double fast_sin_d(const double angle)
    {
      constexpr double A = 4.0/(const_pi *const_pi);
      constexpr double oneminPmin = 1.-0.1952403377008734-0.01915214119105392;

      const double y = A* angle * ( const_pi - angle );
      return y*( oneminPmin + y*( 0.1952403377008734 + y * 0.01915214119105392 ) ) ;
    }

    /**
     * Fast but less accurate sin function for any angle.
     * Implemented by calling fast_sin_d with a mirrored x if needed to
     * forfill the constrained of fast_sin_d to only have values between
     * zero and pi.
     */
    inline double sin(const double raw_angle)
    {
      const double angle = (raw_angle > -const_pi && raw_angle < const_pi)
                           ?
                           raw_angle
                           :
                           FT::fmod(raw_angle + std::copysign(const_pi,raw_angle), const_pi * 2.0) - std::copysign(const_pi,raw_angle);

      if (angle >= 0)
        return fast_sin_d(angle);
      else
        return -fast_sin_d(-angle);
    }

    /**
     * Fast but less accurate cos function for any angle.
     */
    inline double cos(const double angle)
    {
      return FT::sin((const_pi*0.5)-angle);
    }
  }


  template<int dim>
  Point<dim> operator*(const double scalar, const Point<dim> &point);

  template<int dim>
  Point<dim> operator/(const double scalar, const Point<dim> &point);
}
#endif
