#ifndef _glm_h
#define _glm_h
#include <array>
#include <math.h>
#include <limits>

namespace WorldBuilder
{
  /**
   * The code in this file is derived from the OpenGL Mathematics (GLM) library, which
   * is available under two version of MIT Lincense (see below). The derived work in
   * this file is avaible under the same licenses.
   *
   * The changes are mainly specializing and isolation functions which are needed in
   * the World Builder.
   *
   * ================================================================================
   * OpenGL Mathematics (GLM)
   * --------------------------------------------------------------------------------
   * GLM is licensed under The Happy Bunny License and MIT License
   *
   * ================================================================================
   * The Happy Bunny License (Modified MIT License)
   * --------------------------------------------------------------------------------
   * Copyright (c) 2005 - 2014 G-Truc Creation
   *
   * Permission is hereby granted, free of charge, to any person obtaining a copy
   * of this software and associated documentation files (the "Software"), to deal
   * in the Software without restriction, including without limitation the rights
   * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   * copies of the Software, and to permit persons to whom the Software is
   * furnished to do so, subject to the following conditions:
   *
   * The above copyright notice and this permission notice shall be included in
   * all copies or substantial portions of the Software.
   *
   * Restrictions:
   *  By making use of the Software for military purposes, you choose to make a
   *  Bunny unhappy.
   *
   * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   * THE SOFTWARE.
   *
   * ================================================================================
   * The MIT License
   * --------------------------------------------------------------------------------
   * Copyright (c) 2005 - 2014 G-Truc Creation
   *
   * Permission is hereby granted, free of charge, to any person obtaining a copy
   * of this software and associated documentation files (the "Software"), to deal
   * in the Software without restriction, including without limitation the rights
   * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   * copies of the Software, and to permit persons to whom the Software is
   * furnished to do so, subject to the following conditions:
   *
   * The above copyright notice and this permission notice shall be included in
   * all copies or substantial portions of the Software.
   *
   * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   * THE SOFTWARE.
   */
  namespace glm
  {
    namespace quaternion
    {
      struct quat
      {
        quat(double w_, double x_, double y_, double z_)
          :
          w(w_),
          x(x_),
          y(y_),
          z(z_)
        {}


        const quat operator-() const
        {
          return quat(-w, -x, -y, -z);
        }

        const quat operator-(quat &p) const
        {
          return quat(p.w-w, p.x-x, p.y-y, p.z-z);
        }

        quat operator*(double const &s) const
        {
          return quat(
                   w * s, x * s, y * s, z * s);
        }


        quat operator/( double const &s)
        {
          return quat(
                   w / s, x / s, y / s, z / s);
        }

        quat operator+(quat const &p)
        {
          return quat(w + p.w, x + p.x, y + p.y, z + p.z);
        }

        double w, x, y, z;
      };

      inline quat operator*(double const &s,quat const &q)
      {
        return quat(
                 q.w * s, q.x * s, q.y * s, q.z * s);
      }



      inline double dot(quat u, quat v)
      {
        return u.w * v.w + u.x * v.x + u.y * v.y + u.z * v.z;
      }

      inline double mix(double const     &x,
                        double const     &y,
                        double const     &a )
      {
        return x * (1.0 - a) + y * a;
      }

      inline quat quat_cast(std::array<std::array<double,3>,3> const &m)
      {
        double fourXSquaredMinus1 = m[0][0] - m[1][1] - m[2][2];
        double fourYSquaredMinus1 = m[1][1] - m[0][0] - m[2][2];
        double fourZSquaredMinus1 = m[2][2] - m[0][0] - m[1][1];
        double fourWSquaredMinus1 = m[0][0] + m[1][1] + m[2][2];

        int biggestIndex = 0;
        double fourBiggestSquaredMinus1 = fourWSquaredMinus1;
        if (fourXSquaredMinus1 > fourBiggestSquaredMinus1)
          {
            fourBiggestSquaredMinus1 = fourXSquaredMinus1;
            biggestIndex = 1;
          }
        if (fourYSquaredMinus1 > fourBiggestSquaredMinus1)
          {
            fourBiggestSquaredMinus1 = fourYSquaredMinus1;
            biggestIndex = 2;
          }
        if (fourZSquaredMinus1 > fourBiggestSquaredMinus1)
          {
            fourBiggestSquaredMinus1 = fourZSquaredMinus1;
            biggestIndex = 3;
          }

        double biggestVal = sqrt(fourBiggestSquaredMinus1 + static_cast<double>(1)) * static_cast<double>(0.5);
        double mult = static_cast<double>(0.25) / biggestVal;

        switch (biggestIndex)
          {
            case 0:
              return quat {biggestVal, (m[1][2] - m[2][1]) *mult, (m[2][0] - m[0][2]) *mult, (m[0][1] - m[1][0]) *mult};
            case 1:
              return quat {(m[1][2] - m[2][1]) *mult, biggestVal, (m[0][1] + m[1][0]) *mult, (m[2][0] + m[0][2]) *mult};
            case 2:
              return quat {(m[2][0] - m[0][2]) *mult, (m[0][1] + m[1][0]) *mult, biggestVal, (m[1][2] + m[2][1]) *mult};
            case 3:
              return quat {(m[0][1] - m[1][0]) *mult, (m[2][0] + m[0][2]) *mult, (m[1][2] + m[2][1]) *mult, biggestVal};
            default: // Silence a -Wswitch-default warning in GCC. Should never actually get here. Assert is just for sanity.
              throw (false);
              return quat {1, 0, 0, 0};
          }
      }



      inline std::array<std::array<double,3>,3> mat3_cast(quat const &q)
      {
        std::array<std::array<double,3>,3> Result;
        double qxx(q.x * q.x);
        double qyy(q.y * q.y);
        double qzz(q.z * q.z);
        double qxz(q.x * q.z);
        double qxy(q.x * q.y);
        double qyz(q.y * q.z);
        double qwx(q.w * q.x);
        double qwy(q.w * q.y);
        double qwz(q.w * q.z);

        Result[0][0] = double(1) - double(2) * (qyy +  qzz);
        Result[0][1] = double(2) * (qxy + qwz);
        Result[0][2] = double(2) * (qxz - qwy);

        Result[1][0] = double(2) * (qxy - qwz);
        Result[1][1] = double(1) - double(2) * (qxx +  qzz);
        Result[1][2] = double(2) * (qyz + qwx);

        Result[2][0] = double(2) * (qxz + qwy);
        Result[2][1] = double(2) * (qyz - qwx);
        Result[2][2] = double(1) - double(2) * (qxx +  qyy);
        return Result;
      }

      inline quat slerp(quat const &x, quat const &y, double a)
      {
        static_assert(std::numeric_limits<double>::is_iec559, "'slerp' only accept floating-point inputs");

        quat z = y;

        double cosTheta = dot(x, y);

        // If cosTheta < 0, the interpolation will take the long way around the sphere.
        // To fix this, one quat must be negated.
        if (cosTheta < static_cast<double>(0))
          {
            z = -y;

            cosTheta = -cosTheta;
          }

        // Perform a linear interpolation when cosTheta is close to 1 to avoid side effect of sin(angle) becoming a zero denominator
        if (cosTheta > 1.0 - std::numeric_limits<double>::epsilon())
          {
            // Linear interpolation
            return quat(
                     mix(x.w, z.w, a),
                     mix(x.x, z.x, a),
                     mix(x.y, z.y, a),
                     mix(x.z, z.z, a));
          }
        else
          {
            // Essential Mathematics, page 467
            double angle = acos(cosTheta);
            return (sin((1.0 - a) * angle) * x + sin(a * angle) * z) / sin(angle);
          }
      }
    }
  }
}
#endif