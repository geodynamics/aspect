/*
  Copyright (C) 2015 - 2019 by the authors of the ASPECT code.

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

#ifndef _aspect_compat_h
#define _aspect_compat_h

#include <aspect/global.h>

// C++11 related includes.
#include <array>
#include <functional>
#include <memory>

// for std_cxx14::make_unique:
#include <deal.II/base/std_cxx14/memory.h>

#if !DEAL_II_VERSION_GTE(9,2,0)
#include <deal.II/base/table.h>
#include <deal.II/base/function_lib.h>
namespace aspect
{
  namespace Functions
  {
    using namespace dealii;
    using namespace dealii::Functions;
    /**
     * A scalar function that computes its values by (bi-, tri-)linear
     * interpolation from a set of point data that are arranged on a uniformly
     * spaced tensor product mesh. This function is derived from
     * deal.II (dealii/include/deal.II/base/function_lib.h)
     */
    template <int dim>
    class InterpolatedUniformGridData : public dealii::Function<dim>
    {
      public:
        /**
         * Constructor
         * @param interval_endpoints The left and right end points of the
         * (uniformly subdivided) intervals in each of the coordinate directions.
         * @param n_subintervals The number of subintervals in each coordinate
         * direction. A value of one for a coordinate means that the interval is
         * considered as one subinterval consisting of the entire range. A value
         * of two means that there are two subintervals each with one half of the
         * range, etc.
         * @param data_values A dim-dimensional table of data at each of the mesh
         * points defined by the coordinate arrays above. Note that the Table
         * class has a number of conversion constructors that allow converting
         * other data types into a table where you specify this argument.
         */
        InterpolatedUniformGridData(
          const std::array<std::pair<double, double>, dim> &interval_endpoints,
          const std::array<unsigned int, dim>              &n_subintervals,
          const Table<dim, double>                         &data_values);

        /**
         * Compute the value of the function set by bilinear interpolation of the
         * given data set.
         *
         * @param p The point at which the function is to be evaluated.
         * @param component The vector component. Since this function is scalar,
         * only zero is a valid argument here.
         * @return The interpolated value at this point. If the point lies outside
         * the set of coordinates, the function is extended by a constant.
         */
        Tensor<1, dim>
        gradient(const Point<dim> &p,
                 const unsigned int component = 0) const override;

        /**
         * Compute the value of the function set by bilinear interpolation of the
         * given data set.
         *
         * @param p The point at which the function is to be evaluated.
         * @param component The vector component. Since this function is scalar,
         * only zero is a valid argument here.
         * @return The interpolated value at this point. If the point lies outside
         * the set of coordinates, the function is extended by a constant.
         */
        double
        value(const Point<dim> &p,
              const unsigned int component = 0) const override;


      private:
        /**
         * The set of interval endpoints in each of the coordinate directions.
         */
        const std::array<std::pair<double, double>, dim> interval_endpoints;

        /**
         * The number of subintervals in each of the coordinate directions.
         */
        const std::array<unsigned int, dim> n_subintervals;

        /**
         * The data that is to be interpolated.
         */
        const Table<dim, double> data_values;
    };
  }
}

namespace aspect
{
  namespace Functions
  {
    namespace
    {
      // interpolate a data value from a table where ix denotes
      // the (lower) left endpoint of the interval to interpolate
      // in, and p_unit denotes the point in unit coordinates to do so.
      inline
      double
      interpolate(const Table<1, double> &data_values,
                  const TableIndices<1> &ix,
                  const Point<1>         &xi)
      {
        return ((1 - xi[0]) * data_values[ix[0]] +
                xi[0] * data_values[ix[0] + 1]);
      }

      inline
      double
      interpolate(const Table<2, double> &data_values,
                  const TableIndices<2> &ix,
                  const Point<2>         &p_unit)
      {
        return (((1 - p_unit[0]) * data_values[ix[0]][ix[1]] +
                 p_unit[0] * data_values[ix[0] + 1][ix[1]]) *
                (1 - p_unit[1]) +
                ((1 - p_unit[0]) * data_values[ix[0]][ix[1] + 1] +
                 p_unit[0] * data_values[ix[0] + 1][ix[1] + 1]) *
                p_unit[1]);
      }

      inline
      double
      interpolate(const Table<3, double> &data_values,
                  const TableIndices<3> &ix,
                  const Point<3>         &p_unit)
      {
        return ((((1 - p_unit[0]) * data_values[ix[0]][ix[1]][ix[2]] +
                  p_unit[0] * data_values[ix[0] + 1][ix[1]][ix[2]]) *
                 (1 - p_unit[1]) +
                 ((1 - p_unit[0]) * data_values[ix[0]][ix[1] + 1][ix[2]] +
                  p_unit[0] * data_values[ix[0] + 1][ix[1] + 1][ix[2]]) *
                 p_unit[1]) *
                (1 - p_unit[2]) +
                (((1 - p_unit[0]) * data_values[ix[0]][ix[1]][ix[2] + 1] +
                  p_unit[0] * data_values[ix[0] + 1][ix[1]][ix[2] + 1]) *
                 (1 - p_unit[1]) +
                 ((1 - p_unit[0]) * data_values[ix[0]][ix[1] + 1][ix[2] + 1] +
                  p_unit[0] * data_values[ix[0] + 1][ix[1] + 1][ix[2] + 1]) *
                 p_unit[1]) *
                p_unit[2]);
      }


      // Interpolate the gradient of a data value from a table where ix
      // denotes the lower left endpoint of the interval to interpolate
      // in, p_unit denotes the point in unit coordinates, and dx
      // denotes the width of the interval in each dimension.
      inline
      Tensor<1, 1>
      gradient_interpolate(const Table<1, double> &data_values,
                           const TableIndices<1> &ix,
                           const Point<1>         &p_unit,
                           const Point<1>         &dx)
      {
        (void)p_unit;
        Tensor<1, 1> grad;
        grad[0] = (data_values[ix[0] + 1] - data_values[ix[0]]) / dx[0];
        return grad;
      }


      inline
      Tensor<1, 2>
      gradient_interpolate(const Table<2, double> &data_values,
                           const TableIndices<2> &ix,
                           const Point<2>         &p_unit,
                           const Point<2>         &dx)
      {
        Tensor<1, 2> grad;
        double       u00 = data_values[ix[0]][ix[1]],
                     u01       = data_values[ix[0] + 1][ix[1]],
                     u10       = data_values[ix[0]][ix[1] + 1],
                     u11       = data_values[ix[0] + 1][ix[1] + 1];

        grad[0] =
          ((1 - p_unit[1]) * (u01 - u00) + p_unit[1] * (u11 - u10)) / dx[0];
        grad[1] =
          ((1 - p_unit[0]) * (u10 - u00) + p_unit[0] * (u11 - u01)) / dx[1];
        return grad;
      }

      inline
      Tensor<1, 3>
      gradient_interpolate(const Table<3, double> &data_values,
                           const TableIndices<3> &ix,
                           const Point<3>         &p_unit,
                           const Point<3>         &dx)
      {
        Tensor<1, 3> grad;
        double       u000 = data_values[ix[0]][ix[1]][ix[2]],
                     u001       = data_values[ix[0] + 1][ix[1]][ix[2]],
                     u010       = data_values[ix[0]][ix[1] + 1][ix[2]],
                     u100       = data_values[ix[0]][ix[1]][ix[2] + 1],
                     u011       = data_values[ix[0] + 1][ix[1] + 1][ix[2]],
                     u101       = data_values[ix[0] + 1][ix[1]][ix[2] + 1],
                     u110       = data_values[ix[0]][ix[1] + 1][ix[2] + 1],
                     u111       = data_values[ix[0] + 1][ix[1] + 1][ix[2] + 1];

        grad[0] =
          ((1 - p_unit[2]) *
           ((1 - p_unit[1]) * (u001 - u000) + p_unit[1] * (u011 - u010)) +
           p_unit[2] *
           ((1 - p_unit[1]) * (u101 - u100) + p_unit[1] * (u111 - u110))) /
          dx[0];
        grad[1] =
          ((1 - p_unit[2]) *
           ((1 - p_unit[0]) * (u010 - u000) + p_unit[0] * (u011 - u001)) +
           p_unit[2] *
           ((1 - p_unit[0]) * (u110 - u100) + p_unit[0] * (u111 - u101))) /
          dx[1];
        grad[2] =
          ((1 - p_unit[1]) *
           ((1 - p_unit[0]) * (u100 - u000) + p_unit[0] * (u101 - u001)) +
           p_unit[1] *
           ((1 - p_unit[0]) * (u110 - u010) + p_unit[0] * (u111 - u011))) /
          dx[2];

        return grad;
      }
    } // namespace internal

    template <int dim>
    inline
    InterpolatedUniformGridData<dim>::InterpolatedUniformGridData(
      const std::array<std::pair<double, double>, dim> &interval_endpoints,
      const std::array<unsigned int, dim>              &n_subintervals,
      const Table<dim, double>                         &data_values)
      :
      interval_endpoints(interval_endpoints),
      n_subintervals(n_subintervals),
      data_values(data_values)
    {
      for (unsigned int d = 0; d < dim; ++d)
        {
          Assert(n_subintervals[d] >= 1,
                 ExcMessage("There needs to be at least one subinterval in each "
                            "coordinate direction."));
          Assert(interval_endpoints[d].first < interval_endpoints[d].second,
                 ExcMessage("The interval in each coordinate direction needs "
                            "to have positive size"));
          Assert(data_values.size()[d] == n_subintervals[d] + 1,
                 ExcMessage("The data table does not have the correct size."));
        }
    }


    /**
     * This function is derived from
     * deal.II (dealii/source/base/function_lib.cc)
    */
    template <int dim>
    inline
    Tensor<1, dim>
    InterpolatedUniformGridData<dim>::gradient(
      const Point<dim> &p,
      const unsigned int component) const
    {
      (void)component;
      Assert(
        component == 0,
        ExcMessage(
          "This is a scalar function object, the component can only be zero."));

      // find out where this data point lies, relative to the given
      // subdivision points
      TableIndices<dim> ix;
      for (unsigned int d = 0; d < dim; ++d)
        {
          const double delta_x =
            ((this->interval_endpoints[d].second - this->interval_endpoints[d].first) /
             this->n_subintervals[d]);
          if (p[d] <= this->interval_endpoints[d].first)
            ix[d] = 0;
          else if (p[d] >= this->interval_endpoints[d].second - delta_x)
            ix[d] = this->n_subintervals[d] - 1;
          else
            ix[d] = static_cast<unsigned int>(
                      (p[d] - this->interval_endpoints[d].first) / delta_x);
        }

      // now compute the relative point within the interval/rectangle/box
      // defined by the point coordinates found above. truncate below and
      // above to accommodate points that may lie outside the range
      Point<dim> p_unit;
      Point<dim> delta_x;
      for (unsigned int d = 0; d < dim; ++d)
        {
          delta_x[d] =
            ((this->interval_endpoints[d].second - this->interval_endpoints[d].first) /
             this->n_subintervals[d]);
          p_unit[d] = std::max(std::min((p[d] - this->interval_endpoints[d].first -
                                         ix[d] * delta_x[d]) /
                                        delta_x[d],
                                        1.),
                               0.);
        }

      return gradient_interpolate(this->data_values, ix, p_unit, delta_x);
    }

    template <int dim>
    inline
    double
    InterpolatedUniformGridData<dim>::value(const Point<dim> &p,
                                            const unsigned int component) const
    {
      (void)component;
      Assert(
        component == 0,
        ExcMessage(
          "This is a scalar function object, the component can only be zero."));

      // find out where this data point lies, relative to the given
      // subdivision points
      TableIndices<dim> ix;
      for (unsigned int d = 0; d < dim; ++d)
        {
          const double delta_x =
            ((interval_endpoints[d].second - interval_endpoints[d].first) /
             n_subintervals[d]);
          if (p[d] <= interval_endpoints[d].first)
            ix[d] = 0;
          else if (p[d] >= interval_endpoints[d].second - delta_x)
            ix[d] = n_subintervals[d] - 1;
          else
            ix[d] = static_cast<unsigned int>(
                      (p[d] - interval_endpoints[d].first) / delta_x);
        }

      // now compute the relative point within the interval/rectangle/box
      // defined by the point coordinates found above. truncate below and
      // above to accommodate points that may lie outside the range
      Point<dim> p_unit;
      for (unsigned int d = 0; d < dim; ++d)
        {
          const double delta_x =
            ((interval_endpoints[d].second - interval_endpoints[d].first) /
             n_subintervals[d]);
          p_unit[d] = std::max(std::min((p[d] - interval_endpoints[d].first -
                                         ix[d] * delta_x) /
                                        delta_x,
                                        1.),
                               0.);
        }

      return interpolate(data_values, ix, p_unit);
    }
  }

}
#endif

#endif
