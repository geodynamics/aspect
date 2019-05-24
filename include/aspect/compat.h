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

#if !DEAL_II_VERSION_GTE(9,2,0)

#include <deal.II/base/table.h>
#include <deal.II/base/function_lib.h>
namespace aspect
{
  using namespace dealii;
  /**
   * A scalar function that computes its values by (bi-, tri-)linear
   * interpolation from a set of point data that are arranged on a uniformly
   * spaced tensor product mesh. This function is derived from
   * deal.II (dealii/include/deal.II/base/function_lib.h)
   */
  template <int dim>
  class InterpolatedUniformGridData : public dealii::Functions::InterpolatedUniformGridData<dim>
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
      virtual
      Tensor<1, dim>
      gradient(const Point<dim> &p,
               const unsigned int component) const;

  };
}

namespace aspect
{
  template <int dim>
  inline
  InterpolatedUniformGridData<dim>::InterpolatedUniformGridData(
    const std::array<std::pair<double, double>, dim> &interval_endpoints,
    const std::array<unsigned int, dim>              &n_subintervals,
    const Table<dim, double>                         &data_values)
    :
    dealii::Functions::InterpolatedUniformGridData<dim>(&interval_endpoints, &n_subintervals, &data_values)
  {}


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

    // find out where this data point lies
    const TableIndices<dim> ix = table_index_of_point(p);

    Point<dim> dx;
    for (unsigned int d = 0; d < dim; ++d)
      dx[d] = (this->interval_endpoints[d].second - this->interval_endpoints[d].first)/this->n_subintervals[d];

    Point<dim> p_unit;
    for (unsigned int d = 0; d < dim; ++d)
      p_unit[d] =
        std::max(std::min((p[d] - this->coordinate_values[d][ix[d]]) / dx[d], 1.),
                 0.0);

    return gradient_interpolate(this->data_values, ix, p_unit, dx);
  }

}
#endif
