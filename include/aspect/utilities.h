/*
  Copyright (C) 2014 by the authors of the ASPECT code.

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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/


#ifndef __aspect__utilities_h
#define __aspect__utilities_h

#include <deal.II/base/std_cxx1x/array.h>
#include <deal.II/base/point.h>

namespace aspect
{
  /**
   * A namespace for utility functions that might be used in many different
   * places to prevent code duplication.
   */
  namespace Utilities
  {
    using namespace dealii;

    /**
     * Returns spherical coordinates of a cartesian point.
     * The returned array is filled with radius, phi and
     * theta (polar angle). If the dimension is set to 2
     * theta is omitted. Phi is always normalized to
     * [0,2*pi].
     *
     */
    template <int dim>
    std_cxx1x::array<double,dim>
    spherical_coordinates(const Point<dim> &position);

    /**
     * Return the cartesian point of a spherical position
     * defined by radius, phi and theta (polar angle).
     * If the dimension is set to 2 theta is omitted.
     */
    template <int dim>
    Point<dim>
    cartesian_coordinates(const std_cxx1x::array<double,dim> &sposition);

  }
}


#endif
