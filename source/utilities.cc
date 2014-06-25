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

#include <aspect/global.h>
#include <aspect/utilities.h>

namespace aspect
{
  namespace Utilities
  {
    template <int dim>
    std_cxx1x::array<double,dim>
    spherical_coordinates(const Point<dim> &position)
    {
      std_cxx1x::array<double,dim> scoord;

      scoord[0] = std::sqrt(position.norm_square()); // R
      scoord[1] = std::atan2(position(1),position(0)); // Phi
      if (scoord[1] < 0.0) scoord[1] = 2*numbers::PI + scoord[1]; // correct phi to [0,2*pi]
      if (dim==3)
        scoord[2] = std::acos(position(2)/std::sqrt(position.norm_square())); // Theta

      return scoord;
    }

    template <int dim>
    Point<dim>
    cartesian_coordinates(const std_cxx1x::array<double,dim> &scoord)
    {
      Point<dim> ccoord;

      ccoord[0] = scoord[0] * std::sin(scoord[2]) * std::cos(scoord[1]); // X
      ccoord[1] = scoord[0] * std::sin(scoord[2]) * std::sin(scoord[1]); // Y
      ccoord[2] = scoord[0] * std::cos(scoord[2]); // Z
      return ccoord;
    }
  }
}

// explicit instantiation of the functions we implement in this file
namespace aspect
{
  namespace Utilities
  {
#define INSTANTIATE(dim) \
    template \
    std_cxx1x::array<double,dim> \
    spherical_coordinates<dim> (const Point<dim> &); \
    \
    template \
    Point<dim> \
    cartesian_coordinates<dim> (const std_cxx1x::array<double,dim> &);

    ASPECT_INSTANTIATE(INSTANTIATE)
  }
}

