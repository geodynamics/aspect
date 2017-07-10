/*
  Copyright (C) 2011 - 2016 by the authors of the ASPECT code.

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


#include <aspect/global.h>
#include <aspect/prescribed_stokes_solution/circle.h>

namespace aspect
{
  namespace PrescribedStokesSolution
  {
    template <int dim>
    void Circle<dim>::stokes_solution (const Point<dim> &p, Vector<double> &value) const
    {
      value(0) = -p(1);
      value(1) = p(0);
      if (dim == 3)
        value(2) = 0;
      value(dim) = 0;
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace PrescribedStokesSolution
  {
    ASPECT_REGISTER_PRESCRIBED_STOKES_SOLUTION(Circle,
                                               "circle",
                                               "This value describes a vector field that rotates "
                                               "around the z-axis with constant angular velocity "
                                               "(i.e., with a velocity that increases with "
                                               "distance from the axis). The pressure is set "
                                               "to zero.")
  }
}

