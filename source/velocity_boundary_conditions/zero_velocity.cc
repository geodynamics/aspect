/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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


#include <aspect/velocity_boundary_conditions/zero_velocity.h>


namespace aspect
{
  namespace VelocityBoundaryConditions
  {
    template <int dim>
    Tensor<1,dim>
    ZeroVelocity<dim>::
    boundary_velocity (const Point<dim> &) const
    {
      // return a zero tensor regardless of position
      return Tensor<1,dim>();
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace VelocityBoundaryConditions
  {
    ASPECT_REGISTER_VELOCITY_BOUNDARY_CONDITIONS(ZeroVelocity,
                                                 "zero velocity",
                                                 "Implementation of a model in which the boundary "
                                                 "velocity is zero. This is commonly referred to as "
                                                 "a ``stick boundary condition'', indicating that "
                                                 "the material ``sticks'' to the material on the "
                                                 "other side of the boundary.")
  }
}
