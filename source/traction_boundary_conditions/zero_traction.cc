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


#include <aspect/traction_boundary_conditions/zero_traction.h>


namespace aspect
{
  namespace TractionBoundaryConditions
  {
    template <int dim>
    Tensor<1,dim>
    ZeroTraction<dim>::
    traction (const Point<dim> &,
              const Tensor<1,dim> &) const
    {
      // return a zero tensor regardless of position
      return Tensor<1,dim>();
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace TractionBoundaryConditions
  {
    ASPECT_REGISTER_TRACTION_BOUNDARY_CONDITIONS(ZeroTraction,
                                                 "zero traction",
                                                 "Implementation of a model in which the boundary "
                                                 "traction is zero. This is commonly referred to as "
                                                 "an ``open boundary condition'', indicating that "
                                                 "the material experiences no forces in response to "
                                                 "what might exist on the other side of the boundary. "
                                                 "However, this is only true in the case where "
                                                 "hydrostatic pressure is not relevant. If hydrostatic "
                                                 "pressure is not negligible, for example at the sides "
                                                 "of a regional model, the material at the other side "
                                                 "of the boundary does exceed a force, namely the "
                                                 "force normal to the boundary induced by the "
                                                 "hydrostatic pressure.")
  }
}
