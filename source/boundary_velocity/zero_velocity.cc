// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2017 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#include <aspect/boundary_velocity/zero_velocity.h>


namespace aspect
{
  namespace BoundaryVelocity
  {
    template <int dim>
    Tensor<1,dim>
    ZeroVelocity<dim>::
    boundary_velocity (const types::boundary_id ,
                       const Point<dim> &) const
    {
      // return a zero tensor regardless of position
      return Tensor<1,dim>();
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace BoundaryVelocity
  {
    ASPECT_REGISTER_BOUNDARY_VELOCITY_MODEL(ZeroVelocity,
                                            "zero velocity",
                                            "Implementation of a model in which the boundary "
                                            "velocity is zero. This is commonly referred to as "
                                            "a ``stick boundary condition'', indicating that "
                                            "the material ``sticks'' to the material on the "
                                            "other side of the boundary.")
  }
}
