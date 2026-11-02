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


#include <aspect/boundary_traction/zero_traction.h>


namespace aspect
{
  namespace BoundaryTraction
  {
    template <int dim>
    Tensor<1,dim>
    ZeroTraction<dim>::
    boundary_traction (const types::boundary_id,
                       const Point<dim> &,
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
  namespace BoundaryTraction
  {
    ASPECT_REGISTER_BOUNDARY_TRACTION_MODEL(ZeroTraction,
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
