// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2021 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


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

      if (this->get_parameters().include_melt_transport)
        {
          value(dim) = 0;       // fluid pressure
          value(dim+1) = 0;     // compaction pressure

          value(dim+2) = -p(1); // fluid velocity x
          value(dim+3) = p(0);  // fluid velocity y
          if (dim == 3)
            value(dim+4) = 0;

          value(2*dim+2) = 0;   // pressure
        }
      else
        {
          value(dim) = 0;       // pressure
        }
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
