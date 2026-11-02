// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2025 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#ifndef _aspect_boundary_velocity_zero_velocity_h
#define _aspect_boundary_velocity_zero_velocity_h

#include <aspect/boundary_velocity/interface.h>

namespace aspect
{
  namespace BoundaryVelocity
  {
    /**
     * A class that implements zero velocity (stick) boundary conditions.
     *
     * @ingroup BoundaryVelocities
     */
    template <int dim>
    class ZeroVelocity : public Interface<dim>
    {
      public:
        /**
         * Return the boundary velocity as a function of position. For the
         * current class, this function obviously simply returns a zero
         * tensor.
         */
        Tensor<1,dim>
        boundary_velocity (const types::boundary_id boundary_indicator,
                           const Point<dim> &position) const override;
    };
  }
}


#endif
