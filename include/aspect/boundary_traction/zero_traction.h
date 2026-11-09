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


#ifndef _aspect_boundary_traction_zero_traction_h
#define _aspect_boundary_traction_zero_traction_h

#include <aspect/boundary_traction/interface.h>

namespace aspect
{
  namespace BoundaryTraction
  {
    /**
     * A class that implements zero traction boundary conditions. This
     * is equivalent to an open boundary condition in domains where
     * hydrostatic pressure is irrelevant.
     *
     * @ingroup BoundaryTractions
     */
    template <int dim>
    class ZeroTraction : public Interface<dim>
    {
      public:
        /**
         * Return the boundary traction as a function of position. The
         * (outward) normal vector to the domain is also provided as
         * a second argument.
         *
         * For the current class, this function obviously simply returns a zero
         * tensor.
         */
        Tensor<1,dim>
        boundary_traction (const types::boundary_id boundary_indicator,
                           const Point<dim> &position,
                           const Tensor<1,dim> &normal_vector) const override;
    };
  }
}


#endif
