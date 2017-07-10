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


#ifndef _aspect_boundary_traction_zero_traction_h
#define _aspect_boundary_traction_zero_traction_h

#include <aspect/boundary_traction/interface.h>

namespace aspect
{
  namespace BoundaryTraction
  {
    using namespace dealii;

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
        virtual
        Tensor<1,dim>
        boundary_traction (const types::boundary_id boundary_indicator,
                           const Point<dim> &position,
                           const Tensor<1,dim> &normal_vector) const;
    };
  }
}


#endif
