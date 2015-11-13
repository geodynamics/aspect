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


#ifndef __aspect__traction_boundary_conditions_zero_traction_h
#define __aspect__traction_boundary_conditions_zero_traction_h

#include <aspect/traction_boundary_conditions/interface.h>

namespace aspect
{
  namespace TractionBoundaryConditions
  {
    using namespace dealii;

    /**
     * A class that implements zero traction boundary conditions. This
     * is equivalent to an open boundary condition in domains where
     * hydrostatic pressure is irrelevant.
     *
     * @ingroup TractionBoundaryConditionsModels
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
        traction (const Point<dim> &position,
                  const Tensor<1,dim> &normal_vector) const;
    };
  }
}


#endif
