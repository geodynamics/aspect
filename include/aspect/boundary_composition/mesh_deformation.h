/*
  Copyright (C) 2013 - 2023 by the authors of the ASPECT code.

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


#ifndef _aspect_boundary_composition_mesh_deformation_h
#define _aspect_boundary_composition_mesh_deformation_h

#include <aspect/mesh_deformation/interface.h>
#include <aspect/boundary_composition/interface.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace BoundaryComposition
  {
    /**
     * A class that implements a composition boundary condition that
     * is set by the active mesh deformation plugins. Their returned
     * values at a given point on a boundary are summed.
     *
     * @ingroup BoundaryCompositions
     */
    template <int dim>
    class MeshDeformation : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Initialization function. This function is called once at the
         * beginning of the program after parse_parameters is run.
         */
        void
        initialize () override;

        /**
         * This function returns the composition value returned by the
         * mesh deformation interface, which is the sum of the boundary compositions
         * provided by the active mesh deformation objects.
         */
        double boundary_composition (const types::boundary_id boundary_indicator,
                                     const Point<dim> &position,
                                     const unsigned int compositional_field) const override;
    };
  }
}


#endif
