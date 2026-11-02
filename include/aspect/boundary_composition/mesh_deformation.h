// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2013-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


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
