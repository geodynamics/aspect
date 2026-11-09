// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2019 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#ifndef _aspect_mesh_refinement_slope_h
#define _aspect_mesh_refinement_slope_h

#include <aspect/mesh_refinement/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace MeshRefinement
  {

    /**
     * A class that implements a mesh refinement criterion that refines the
     * mesh with a deforming surface.  Specifically, it calculates the slope of the
     * surface by comparing the local normal and gravity vectors.
     * Cells with steeper slopes get refined. This is useful for cases where
     * there is steep topography which needs to be tracked better.
     *
     * @ingroup MeshRefinement
     */
    template <int dim>
    class Slope : public Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:

        /**
         * Execute this mesh refinement criterion.
         *
         * @param[out] error_indicators A vector that for every active cell of
         * the current mesh (which may be a partition of a distributed mesh)
         * provides an error indicator. This vector will already have the
         * correct size when the function is called.
         */
        void
        execute (Vector<float> &error_indicators) const override;
    };
  }
}

#endif
