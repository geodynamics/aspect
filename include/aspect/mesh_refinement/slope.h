/*
  Copyright (C) 2011 - 2017 by the authors of the ASPECT code.

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
        virtual
        void
        execute (Vector<float> &error_indicators) const;
    };
  }
}

#endif
