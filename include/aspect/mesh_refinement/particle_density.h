/*
  Copyright (C) 2015 - 2016 by the authors of the ASPECT code.

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


#ifndef _aspect_mesh_refinement_particle_density_h
#define _aspect_mesh_refinement_particle_density_h

#include <aspect/mesh_refinement/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace MeshRefinement
  {

    /**
     * A mesh refinement criterion that computes
     * refinement indicators based on the density
     * of particles. In practice this plugin
     * equilibrates the number of particles per cell,
     * leading to fine cells in high particle density regions
     * and coarse cells in low particle density regions.
     * This plugin is mostly useful for models with inhomogeneous
     * particle density, e.g. when tracking an initial interface
     * with a high particle density, or when the spatial particle
     * density denotes the region of interest. Additionally, this
     * plugin tends to balance the computational load between
     * processes in parallel computations, because the particle
     * and mesh density is more aligned.
     *
     * @ingroup MeshRefinement
     */
    template <int dim>
    class ParticleDensity : public Interface<dim>,
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
