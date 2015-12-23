/*
  Copyright (C) 2015 by the authors of the ASPECT code.

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


#include <aspect/mesh_refinement/particle_density.h>

#include <aspect/postprocess/tracers.h>
#include <aspect/simulator.h>

namespace aspect
{
  namespace MeshRefinement
  {
    template <int dim>
    void
    ParticleDensity<dim>::execute(Vector<float> &indicators) const
    {
      const Postprocess::Tracers<dim> *tracer_postprocessor = this->template find_postprocessor<Postprocess::Tracers<dim> >();

      AssertThrow(tracer_postprocessor != 0,
                  ExcMessage("The mesh refinement plugin 'tracer distribution' requires the "
                             "postprocessor plugin 'tracers' to be selected. Please activate the "
                             "tracers or deactivate this mesh refinement plugin."));

      const std::multimap<Particle::types::LevelInd, Particle::Particle<dim> > *particles = &tracer_postprocessor->get_particle_world().get_particles();

      unsigned int i = 0;
      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();
      for (; cell!=endc; ++cell,++i)
        if (cell->is_locally_owned())
          {
            const Particle::types::LevelInd cell_index (cell->level(),cell->index());
            const unsigned int n_tracers = particles->count(cell_index);

            // Note that measure() only gives a linear approximation of the cell's volume.
            // This is unaccurate for curved cells, but the exact measure is likely not
            // important for this plugin since the cell distortion is usually small.
            indicators(i) = n_tracers / cell->measure();
          }
      return;
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MeshRefinement
  {
    ASPECT_REGISTER_MESH_REFINEMENT_CRITERION(ParticleDensity,
                                              "particle density",
                                              "A mesh refinement criterion that computes "
                                              "refinement indicators that equal the areal (in 2d) "
                                              "or volumetric (in 3d) density of particles in this cell. "
                                              "This plugin is useful for models with inhomogeneous "
                                              "particle density, e.g. when tracking an initial interface "
                                              "with a high particle density, or when the spatial particle "
                                              "density denotes the region of interest. Additionally, this "
                                              "plugin tends to balance the computational load between "
                                              "processes in parallel computations, because the particle "
                                              "number per cell is more similar.")
  }
}
