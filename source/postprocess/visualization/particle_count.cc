/*
  Copyright (C) 2016 by the authors of the ASPECT code.

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


#include <aspect/postprocess/visualization/particle_count.h>

#include <aspect/postprocess/tracers.h>
#include <aspect/particle/world.h>
#include <aspect/simulator.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      std::pair<std::string, Vector<float> *>
      ParticleCount<dim>::execute() const
      {
        const Postprocess::Tracers<dim> *tracer_postprocessor = this->template find_postprocessor<Postprocess::Tracers<dim> >();

        AssertThrow(tracer_postprocessor != 0,
                    ExcMessage("The <tracers> postprocessor was not found in the list of "
                               "active postprocessors. You need to select this postprocessor to "
                               "be able to select the <particle count> visualization plugin."));

        const std::multimap<aspect::Particle::types::LevelInd, aspect::Particle::Particle<dim> > particles =
          tracer_postprocessor->get_particle_world().get_particles();

        std::pair<std::string, Vector<float> *>
        return_value ("particles_per_cell",
                      new Vector<float>(this->get_triangulation().n_active_cells()));

        // loop over all of the cells and count particles in them
        typename DoFHandler<dim>::active_cell_iterator
        cell = this->get_dof_handler().begin_active(),
        endc = this->get_dof_handler().end();

        unsigned int cell_index = 0;
        for (; cell!=endc; ++cell,++cell_index)
          if (cell->is_locally_owned())
            {
              const aspect::Particle::types::LevelInd current_cell (cell->level(),cell->index());

              (*return_value.second)(cell_index) = static_cast<float> (particles.count(current_cell));
            }

        return return_value;
      }

    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(ParticleCount,
                                                  "particle count",
                                                  "A visualization output object that generates output "
                                                  "about the number of particles per cell.")
    }
  }
}
