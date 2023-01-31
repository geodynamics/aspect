/*
  Copyright (C) 2016 - 2022 by the authors of the ASPECT code.

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


#include <aspect/postprocess/visualization/particle_count.h>

#include <aspect/postprocess/particles.h>
#include <aspect/particle/world.h>
#include <aspect/simulator.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      ParticleCount<dim>::
      ParticleCount ()
        :
        CellDataVectorCreator<dim>("")  // no physical units
      {}



      template <int dim>
      std::pair<std::string, Vector<float> *>
      ParticleCount<dim>::execute() const
      {
        const Postprocess::Particles<dim> &particle_postprocessor =
          this->get_postprocess_manager().template get_matching_postprocessor<Postprocess::Particles<dim>>();

        const Particle::ParticleHandler<dim> &particle_handler =
          particle_postprocessor.get_particle_world().get_particle_handler();

        std::pair<std::string, Vector<float> *>
        return_value ("particles_per_cell",
                      new Vector<float>(this->get_triangulation().n_active_cells()));

        // loop over all of the cells and count particles in them
        for (const auto &cell : this->get_dof_handler().active_cell_iterators())
          if (cell->is_locally_owned())
            {
              (*return_value.second)(cell->active_cell_index()) = static_cast<float> (particle_handler.n_particles_in_cell(cell));
            }

        return return_value;
      }



      template <int dim>
      std::list<std::string>
      ParticleCount<dim>::required_other_postprocessors() const
      {
        return std::list<std::string> (1, "particles");
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
                                                  "about the number of particles per cell."
                                                  "\n\n"
                                                  "Physical units: None.")
    }
  }
}
