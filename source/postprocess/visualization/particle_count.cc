/*
  Copyright (C) 2016 - 2024 by the authors of the ASPECT code.

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
#include <aspect/particle/manager.h>


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
      std::pair<std::string, std::unique_ptr<Vector<float>>>
      ParticleCount<dim>::execute() const
      {
        std::pair<std::string, std::unique_ptr<Vector<float>>>
        return_value ("particles_per_cell",
                      std::make_unique<Vector<float>>(this->get_triangulation().n_active_cells()));

        // loop over all of the cells and count particles in them
        for (const auto &cell : this->get_dof_handler().active_cell_iterators())
          if (cell->is_locally_owned())
            {
              unsigned int n_particles_in_cell = 0;
              for (unsigned int particle_manager_index = 0;
                   particle_manager_index < this->n_particle_managers();
                   ++particle_manager_index)
                n_particles_in_cell += this->get_particle_manager(particle_manager_index).get_particle_handler().n_particles_in_cell(cell);

              (*return_value.second)(cell->active_cell_index()) = static_cast<float>(n_particles_in_cell);
            }

        return return_value;
      }



      template <int dim>
      std::list<std::string>
      ParticleCount<dim>::required_other_postprocessors() const
      {
        return {"particles"};
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
