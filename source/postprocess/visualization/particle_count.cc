// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2016-2024 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


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
