/*
  Copyright (C) 2017 - 2024 by the authors of the ASPECT code.

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

#include <aspect/particle/interpolator/harmonic_average.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/base/signaling_nan.h>

namespace aspect
{
  namespace Particle
  {
    namespace Interpolator
    {
      template <int dim>
      std::vector<std::vector<double>>
      HarmonicAverage<dim>::properties_at_points(const ParticleHandler<dim> &particle_handler,
                                                 const std::vector<Point<dim>> &positions,
                                                 const ComponentMask &selected_properties,
                                                 const typename parallel::distributed::Triangulation<dim>::active_cell_iterator &cell) const
      {
        const typename ParticleHandler<dim>::particle_iterator_range particle_range =
          particle_handler.particles_in_cell(cell);

        const unsigned int n_particles = std::distance(particle_range.begin(),particle_range.end());
        const unsigned int n_particle_properties = particle_handler.n_properties_per_particle();

        std::vector<double> cell_properties (n_particle_properties,numbers::signaling_nan<double>());

        for (unsigned int i = 0; i < n_particle_properties; ++i)
          if (selected_properties[i])
            cell_properties[i] = 0.0;

        if (n_particles > 0)
          {
            for (const auto &particle : particle_range)
              {
                const ArrayView<const double> &particle_properties = particle.get_properties();

                for (unsigned int i = 0; i < particle_properties.size(); ++i)
                  if (selected_properties[i])
                    {
                      AssertThrow(particle_properties[i] > 0,
                                  ExcMessage ("All particle property values must be greater than 0 for harmonic averaging!"));
                      cell_properties[i] += 1/particle_properties[i];
                    }
              }

            for (unsigned int i = 0; i < n_particle_properties; ++i)
              if (selected_properties[i])
                {
                  cell_properties[i] = n_particles/cell_properties[i];
                }
          }
        // If there are no particles in this cell use the average of the
        // neighboring cells.
        else
          {
            std::vector<typename parallel::distributed::Triangulation<dim>::active_cell_iterator> neighbors;
            GridTools::get_active_neighbors<parallel::distributed::Triangulation<dim>>(cell,neighbors);

            unsigned int non_empty_neighbors = 0;
            for (unsigned int i=0; i<neighbors.size(); ++i)
              {
                // Only recursively call this function if the neighbor cell contains
                // particles (else we end up in an endless recursion)
                if (particle_handler.n_particles_in_cell(neighbors[i]) == 0)
                  continue;

                const std::vector<double> neighbor_properties = properties_at_points(particle_handler,
                                                                                     std::vector<Point<dim>> (1,neighbors[i]->center(true,false)),
                                                                                     selected_properties,
                                                                                     neighbors[i])[0];

                for (unsigned int i = 0; i < n_particle_properties; ++i)
                  if (selected_properties[i])
                    {
                      AssertThrow(neighbor_properties[i] > 0,
                                  ExcMessage ("All particle property values must be greater than 0 for harmonic averaging!"));
                      cell_properties[i] += 1/neighbor_properties[i];
                    }

                ++non_empty_neighbors;
              }

            if (!allow_cells_without_particles)
              {
                AssertThrow(non_empty_neighbors != 0,
                            ExcMessage("A cell and all of its neighbors do not contain any particles. "
                                       "The `harmonic average' interpolation scheme does not support this case unless specified "
                                       "in Allow cells without particles."));
              }

            for (unsigned int i = 0; i < n_particle_properties; ++i)
              {
                if (selected_properties[i] && non_empty_neighbors !=0)
                  cell_properties[i] = non_empty_neighbors/cell_properties[i];
                // Assume property is zero for any areas with no particles
                else if (allow_cells_without_particles && non_empty_neighbors == 0)
                  cell_properties[i] = 0;
              }

          }

        return std::vector<std::vector<double>> (positions.size(),cell_properties);
      }



      template <int dim>
      void
      HarmonicAverage<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.declare_entry ("Allow cells without particles", "false",
                           Patterns::Bool (),
                           "By default, every cell needs to contain particles to use this interpolator "
                           "plugin. If this parameter is set to true, cells are allowed to have no particles. "
                           "In case both the current cell and its neighbors are empty, "
                           "the interpolator will return 0 for the current cell's properties.");
      }



      template <int dim>
      void
      HarmonicAverage<dim>::parse_parameters (ParameterHandler &prm)
      {
        allow_cells_without_particles = prm.get_bool("Allow cells without particles");
      }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Particle
  {
    namespace Interpolator
    {
      ASPECT_REGISTER_PARTICLE_INTERPOLATOR(HarmonicAverage,
                                            "harmonic average",
                                            "Return the harmonic average of all particle properties in the "
                                            "given cell. If the cell contains no particles, return the "
                                            "harmonic average of the properties in the neighboring cells. "
                                            "In case the neighboring cells are also empty, and 'Allow cells "
                                            "without particles' is set to true, the interpolator returns 0. "
                                            "Otherwise, an exception is thrown. ")
    }
  }
}
