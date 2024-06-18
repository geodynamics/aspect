/*
  Copyright (C) 2024 by the authors of the ASPECT code.

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

#include <aspect/particle/interpolator/distance_weighted_average.h>

#include <deal.II/base/signaling_nan.h>
#include <deal.II/grid/grid_tools.h>

namespace aspect
{
  namespace Particle
  {
    namespace Interpolator
    {
      namespace internal
      {
        double weight (const double distance, const double interpolation_range)
        {
          // linear weight function (hat function)
          return std::max(1.0 - (distance/interpolation_range),0.0);
        }
      }



      template <int dim>
      void
      DistanceWeightedAverage<dim>::initialize()
      {
        grid_cache = std::make_unique<GridTools::Cache<dim>>(this->get_triangulation(), this->get_mapping());
      }



      template <int dim>
      std::vector<std::vector<double>>
      DistanceWeightedAverage<dim>::properties_at_points(const ParticleHandler<dim> &particle_handler,
                                                         const std::vector<Point<dim>> &positions,
                                                         const ComponentMask &selected_properties,
                                                         const typename parallel::distributed::Triangulation<dim>::active_cell_iterator &cell) const
      {
        const unsigned int n_interpolate_positions = positions.size();
        const unsigned int n_particle_properties = particle_handler.n_properties_per_particle();

        // Create with signaling NaNs
        std::vector<std::vector<double>> cell_properties(n_interpolate_positions,
                                                          std::vector<double>(n_particle_properties,
                                                                              numbers::signaling_nan<double>()));

        // Set requested properties to 0.0
        for (unsigned int index_positions = 0; index_positions < n_interpolate_positions; ++index_positions)
          for (unsigned int index_properties = 0; index_properties < n_particle_properties; ++index_properties)
            if (selected_properties[index_properties])
              cell_properties[index_positions][index_properties] = 0.0;

        std::set<typename Triangulation<dim>::active_cell_iterator> cell_and_neighbors;

        const auto &vertex_to_cell_map = grid_cache->get_vertex_to_cell_map();

        for (const auto v : cell->vertex_indices())
          {
            const unsigned int vertex_index = cell->vertex_index(v);
            cell_and_neighbors.insert(vertex_to_cell_map[vertex_index].begin(),
                                      vertex_to_cell_map[vertex_index].end());
          }

        // Average over all particles that are within half a cell diameter.
        // This distance strikes a balance between required number of particles per
        // cell and accuracy of the interpolation. This also assumes we can find
        // most particles within this distance in neighbor cells. If we do not find
        // some particles in range (e.g. because the neighbor is refined) this
        // will slightly affect the accuracy of the interpolation, but we accept that.
        //
        // TODO: This could be made dependent on the number of particles per cell, the more
        // particles, the smaller the interpolation range to increase accuracy.
        const double interpolation_range = 0.5 * cell->diameter();

        std::vector<double> integrated_weight(n_interpolate_positions,0.0);

        for (const auto &current_cell: cell_and_neighbors)
          {
            const typename ParticleHandler<dim>::particle_iterator_range particle_range =
              particle_handler.particles_in_cell(current_cell);

            for (const auto &particle: particle_range)
              {
                const ArrayView<const double> particle_properties = particle.get_properties();
                unsigned int index_positions = 0;

                for (const auto &interpolation_point: positions)
                  {
                    const double distance = particle.get_location().distance(interpolation_point);
                    const double weight = internal::weight(distance, interpolation_range);

                    for (unsigned int index_properties = 0; index_properties < particle_properties.size(); ++index_properties)
                      if (selected_properties[index_properties])
                        cell_properties[index_positions][index_properties] += weight * particle_properties[index_properties];

                    integrated_weight[index_positions] += weight;
                    ++index_positions;
                  }
              }
          }

        for (unsigned int index_positions = 0; index_positions < n_interpolate_positions; ++index_positions)
          {
            AssertThrow(integrated_weight[index_positions] > 0.0, ExcInternalError());

            for (unsigned int index_properties = 0; index_properties < n_particle_properties; ++index_properties)
              if (selected_properties[index_properties])
                cell_properties[index_positions][index_properties] /= integrated_weight[index_positions];
          }

        return cell_properties;
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
      ASPECT_REGISTER_PARTICLE_INTERPOLATOR(DistanceWeightedAverage,
                                            "distance weighted average",
                                            "Interpolates particle properties onto a vector of points using a "
                                            "distance weighed averaging method.")
    }
  }
}
