/*
 Copyright (C) 2017 by the authors of the ASPECT code.

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

#include <aspect/particle/interpolator/nearest_neighbor.h>
#include <aspect/postprocess/particles.h>
#include <aspect/simulator.h>

#include <deal.II/grid/grid_tools.h>

namespace aspect
{
  namespace Particle
  {
    namespace Interpolator
    {
      template <int dim>
      std::vector<std::vector<double> >
      NearestNeighbor<dim>::properties_at_points(const std::multimap<types::LevelInd, Particle<dim> > &particles,
                                                 const std::vector<Point<dim> > &positions,
                                                 const ComponentMask &selected_properties,
                                                 const typename parallel::distributed::Triangulation<dim>::active_cell_iterator &cell) const
      {
        const Postprocess::Particles<dim> *particle_postprocessor = this->template find_postprocessor<Postprocess::Particles<dim> >();

        typename parallel::distributed::Triangulation<dim>::active_cell_iterator found_cell;

        if (cell == typename parallel::distributed::Triangulation<dim>::active_cell_iterator())
          {
            // We can not simply use one of the points as input for find_active_cell_around_point
            // because for vertices of mesh cells we might end up getting ghost_cells as return value
            // instead of the local active cell. So make sure we are well in the inside of a cell.
            Assert(positions.size() > 0,
                   ExcMessage("The particle property interpolator was not given any "
                              "positions to evaluate the particle properties at."));

            const Point<dim> approximated_cell_midpoint = std::accumulate (positions.begin(), positions.end(), Point<dim>())
                                                          / static_cast<double> (positions.size());

            found_cell =
              (GridTools::find_active_cell_around_point<> (this->get_mapping(),
                                                           this->get_triangulation(),
                                                           approximated_cell_midpoint)).first;
          }
        else
          found_cell = cell;

        const types::LevelInd cell_index = std::make_pair(found_cell->level(),found_cell->index());

        const std::pair<typename std::multimap<types::LevelInd, Particle<dim> >::const_iterator,
              typename std::multimap<types::LevelInd, Particle<dim> >::const_iterator> particle_range =
                (cell->is_locally_owned())
                ?
                particles.equal_range(cell_index)
                :
                particle_postprocessor->get_particle_world().get_ghost_particles().equal_range(cell_index);

        const unsigned int n_particles = std::distance(particle_range.first,particle_range.second);
        const unsigned int n_particle_properties = particles.begin()->second.get_properties().size();

        std::vector<double> temp(n_particle_properties, 0.0);
        std::vector<std::vector<double> > point_properties(positions.size(), temp);

        for (unsigned int pos_idx=0; pos_idx < positions.size(); ++pos_idx)
          {
            double minimum_distance = std::numeric_limits<double>::max();
            if (n_particles > 0)
              {
                Particle<dim> nearest_neighbor;
                for (typename std::multimap<types::LevelInd, Particle<dim> >::const_iterator particle = particle_range.first;
                     particle != particle_range.second; ++particle)
                  {
                    const double dist = (positions[pos_idx] - particle->second.get_location()).norm_square();
                    if (dist < minimum_distance)
                      {
                        minimum_distance = dist;
                        nearest_neighbor = particle->second;
                      }
                  }
                const dealii::ArrayView<const double> neighbor_props = nearest_neighbor.get_properties();
                for (unsigned int i = 0; i < n_particle_properties; ++i)
                  if (selected_properties[i])
                    point_properties[pos_idx][i] = neighbor_props[i];
              }
            else
              {
                std::vector<typename parallel::distributed::Triangulation<dim>::active_cell_iterator> neighbors;
                GridTools::get_active_neighbors<parallel::distributed::Triangulation<dim> >(found_cell,neighbors);

                int nearest_neighbor_cell = -std::numeric_limits<int>::max();
                for (unsigned int i=0; i<neighbors.size(); ++i)
                  {
                    // Only recursively call this function if the neighbor cell contains
                    // particles (else we end up in an endless recursion)
                    if ((neighbors[i]->is_locally_owned())
                        && (particles.count(std::make_pair(neighbors[i]->level(),neighbors[i]->index())) == 0))
                      continue;
                    else if ((!neighbors[i]->is_locally_owned())
                             && (particle_postprocessor->get_particle_world().get_ghost_particles().count(
                                   std::make_pair(neighbors[i]->level(),neighbors[i]->index())) == 0))
                      continue;

                    const double dist = (positions[pos_idx] - neighbors[i]->center()).norm_square();
                    if (dist < minimum_distance)
                      {
                        minimum_distance = dist;
                        nearest_neighbor_cell = i;
                      }
                  }

                point_properties[pos_idx] = properties_at_points(particles,
                                                                 std::vector<Point<dim> > (1,positions[pos_idx]),
                                                                 selected_properties,
                                                                 neighbors[nearest_neighbor_cell])[0];
              }

            AssertThrow(minimum_distance < std::numeric_limits<double>::max(),
                        ExcMessage("Failed to find any neighbor particles."));
          }

        return point_properties;
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
      ASPECT_REGISTER_PARTICLE_INTERPOLATOR(NearestNeighbor,
                                            "nearest neighbor",
                                            "Return the properties of the nearest neighboring particle "
                                            "in the current cell, or nearest particle in nearest neighboring "
                                            "cell if current cell is empty.")
    }
  }
}
