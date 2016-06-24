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

#include <aspect/particle/interpolator/first_particle.h>

#include <deal.II/grid/grid_tools.h>

namespace aspect
{
  namespace Particle
  {
    namespace Interpolator
    {
      /**
       * Return the properties of the first tracer of the cell containing the
       * given positions.
       */
      template <int dim>
      std::vector<std::vector<double> >
      FirstParticle<dim>::properties_at_points(const std::multimap<types::LevelInd, Particle<dim> > &particles,
                                               const std::vector<Point<dim> > &positions,
                                               const typename parallel::distributed::Triangulation<dim>::cell_iterator &cell) const
      {
        std::vector<std::vector<double> > properties(positions.size());

        // If the caller has provided a cell, we assume all positions are in that cell
        typename parallel::distributed::Triangulation<dim>::cell_iterator found_cell = cell;

        for (unsigned int point = 0; point < positions.size(); ++point)
          {
            // If the caller has not provided a cell, find the cell for every position
            if (found_cell == typename parallel::distributed::Triangulation<dim>::cell_iterator())
              {
                found_cell = (GridTools::find_active_cell_around_point<> (this->get_mapping(), this->get_triangulation(), positions[point])).first;
              }

            const types::LevelInd cell_index = std::make_pair<unsigned int, unsigned int> (cell->level(),cell->index());

            AssertThrow(particles.find(cell_index) != particles.end(),
                        ExcMessage("At least one cell contained no particles. The 'First "
                                   "particle' interpolation scheme does not support this case. "));

            // Find will only return the first particle it finds in that particular cell,
            // it is *not* the closest particle to the given position.
            const Particle<dim> particle = particles.find(cell_index)->second;

            properties[point] = particle.get_properties();
          }

        return properties;
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
      ASPECT_REGISTER_PARTICLE_INTERPOLATOR(FirstParticle,
                                            "first particle",
                                            "Return the properties of the first tracer in the given cell.")
    }
  }
}
