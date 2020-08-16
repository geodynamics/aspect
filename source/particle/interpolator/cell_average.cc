/*
  Copyright (C) 2016 - 2019 by the authors of the ASPECT code.

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

#include <aspect/particle/interpolator/cell_average.h>
#include <aspect/postprocess/particles.h>
#include <aspect/simulator.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/base/signaling_nan.h>

namespace aspect
{
  namespace Particle
  {
    namespace Interpolator
    {
      template <int dim>
      std::vector<std::vector<double> >
      CellAverage<dim>::properties_at_points(const ParticleHandler<dim> &particle_handler,
                                             const std::vector<Point<dim> > &positions,
                                             const ComponentMask &selected_properties,
                                             const typename parallel::distributed::Triangulation<dim>::active_cell_iterator &cell) const
      {
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

        const typename ParticleHandler<dim>::particle_iterator_range particle_range =
          particle_handler.particles_in_cell(found_cell);

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
                    cell_properties[i] += particle_properties[i];
              }

            for (unsigned int i = 0; i < n_particle_properties; ++i)
              if (selected_properties[i])
                cell_properties[i] /= n_particles;
          }
        // If there are no particles in this cell use the average of the
        // neighboring cells.
        else
          {
            std::vector<typename parallel::distributed::Triangulation<dim>::active_cell_iterator> neighbors;
            GridTools::get_active_neighbors<parallel::distributed::Triangulation<dim> >(found_cell,neighbors);

            unsigned int non_empty_neighbors = 0;
            for (unsigned int i=0; i<neighbors.size(); ++i)
              {
                // Only recursively call this function if the neighbor cell contains
                // particles (else we end up in an endless recursion)
                if (particle_handler.n_particles_in_cell(neighbors[i]) == 0)
                  continue;

                const std::vector<double> neighbor_properties = properties_at_points(particle_handler,
                                                                                     std::vector<Point<dim> > (1,neighbors[i]->center(true,false)),
                                                                                     selected_properties,
                                                                                     neighbors[i])[0];

                for (unsigned int i = 0; i < n_particle_properties; ++i)
                  if (selected_properties[i])
                    cell_properties[i] += neighbor_properties[i];

                ++non_empty_neighbors;
              }

            if (!allow_cells_without_particles)
              {
                AssertThrow(non_empty_neighbors != 0,
                            ExcMessage("A cell and all of its neighbors do not contain any particles. "
                                       "The `cell average' interpolation scheme does not support this case unless specified "
                                       "in Allow cells without particles."));
              }

            for (unsigned int i = 0; i < n_particle_properties; ++i)
              {
                if (selected_properties[i] && non_empty_neighbors != 0)
                  cell_properties[i] /= non_empty_neighbors;
                // Assume property is zero for any areas with no particles
                else if (allow_cells_without_particles && non_empty_neighbors == 0)
                  cell_properties[i] = 0;
              }

          }

        return std::vector<std::vector<double> > (positions.size(),cell_properties);
      }



      template <int dim>
      void
      CellAverage<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Particles");
          {
            prm.declare_entry ("Allow cells without particles", "false",
                               Patterns::Bool (),
                               "By default, every cell needs to contain particles to use this interpolator "
                               "plugin. If this parameter is set to true, cells are allowed to have no particles, "
                               "in which case the interpolator will return 0 for the cell's properties.");
          }
          prm.leave_subsection ();
        }
        prm.leave_subsection ();
      }



      template <int dim>
      void
      CellAverage<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Particles");
          {
            allow_cells_without_particles = prm.get_bool("Allow cells without particles");
          }
          prm.leave_subsection ();
        }
        prm.leave_subsection ();
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
      ASPECT_REGISTER_PARTICLE_INTERPOLATOR(CellAverage,
                                            "cell average",
                                            "Return the average of all particle properties in the given cell.")
    }
  }
}
