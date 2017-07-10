/*
 Copyright (C) 2016 - 2017 by the authors of the ASPECT code.

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

#ifndef _aspect_particle_interpolator_cell_average_h
#define _aspect_particle_interpolator_cell_average_h

#include <aspect/particle/interpolator/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace Particle
  {
    namespace Interpolator
    {
      /**
       * Return the averaged properties of all particles on the given cell.
       *
       * @ingroup ParticleInterpolators
       */
      template <int dim>
      class CellAverage : public Interface<dim>, public aspect::SimulatorAccess<dim>
      {
        public:
          /**
           * Return the cell-wise averaged properties of all particles of the cell containing the
           * given positions.
           */
          virtual
          std::vector<std::vector<double> >
          properties_at_points(const std::multimap<types::LevelInd, Particle<dim> > &particles,
                               const std::vector<Point<dim> > &positions,
                               const ComponentMask &selected_properties,
                               const typename parallel::distributed::Triangulation<dim>::active_cell_iterator &cell) const;

          // avoid -Woverloaded-virtual:
          using Interface<dim>::properties_at_points;
      };
    }
  }
}

#endif
