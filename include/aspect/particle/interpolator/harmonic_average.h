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
 along with ASPECT; see the file LICENSE.  If not see
 <http://www.gnu.org/licenses/>.
 */

#ifndef _aspect_particle_interpolator_harmonic_average_h
#define _aspect_particle_interpolator_harmonic_average_h

#include <aspect/particle/interpolator/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace Particle
  {
    namespace Interpolator
    {
      /**
       * Return the harmonic averaged properties of all particles on the given cell.
       *
       * @ingroup ParticleInterpolators
       */
      template <int dim>
      class HarmonicAverage : public Interface<dim>, public aspect::SimulatorAccess<dim>
      {
        public:
          /**
           * Return the cell-wise harmonic averaged properties of all particles of the cell containing the
           * given positions.
           */
          virtual
          std::vector<std::vector<double> >
          properties_at_points(const ParticleHandler<dim> &particle_handler,
                               const std::vector<Point<dim> > &positions,
                               const ComponentMask &selected_properties,
                               const typename parallel::distributed::Triangulation<dim>::active_cell_iterator &cell) const;

          // avoid -Woverloaded-virtual:
          using Interface<dim>::properties_at_points;

          /**
           * @copydoc Interface<dim>::declare_parameters()
           */
          static
          void
          declare_parameters (ParameterHandler &prm);

          /**
           * @copydoc Interface<dim>::parse_parameters()
           */
          virtual
          void
          parse_parameters (ParameterHandler &prm);

        private:
          /**
           * Generally, particles need to be distributed throughout the
           * model domain when being applied to a field. This parameter
           * allows particles to be applied as a field even if some cells
           * have no particles.
           */
          bool allow_cells_without_particles;
      };
    }
  }
}

#endif
