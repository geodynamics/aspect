// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2016-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

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
          std::vector<std::vector<double>>
          properties_at_points(const ParticleHandler<dim> &particle_handler,
                               const std::vector<Point<dim>> &positions,
                               const ComponentMask &selected_properties,
                               const typename parallel::distributed::Triangulation<dim>::active_cell_iterator &cell) const override;

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
          void
          parse_parameters (ParameterHandler &prm) override;

        private:
          /**
           * By default, every cell needs to contain particles to use this interpolator
           * plugin. If this parameter is set to true, cells are allowed to have no particles,
           * in which case the interpolator will return 0 for the cell's properties.
           *
           * This variable is read from the parameter file through a parameter
           * called 'Allow cells without particles'.
           */
          bool allow_cells_without_particles;
      };
    }
  }
}

#endif
