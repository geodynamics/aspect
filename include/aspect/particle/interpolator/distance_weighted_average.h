// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2024 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

#ifndef _aspect_particle_interpolator_distance_weighted_average_h
#define _aspect_particle_interpolator_distance_weighted_average_h

#include <aspect/particle/interpolator/interface.h>
#include <aspect/simulator_access.h>

#include <deal.II/grid/grid_tools_cache.h>

namespace aspect
{
  namespace Particle
  {
    namespace Interpolator
    {
      /**
       * Return the interpolated properties of all particles of the given cell using
       * a distance weighted averaging. The weight function is a hat function.
       *
       * @ingroup ParticleInterpolators
       */
      template <int dim>
      class DistanceWeightedAverage : public Interface<dim>, public aspect::SimulatorAccess<dim>
      {
        public:
          /**
           * Initialize function. Called once at the beginning of the model.
           */
          void initialize() override;

          /**
           * Return the cell-wise evaluated particle properties at the given @p positions.
           */
          std::vector<std::vector<double>>
          properties_at_points(const ParticleHandler<dim> &particle_handler,
                               const std::vector<Point<dim>> &positions,
                               const ComponentMask &selected_properties,
                               const typename parallel::distributed::Triangulation<dim>::active_cell_iterator &cell) const override;

        private:
          /**
           * Cached information that stores information about the grid so that we
           * do not need to recompute it every time properties_at_points() is called.
           */
          std::unique_ptr<GridTools::Cache<dim>> grid_cache;
      };
    }
  }
}

#endif
