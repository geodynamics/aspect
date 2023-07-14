/*
 Copyright (C) 2017 - 2022 by the authors of the ASPECT code.

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

#ifndef _aspect_particle_interpolator_bilinear_least_squares_h
#define _aspect_particle_interpolator_bilinear_least_squares_h

#include <aspect/particle/interpolator/interface.h>
#include <aspect/simulator_access.h>

#include <aspect/particle/interpolator/cell_average.h>

namespace aspect
{
  namespace Particle
  {
    namespace Interpolator
    {
      namespace internal
      {
        bool string_to_bool(const std::string &s);
      }

      /**
       * Evaluate the properties of all particles of the given cell
       * using a least squares projection onto the set of bilinear
       * (or, in 3d, trilinear) functions.
       *
       * @ingroup ParticleInterpolators
       */
      template <int dim>
      class BilinearLeastSquares : public Interface<dim>, public aspect::SimulatorAccess<dim>
      {
        public:
          /**
           * Return the cell-wise evaluated properties of the bilinear least squares function at the positions.
           */
          std::vector<std::vector<double>>
          properties_at_points(const ParticleHandler<dim> &particle_handler,
                               const std::vector<Point<dim>> &positions,
                               const ComponentMask &selected_properties,
                               const typename parallel::distributed::Triangulation<dim>::active_cell_iterator &cell) const override;

          // avoid -Woverloaded-virtual:
          using Interface<dim>::properties_at_points;

          /**
           * Declare the parameters this class takes through input files.
           */
          static
          void
          declare_parameters (ParameterHandler &prm);

          /**
           * Read the parameters this class declares from the parameter file.
           */
          void
          parse_parameters (ParameterHandler &prm) override;

        private:
          /**
           * A component mask that determines whether a limiting scheme is
           * used for each interpolated property. The limiting scheme
           * prevents overshoot and undershoot of interpolated particle
           * properties based on the local max and min of the particle
           * properties in that cell (i.e. the interpolated properties
           * will never exxceed the max and min of the properties on the particles).
           */
          ComponentMask use_linear_least_squares_limiter;

          /**
           * A component mask that determines whether a boundary condition
           * can be extrapolated for use in the limiting scheme. If boundary
           * extrapolation is enabled for a given property index, then the
           * limiter should be as well. Boundary extrapolation should help
           * the accuracy of properties that are smooth, although it can allow
           * undershoots and overshoots to occur if used with characteristic
           * functions or functions with discontinuities near a model boundary.
           */
          ComponentMask use_boundary_extrapolation;

          /**
           * Fallback method if there are too few particles in a cell to
           * perform a bilinear least squares interpolation.
           */
          Interpolator::CellAverage<dim> fallback_interpolator;
      };
    }
  }
}

#endif
