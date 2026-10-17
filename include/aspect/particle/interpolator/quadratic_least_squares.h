// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2019-2024 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

#ifndef _aspect_particle_interpolator_quadratic_least_squares_h
#define _aspect_particle_interpolator_quadratic_least_squares_h

#include <aspect/particle/interpolator/interface.h>
#include <aspect/particle/interpolator/cell_average.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace Particle
  {
    namespace Interpolator
    {
      /**
       * Return the properties of all particles of the given cell by using a least squares
       * projection onto the space of quadratic functions.
       * For more details on the implementation and properties of this plugin see:
       * Mack Gregory and Elbridge Gerry Puckett (2020),
       * "Improving the Accuracy and Efficiency of Hybrid Finite
       * Element / Particle-In-Cell Methods for Modeling Geologic Processes",
       * Abstract Submission ID# 743654
       * https://agu.confex.com/agu/fm20/prelim.cgi/Paper/743654
       * Session: A002 - Addressing Challenges for the Next Generation of Earth System Models
       * Submitted to 2020 AGU Fall Meeting, San Francisco, CA Dec. 7-11, 2020
       *
       * @ingroup ParticleInterpolators
       */
      template <int dim>
      class QuadraticLeastSquares : public Interface<dim>, public aspect::SimulatorAccess<dim>
      {
        public:
          /**
           * Return the cell-wise evaluated properties of the quadratic least squares function at the positions.
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
           * Use the cell average interpolator in case the quadratic least
           * squares interpolator fails due to a lack of particles.
           */
          Interpolator::CellAverage<dim> fallback_interpolator;

          /**
           * Enables a limiting scheme that prevents overshoot and
           * undershoot of interpolated particles based upon the range
           * of property values found in the area of the cell.
           */
          ComponentMask use_quadratic_least_squares_limiter;

          /**
           * Enables a linear extrapolation of boundary values for the limiting scheme
           */
          ComponentMask use_boundary_extrapolation;

          /**
           * Calculate the value of the interpolation function at a given position
           */
          double evaluate_interpolation_function(const Vector<double> &coefficients, const Point<dim> &position) const;

          /**
           * Update the bounds of where the plane reaches by checking whether each of the critical points
           * are in the cell and evauating their value.
           */
          std::pair<double, double> get_interpolation_bounds(const dealii::Vector<double> &coefficients) const;

          /*
           * Find all points that may contain the minimum or maximum values of the interpolation in the cell.
           */
          std::vector<dealii::Point<dim>> get_critical_points(const dealii::Vector<double> &coefficients) const;

      };
    }
  }
}

#endif
