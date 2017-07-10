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

#ifndef _aspect_particle_interpolator_bilinear_least_squares_h
#define _aspect_particle_interpolator_bilinear_least_squares_h

#include <aspect/particle/interpolator/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace Particle
  {
    namespace Interpolator
    {
      /**
       * Return the interpolated properties of all particles of the given cell using bilinear least squares method.
       * Currently, only the two dimensional model is supported.
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
          virtual
          std::vector<std::vector<double> >
          properties_at_points(const std::multimap<types::LevelInd, Particle<dim> > &particles,
                               const std::vector<Point<dim> > &positions,
                               const ComponentMask &selected_properties,
                               const typename parallel::distributed::Triangulation<dim>::active_cell_iterator &cell) const;

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
          virtual
          void
          parse_parameters (ParameterHandler &prm);

        private:
          /**
           * Variables related to a limiting scheme that prevents overshoot and
           * undershoot of interpolated particle properties based on global max
           * and global min for each propery.
           */
          bool use_global_valued_limiter;

          /**
           * For each interpolated particle property, a global max and global
           * min are stored as elements of vectors.
           */
          std::vector<double> global_maximum_particle_properties;
          std::vector<double> global_minimum_particle_properties;
      };
    }
  }
}

#endif
