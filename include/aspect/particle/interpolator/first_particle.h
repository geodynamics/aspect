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

#ifndef __aspect__particle_interpolator_first_h
#define __aspect__particle_interpolator_first_h

#include <aspect/particle/interpolator/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace Particle
  {
    namespace Interpolator
    {
      /**
       * Return the properties of the first tracer of the given cell.
       *
       * @ingroup ParticleInterpolators
       */
      template <int dim>
      class FirstParticle : public Interface<dim>, public aspect::SimulatorAccess<dim>
      {
        public:
          /**
           * Perform an interpolation of the properties of the particles in
           * this cell onto a vector of positions in this cell.
           * Implementations of this function must return a vector of a vector
           * of doubles which contains a somehow computed
           * value of all particle properties at all given positions.
           *
           * @param [in] particles Reference to the particle map.
           * @param [in] positions The vector of positions where the properties
           * should be evaluated.
           * @param [in] cell An optional iterator to the cell containing the
           * particles. Not all callers will know the cell of the particles,
           * but providing the cell when known speeds up the interpolation
           * significantly.
           * @return A vector with as many entries as @p positions. Every entry
           * is a vector of interpolated tracer properties at this position.
           */
          virtual
          std::vector<std::vector<double> >
          properties_at_points(const std::multimap<types::LevelInd, Particle<dim> > &particles,
                               const std::vector<Point<dim> > &positions,
                               const typename parallel::distributed::Triangulation<dim>::cell_iterator &cell = typename parallel::distributed::Triangulation<dim>::cell_iterator()) const;
      };

    }
  }
}

#endif
