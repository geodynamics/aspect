/*
 Copyright (C) 2016 - 2022 by the authors of the ASPECT code.

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

#ifndef _aspect_particle_generator_quadrature_points_h
#define _aspect_particle_generator_quadrature_points_h

#include <aspect/particle/generator/interface.h>

namespace aspect
{
  namespace Particle
  {
    namespace Generator
    {
      /**
       * Generates particles at the quadrature points of each active cell of the triangulation.
       * Here, Gauss quadrature of degree, (velocity_degree + 1), is used similarly to the assembly of the
       * Stokes matrix.
       * @ingroup ParticleGenerators
       */
      template <int dim>
      class QuadraturePoints : public Interface<dim>
      {
        public:
          /**
           * Generates particles at the quadrature points of each active cell of the triangulation.
           * Here, Gauss quadrature of degree, (velocity_degree + 1), is used similarly to the assembly of the
           * Stokes matrix.
           *
           * @param [in,out] particle_handler The particle handler into which
           * the generated particles should be inserted.
           */
          void
          generate_particles(Particles::ParticleHandler<dim> &particle_handler) override;
      };
    }
  }
}

#endif
