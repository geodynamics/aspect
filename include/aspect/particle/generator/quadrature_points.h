// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2016-2024 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

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
