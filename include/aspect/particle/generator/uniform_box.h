// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2015-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

#ifndef _aspect_particle_generator_uniform_box_h
#define _aspect_particle_generator_uniform_box_h

#include <aspect/particle/generator/interface.h>

namespace aspect
{
  namespace Particle
  {
    namespace Generator
    {
      /**
       * Generate a uniform distribution of particles in a box region in the
       * model domain. Uniform here means the particles will be generated with
       * an equal spacing in each spatial dimension. Note that in order
       * to produce a regular distribution the number of generated
       * particles might not exactly match the one specified in the
       * input file.
       *
       * @ingroup ParticleGenerators
       */
      template <int dim>
      class UniformBox : public Interface<dim>
      {
        public:
          /**
           * Generate a uniformly randomly distributed set of particles in a
           * box-like subdomain of the global domain.
           *
           * @param [in,out] particle_handler The particle handler into which
           * the generated particles should be inserted.
           */
          void
          generate_particles(Particles::ParticleHandler<dim> &particle_handler) override;

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
           * Number of initial particles to create.
           *
           * This variable is read from the parameter file through a parameter
           * called 'Number of particles'.
           */
          types::particle_index n_particles;

          /**
           * The minimum coordinates of the particle region, i.e. one corner of
           * the n-dimensional box region in which particles are generated.
           * This variable is read from the parameter file through parameters
           * called 'Minimum x', 'Minimum y', and 'Minimum z'.
           */
          Point<dim> P_min;

          /**
           * The maximum coordinates of the particle region, i.e. the opposite
           * corner of the n-dimensional box region in which particles are
           * generated.
           *
           * This variable is read from the parameter file through parameters
           * called 'Maximum x', 'Maximum y', and 'Maximum z'.
           */
          Point<dim> P_max;
      };

    }
  }
}

#endif
