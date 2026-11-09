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

#ifndef _aspect_particle_generator_uniform_radial_h
#define _aspect_particle_generator_uniform_radial_h

#include <aspect/particle/generator/interface.h>

#include <array>

namespace aspect
{
  namespace Particle
  {
    namespace Generator
    {
      /**
       * Generate a uniform radial distribution of particles over the entire
       * simulation domain. Uniform here means
       * the particles will be generated with an equal spacing in
       * each spherical spatial dimension, i.e. the particles are
       * created at coordinates that increase linearly with equal
       * spacing in radius, colatitude and longitude around a
       * certain center point. Note that in order
       * to produce a regular distribution the number of generated
       * particles might not exactly match the one specified in the
       * input file.
       *
       * @ingroup ParticleGenerators
       */
      template <int dim>
      class UniformRadial : public Interface<dim>
      {
        public:
          /**
           * Generate a uniformly distributed set of particles in a
           * circular or spherical subdomain of the global domain.
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
           * The minimum spherical coordinates of the particle region, i.e.
           * the first radius, colatitude and longitude from the given
           * center position P_center where particles are generated.
           *
           * This variable is read from the parameter file through parameters
           * called 'Minimum radius', 'Minimum longitude', and 'Minimum latitude'.
           */
          std::array<double,dim> P_min;

          /**
           * The maximum spherical coordinates of the particle region, i.e.
           * the last radius, colatitude and longitude from the given
           * center position P_center where particles are generated.
           * This variable is read from the parameter file through parameters
           * called 'Maximum radius', 'Maximum longitude', and 'Maximum latitude'.
           */
          std::array<double,dim> P_max;

          /**
           * The center of the particle region. Defaults to the origin.
           * This variable is read from the parameter file through parameters
           * called 'Center x', 'Center y', and 'Center z'.
           */
          Point<dim> P_center;

          /**
           * The number of radial layers of particles that will be generated.
           * In particular this parameter determines the radial spacing between
           * particle layers as Pmax[0] - P_min[0] / radial_layers.
           * This variable is read from the parameter file through a parameter
           * called 'Radial layers'.
           */
          unsigned int radial_layers;

          /**
           * Number of particles to create in total. To preserve the equal
           * spacing (in spherical coordinates) between particles, the number
           * of actually generated
           * particles can differ slightly from this number.
           *
           * This variable is read from the parameter file through a parameter
           * called 'Number of particles'.
           */
          types::particle_index n_particles;
      };

    }
  }
}

#endif
