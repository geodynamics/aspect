/*
 Copyright (C) 2015 - 2022 by the authors of the ASPECT code.

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
           */
          std::array<double,dim> P_min;

          /**
           * The maximum spherical coordinates of the particle region, i.e.
           * the last radius, colatitude and longitude from the given
           * center position P_center where particles are generated.
           */
          std::array<double,dim> P_max;

          /**
           * The center of the particle region. Defaults to the origin.
           */
          Point<dim> P_center;

          /**
           * The number of radial layers of particles that will be generated.
           * In particular this parameter determines the radial spacing between
           * particle layers as Pmax[0] - P_min[0] / radial_layers.
           */
          unsigned int radial_layers;

          /**
           * Number of particles to create in total. To preserve the equal
           * spacing (in spherical coordinates) between particles, the number
           * of actually generated
           * particles can differ slightly from this number.
           */
          types::particle_index n_particles;
      };

    }
  }
}

#endif
