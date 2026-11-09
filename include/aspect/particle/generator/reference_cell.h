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

#ifndef _aspect_particle_generator_reference_cell_h
#define _aspect_particle_generator_reference_cell_h

#include <aspect/particle/generator/interface.h>

namespace aspect
{
  namespace Particle
  {
    namespace Generator
    {
      /**
       * Generate a uniform distribution of particles in the unit cell and transforms
       * each of the particles back to real region in the model domain.
       * Uniform here means the particles will be generated with
       * an equal spacing in each spatial dimension.
       *
       * @ingroup ParticleGenerators
       */
      template <int dim>
      class ReferenceCell : public Interface<dim>
      {
        public:
          /**
           * Generate a uniform distribution of particles in the unit cell and transforms
           * each of the particles back to real region in the model domain.
           * Uniform here means the particles will be generated with
           * an equal spacing in each spatial dimension.
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
           * Calculates particle position in the unit cell.
           *
           * @return A vector container of points each of size dim.
           */
          std::vector<Point<dim>> generate_particle_positions_in_unit_cell();

          /**
           * Number of particles to create for each spatial dimension as
           * specified within the parameter file.
           */
          std::vector<unsigned int> number_of_particles;
      };

    }
  }
}

#endif
