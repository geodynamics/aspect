/*
 Copyright (C) 2016 by the authors of the ASPECT code.

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
           * @param [in,out] particles A multimap between cells and their
           * particles. This map will be filled in this function.
           */
          virtual
          void
          generate_particles(std::multimap<types::LevelInd, Particle<dim> > &particles);

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
           * Calculates particle position in the unit cell.
           *
           * @return A vector container of points each of size dim.
           */
          std::vector<Point<dim> > generate_particle_positions_in_unit_cell();

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
