/*
 Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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

#ifndef __aspect__particle_generator_uniform_box_h
#define __aspect__particle_generator_uniform_box_h

#include <aspect/particle/generator/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace Particle
  {
    namespace Generator
    {
      /**
       *  Generate an uniform distribution of particles in a box region in the
       *  model domain.
       */
      template <int dim>
      class UniformBox : public Interface<dim>, public SimulatorAccess<dim>
      {
        public:
          /**
           * Constructor.
           */
          UniformBox();

          /**
           * Generate a uniformly randomly distributed set of particles in the current triangulation.
           */
          virtual
          void
          generate_particles(const double total_num_particles,
                             Particle::World<dim> &world);

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
           * The minimum coordinates of the tracer region.
           */
          Point<dim> P_min;

          /**
           * The maximum coordinates of the tracer region.
           */
          Point<dim> P_max;

          /**
           * Generate a particle at the specified position and with the
           * specified id.
           */
          void
          generate_particle(const Point<dim> &position,
                            const unsigned int id,
                            Particle::World<dim> &world);
      };

    }
  }
}

#endif
