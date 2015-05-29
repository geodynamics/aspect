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

#ifndef __aspect__particle_generator_uniform_radial_h
#define __aspect__particle_generator_uniform_radial_h

#include <aspect/particle/generator/interface.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/std_cxx11/array.h>

namespace aspect
{
  namespace Particle
  {
    namespace Generator
    {
      using namespace dealii;

      /**
       * Generate uniform radial distribution of particles over entire simulation domain.
       */
      template <int dim>
      class UniformRadial : public Interface<dim>, public SimulatorAccess<dim>
      {
        public:
          /**
           * Constructor.
           */
          UniformRadial();

          /**
           * Generate a uniformly distributed set of particles in the current circular or
           * spherical domain
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
           * The minimum spherical coordinates of the tracer region.
           */
          std_cxx11::array<double,dim> P_min;

          /**
           * The maximum spherical coordinates of the tracer region.
           */
          std_cxx11::array<double,dim> P_max;

          /**
           * The center of the tracer region. Defaults to the origin.
           */
          Point<dim> P_center;

          /**
           * The number of radial layers of particles that will be generated.
           */
          unsigned int radial_layers;
      };

    }
  }
}

#endif
