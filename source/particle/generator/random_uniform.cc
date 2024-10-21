/*
  Copyright (C) 2015 - 2024 by the authors of the ASPECT code.

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

#include <aspect/particle/generator/random_uniform.h>

#include <deal.II/base/function.h>

namespace aspect
{
  namespace Particle
  {
    namespace Generator
    {
      template <int dim>
      void
      RandomUniform<dim>::generate_particles(Particles::ParticleHandler<dim> &particle_handler)
      {
        Particles::Generators::probabilistic_locations(this->get_triangulation(),
                                                       Functions::ConstantFunction<dim>(1.0),
                                                       random_cell_selection,
                                                       n_particles,
                                                       particle_handler,
                                                       this->get_mapping(),
                                                       random_number_seed);
      }



      template <int dim>
      void
      RandomUniform<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Generator");
        {
          prm.enter_subsection("Random uniform");
          {
            prm.declare_entry ("Number of particles", "1000",
                               Patterns::Double (0.),
                               "Total number of particles to create (not per processor or per element). "
                               "The number is parsed as a floating point number (so that one can "
                               "specify, for example, '1e4' particles) but it is interpreted as "
                               "an integer, of course.");

            prm.declare_entry ("Random cell selection", "true",
                               Patterns::Bool(),
                               "If true, particle numbers per cell are calculated randomly "
                               "according to their respective probability density. "
                               "This means particle numbers per cell can deviate statistically from "
                               "the integral of the probability density. If false, "
                               "first determine how many particles each cell should have "
                               "based on the integral of the density over each of the cells, "
                               "and then once we know how many particles we want on each cell, "
                               "choose their locations randomly within each cell.");

            prm.declare_entry ("Random number seed", "5432",
                               Patterns::Integer(0),
                               "The seed for the random number generator that controls "
                               "the particle generation. Keep constant to generate "
                               "identical particle distributions in subsequent model "
                               "runs. Change to get a different distribution. In parallel "
                               "computations the seed is further modified on each process "
                               "to ensure different particle patterns on different "
                               "processes. Note that the number of particles per processor "
                               "is not affected by the seed.");
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }


      template <int dim>
      void
      RandomUniform<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Generator");
        {
          prm.enter_subsection("Random uniform");
          {
            n_particles = static_cast<types::particle_index>(prm.get_double ("Number of particles"));
            random_cell_selection = prm.get_bool("Random cell selection");
            random_number_seed = prm.get_integer("Random number seed");
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Particle
  {
    namespace Generator
    {
      ASPECT_REGISTER_PARTICLE_GENERATOR(RandomUniform,
                                         "random uniform",
                                         "Generates a random uniform distribution of "
                                         "particles over the entire simulation domain. "
                                         "This generator can be understood as the special case "
                                         "of the 'probability density function' generator where "
                                         "the probability density is constant over the domain.")
    }
  }
}
