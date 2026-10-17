// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2015-2024 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

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
