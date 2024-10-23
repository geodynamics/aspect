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

#include <aspect/particle/generator/probability_density_function.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/geometry_info.h>

namespace aspect
{
  namespace Particle
  {
    namespace Generator
    {
      template <int dim>
      void
      ProbabilityDensityFunction<dim>::generate_particles(Particles::ParticleHandler<dim> &particle_handler)
      {
        Particles::Generators::probabilistic_locations(this->get_triangulation(),
                                                       function,
                                                       random_cell_selection,
                                                       n_particles,
                                                       particle_handler,
                                                       this->get_mapping(),
                                                       random_number_seed);
      }



      template <int dim>
      void
      ProbabilityDensityFunction<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Generator");
        {
          prm.enter_subsection("Probability density function");
          {
            Functions::ParsedFunction<dim>::declare_parameters (prm, 1);

            // This parameter overwrites one of the parameters in ParsedFunction
            // with a new default value of 1.0. The original default value was 0
            // everywhere, which is not allowed.
            prm.declare_entry("Function expression",
                              "1.0",
                              Patterns::Anything(),
                              "The formula that denotes the spatially variable probability density function. "
                              "This expression "
                              "may contain any of the usual operations such as addition or "
                              "multiplication, as well as all of the common functions such as "
                              "`sin' or `cos'. In addition, it may contain expressions like "
                              "`if(x>0, 1, 0)' where the expression evaluates to the second "
                              "argument if the first argument is true, and to the third argument "
                              "otherwise; this example would result in no particles at all in that "
                              "part of the domain where $x==0$, and a constant particle density "
                              "in the rest of the domain. For a full overview of possible expressions accepted "
                              "see the documentation of the muparser library at http://muparser.beltoforion.de/. "
                              "Note that the function has to be non-negative everywhere in the domain, "
                              "and needs to be positive in at least some parts of the domain.");

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
      ProbabilityDensityFunction<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Generator");
        {
          prm.enter_subsection("Probability density function");
          {
            n_particles = static_cast<types::particle_index>(prm.get_double ("Number of particles"));
            random_cell_selection = prm.get_bool("Random cell selection");
            random_number_seed = prm.get_integer("Random number seed");

            try
              {
                function.parse_parameters (prm);
              }
            catch (...)
              {
                std::cerr << "ERROR: FunctionParser failed to parse\n"
                          << "\t'Particle.Generator.Probability density function'\n"
                          << "with expression\n"
                          << "\t'" << prm.get("Function expression") << "'";
                throw;
              }
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
      ASPECT_REGISTER_PARTICLE_GENERATOR(ProbabilityDensityFunction,
                                         "probability density function",
                                         "Generate a random distribution of "
                                         "particles over the entire simulation domain. "
                                         "The probability density is prescribed in the "
                                         "form of a user-prescribed function. The "
                                         "format of this function follows the syntax "
                                         "understood by the muparser library, see "
                                         "{ref}\\`sec:run-aspect:parameters-overview:muparser-format\\`. The "
                                         "return value of the function is always "
                                         "checked to be a non-negative probability "
                                         "density but it can be zero in "
                                         "parts of the domain.")
    }
  }
}
