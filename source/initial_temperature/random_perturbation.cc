/*
  Copyright (C) 2017 by the authors of the ASPECT code.

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


#include <aspect/initial_temperature/random_perturbation.h>

#include <time.h>

namespace aspect
{
  namespace InitialTemperature
  {
    template <int dim>
    double
    RandomPerturbation<dim>::
    initial_temperature (const Point<dim> &/*position*/) const
    {
      // Uniform distribution on the interval [0,1]. This
      // will be used to generate the random temperature perturbation.
      boost::uniform_01<double> uniform_distribution_01;

      return magnitude * 2.0 * (uniform_distribution_01(random_number_generator) - 0.5);
    }

    template <int dim>
    void
    RandomPerturbation<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Initial temperature model");
      {
        prm.enter_subsection("Random perturbation");
        {
          prm.declare_entry ("Magnitude", "1.0",
                             Patterns::Double (0),
                             "The magnitude of the random perturbation.");
          prm.declare_entry ("Use random seed", "false",
                             Patterns::Bool(),
                             "Whether to use a random seed for the random "
                             "number generator. This parameter controls whether "
                             "this plugin generates different or identical "
                             "perturbations for different model runs.");
          prm.declare_entry ("Random seed", "5432",
                             Patterns::Integer(0),
                             "The initial seed for the random perturbation. If "
                             "'Use random seed' is set to 'false' model runs "
                             "with identical seed will lead to identical "
                             "perturbations, while model runs with different "
                             "seeds will generate different perturbations.");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }


    template <int dim>
    void
    RandomPerturbation<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Initial temperature model");
      {
        prm.enter_subsection("Random perturbation");
        {
          magnitude = prm.get_double ("Magnitude");

          if (prm.get_bool("Use random seed"))
            random_number_generator.seed(time(NULL));
          else
            random_number_generator.seed(prm.get_integer("Random seed"));
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace InitialTemperature
  {
    ASPECT_REGISTER_INITIAL_TEMPERATURE_MODEL(RandomPerturbation,
                                              "random perturbation",
                                              "An initial tempeature anomaly that perturbs the temperature "
                                              "following a random noise with uniform distribution and user "
                                              "specified magnitude.")
  }
}
