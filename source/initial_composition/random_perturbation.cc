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


#include <aspect/initial_composition/random_perturbation.h>

// #include <boost/functional/hash.hpp>
#include <random>
#include <time.h>

namespace aspect
{
  namespace InitialComposition
  {
    template <int dim>
    void
    RandomPerturbation<dim>::
    initialize()
    {
      const unsigned int my_rank = Utilities::MPI::this_mpi_process(this->get_mpi_communicator());
      if (use_random_seed)
        random_number_generator.seed(time(NULL)+my_rank);
      else
        random_number_generator.seed(9088+my_rank);
    }

    template <int dim>
    double
    RandomPerturbation<dim>::
    initial_composition (const Point<dim> &position,  const unsigned int composition_index) const
    {
      // Uniform distribution on the interval [-magnitude,magnitude). This
      // will be used to generate the random temperature perturbation.
      std::uniform_real_distribution<double> uniform_distribution(-magnitude,magnitude);
      return uniform_distribution(random_number_generator);
    }

    template <int dim>
    void
    RandomPerturbation<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Initial composition model");
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
                             "perturbations for subsequent model runs of"
                             "the same setup.");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }


    template <int dim>
    void
    RandomPerturbation<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Initial composition model");
      {
        prm.enter_subsection("Random perturbation");
        {
          magnitude = prm.get_double ("Magnitude");
          use_random_seed = prm.get_bool("Use random seed");
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
  namespace InitialComposition
  {
    ASPECT_REGISTER_INITIAL_COMPOSITION_MODEL(RandomPerturbation,
                                              "random perturbation",
                                              "An initial composition anomaly that perturbs the composition "
                                              "following a random noise with uniform distribution and user "
                                              "specified magnitude.")
  }
}
