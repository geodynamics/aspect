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

namespace aspect
{
  namespace InitialComposition
  {
    template <int dim>
    RandomPerturbation<dim>::
    RandomPerturbation ()
      :
      random_number_generator(5432)
    {
    }

    template <int dim>
    double
    RandomPerturbation<dim>::
    initial_composition (const Point<dim> &/*position*/,
                         const unsigned int /*compositional_index*/) const
    {
      // Uniform distribution on the interval [0,1]. This
      // will be used to generate the random composition perturbation.
      boost::uniform_01<double> uniform_distribution_01;

      return magnitude * 2.0 * (uniform_distribution_01(random_number_generator) - 0.5);
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
