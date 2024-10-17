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


#ifndef _aspect_initial_temperature_random_perturbation_h
#define _aspect_initial_temperature_random_perturbation_h

#include <aspect/initial_temperature/interface.h>
#include <aspect/simulator_access.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#include <random>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

namespace aspect
{
  namespace InitialTemperature
  {
    using namespace dealii;

    /**
     * A class that describes a perturbation for the initial temperature field
     * for any geometry model or dimension. The perturbation follows
     * a random noise of a prescribed magnitude.
     *
     * @ingroup InitialTemperatures
     */
    template <int dim>
    class RandomPerturbation : public Interface<dim>, public aspect::SimulatorAccess<dim>
    {
      public:
        void initialize () override;
        /**
         * Return the initial temperature as a function of position.
         */

        double initial_temperature (const Point<dim> &position) const override;

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
         * The maximal magnitude of the random noise.
         */
        double magnitude;

        /**
         * Whether to use a random seed for the random
         * number generator. This parameter controls whether
         * this plugin generates different or identical
         * perturbations for subsequent model runs of
         * the same setup.
         */
        bool use_random_seed;

        /**
         * A random number generator used by this class to get
         * random temperature perturbations. Should return the same
         * random numbers every time since it is seeded with one.
         */
        mutable std::mt19937 random_number_generator;

    };
  }
}

#endif
