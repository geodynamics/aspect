/*
  Copyright (C) 2011 - 2023 by the authors of the ASPECT code.

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


#ifndef _aspect_initial_temperature_random_gaussian_perturbation_h
#define _aspect_initial_temperature_random_gaussian_perturbation_h

#include <aspect/initial_temperature/interface.h>


namespace aspect
{
  namespace InitialTemperature
  {
    /**
     * A class that describes a perturbation to a zero temperature field
     * by placing several Gaussians with a given magnitude and width at
     * random locations throughout the model domain.
     *
     * @ingroup InitialTemperatures
     */
    template <int dim>
    class RandomGaussianPerturbation : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Initialization function. This function is called once at the
         * beginning of the program and computes the locations and
         * magnitudes of all perturbations.
         */
        void
        initialize () override;

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
         * The number of the random Gaussian perturbations.
         */
        unsigned int n_perturbations;

        /**
         * The maximal magnitude of the random Gaussian perturbations.
         */
        double max_magnitude;

        /**
         * The width of the random Gaussian perturbations.
         */
        double width;

        /**
         * Vectors that store all (random) locations and magnitudes of
         * the perturbations.
         */
        std::vector<Point<dim>> perturbation_centers;
        std::vector<double> perturbation_magnitudes;
    };
  }
}

#endif
