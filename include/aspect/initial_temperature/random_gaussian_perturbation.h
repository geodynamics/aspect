// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


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
         *
         * This variable is read from the parameter file through a parameter called 'Number of perturbations'.
         */
        unsigned int n_perturbations;

        /**
         * The maximal magnitude of the random Gaussian perturbations.
         * This variable is read from the parameter file through a parameter called 'Maximum magnitude'.
         */
        double max_magnitude;

        /**
         * The width of the random Gaussian perturbations.
         *
         * This variable is read from the parameter file through a parameter called 'Width'.
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
