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


#ifndef _aspect_initial_temperature_harmonic_perturbation_h
#define _aspect_initial_temperature_harmonic_perturbation_h

#include <aspect/initial_temperature/interface.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace InitialTemperature
  {
    /**
     * A class that describes a perturbed initially constant temperature field
     * for any geometry model or dimension in shape of a harmonic function.
     * For 3D spherical shell models this is achieved by using spherical
     * harmonics, in any other case sine function are scaled to fit the model
     * geometry.
     *
     * @ingroup InitialTemperatures
     */
    template <int dim>
    class HarmonicPerturbation : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
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
         * The radial/depth wave number of the harmonic perturbation. All wave
         * number variables are in fact twice the wave number in a
         * mathematical sense. This allows the user to prescribe a single
         * up-/downswing or half periods.
         *
         * This variable is read from the parameter file through a parameter called 'Vertical wave number'.
         */
        int vertical_wave_number;

        /**
         * The lateral wave number  of the harmonic perturbation in the first
         * dimension. This is the only lateral wave number in 2D and equals
         * the degree of the spherical harmonics in a 3D spherical shell.
         * This variable is read from the parameter file through a parameter called 'Lateral wave number one'.
         */
        int lateral_wave_number_1;

        /**
         * The lateral wave number of the harmonic perturbation in the second
         * dimension. This is not used in 2D and equals the order of the
         * spherical harmonics in a 3D spherical shell.
         * This variable is read from the parameter file through a parameter called 'Lateral wave number two'.
         */
        int lateral_wave_number_2;

        /**
         * The maximal magnitude of the harmonic perturbation.
         * This variable is read from the parameter file through a parameter called 'Magnitude'.
         */
        double magnitude;

        /**
         * The background temperature the harmonic perturbation is applied on
         * in an incompressible material model. In case of a compressible
         * material model the perturbation is applied on top of an adiabatic
         * profile and this variable is not used at all.
         *
         * This variable is read from the parameter file through a parameter called 'Reference temperature'.
         */
        double reference_temperature;
    };
  }
}

#endif
