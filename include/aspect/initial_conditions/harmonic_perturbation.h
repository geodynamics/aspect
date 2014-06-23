/*
  Copyright (C) 2011, 2012 by the authors of the ASPECT code.

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


#ifndef __aspect__initial_conditions_harmonic_perturbation_h
#define __aspect__initial_conditions_harmonic_perturbation_h

#include <aspect/initial_conditions/interface.h>
#include <aspect/simulator.h>


namespace aspect
{
  namespace InitialConditions
  {
    using namespace dealii;

    /**
     * A class that describes a perturbed initially constant temperature field
     * for any geometry model or dimension in shape of a harmonic function.
     * For 3D spherical shell models this is achieved by using spherical
     * harmonics, in any other case sine function are scaled to fit the model
     * geometry.
     *
     * @ingroup InitialConditionsModels
     */
    template <int dim>
    class HarmonicPerturbation : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Return the initial temperature as a function of position.
         */
        virtual
        double initial_temperature (const Point<dim> &position) const;

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);


      private:

        /**
         * The radial/depth wave number of the harmonic perturbation. All wave
         * number variables are in fact twice the wave number in a
         * mathematical sense. This allows the user to prescribe a single up-
         * / downswing or half periods.
         */
        int vertical_wave_number;

        /**
         * The lateral wave number  of the harmonic perturbation in the first
         * dimension. This is the only lateral wave number in 2D and equals
         * the degree of the spherical harmonics in a 3D spherical shell.
         */
        int lateral_wave_number_1;

        /**
         * The lateral wave number of the harmonic perturbation in the second
         * dimension. This is not used in 2D and equals the order of the
         * spherical harmonics in a 3D spherical shell.
         */
        int lateral_wave_number_2;

        /**
         * The maximal magnitude of the harmonic perturbation.
         */
        double magnitude;

        /**
         * The background temperature the harmonic perturbation is applied on
         * in an incompressible material model. In case of a compressible
         * material model the perturbation is applied on top of an adiabatic
         * profile and this variable is not used at all.
         */
        double reference_temperature;
    };
  }
}

#endif
