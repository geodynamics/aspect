/*
  Copyright (C) 2011 - 2016 by the authors of the ASPECT code.

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


#ifndef _aspect_initial_temperature_harmonic_perturbation_citcoms_h
#define _aspect_initial_temperature_harmonic_perturbation_citcoms_h

#include <aspect/initial_temperature/interface.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace InitialTemperature
  {
    using namespace dealii;
    /**
     * This setup is a different implementation than the one in the main code
     * and it is based on the Arrial et al. (2014) setup. It generates an
     * initial constant temperature field which is perturbed following a
     * spherical harmonic function in lateral and radial direction. This
     * setup can only be used for a hollow sphere.
     *
     * @ingroup InitialTemperatures
     */
    template <int dim>
    class HarmonicPerturbationCitcomS : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
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
         * number variables are in fact twice the wave number in a mathematical
         * sense. This allows the user to prescribe a single up-/ downswing or
         * half periods.
         */
        int vertical_wave_number;

        /**
         * The degree of the spherical harmonics in a 3D spherical shell
         * for the lateral wave number of the axisymmetric harmonic
         * perturbation.
         */
        int lateral_wave_number_l1;

        /**
         * The order of the spherical harmonics in a 3D spherical shell
         * for the lateral wave number of the axysymmetric harmonic
         * perturbation.
         */
        int lateral_wave_number_m1;

        /**
         * The degree of the spherical harmonics in a 3D spherical shell
         * for the lateral wave number of the nonaxisymmetric harmonic
         * perturbation.
         */
        int lateral_wave_number_l2;

        /**
         * The order of the spherical harmonics in a 3D spherical shell
         * for the lateral wave number of the nonaxisymmetric harmonic
         * perturbation.
         */
        int lateral_wave_number_m2;

        /**
         * The maximal magnitude of the spherical harmonic perturbation.
         */
        double magnitude;

        /**
         * A perturbation parameter has been introduced to slowly pertrub the
         * amplitude of the nonaxisymmetric mode..
         */
        double delta;
    };
  }
}

#endif
