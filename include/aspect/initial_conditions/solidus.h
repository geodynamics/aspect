/*
  Copyright (C) 2014 by the authors of the ASPECT code.

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



#ifndef __aspect__initial_conditions_solidus_h
#define __aspect__initial_conditions_solidus_h

#include <aspect/initial_conditions/interface.h>
#include <aspect/simulator.h>

namespace aspect
{


  namespace InitialConditions
  {
    /**
     * Data class to handle the melting curve.
     */
    class MeltingCurve
    {
      public:

        /**
         * Read the data file into the class.
         */
        void read(const std::string &filename);

        /**
         * Get the melting temperature.
         */
        double T(const double p, const double radius) const;

        /**
         * Is the melting curve denpendent on radius. The melting curve can be
         * dependent on radius or pressure.
         */
        bool is_radius;

        /**
         * Number of data points in the melting curve data.
         */
        unsigned int n_points;
      private:
        /**
         * Data array for temperature.
         */
        std::vector<double> T_array;

        /**
         * Data array for pressure/radius.
         */
        std::vector<double> P_or_R_array;

        /**
         * Name of the data file.
         */
        std::string data_filename;
    };


    /**
     * A class that implements temperature initial conditions based on solidus
     * provided by a data file.
     *
     * @ingroup InitialConditionsModels
     */
    template <int dim>
    class Solidus : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */
        Solidus ();

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
         * Returns spherical coordinates of a cartesian position.
         */
        const Tensor<1,dim>
        spherical_surface_coordinates(const Tensor<1,dim> &position) const;

        /**
         * Lithosphere thickness.
         */
        double       litho_thick;

        /**
         * Magnitude of temperature perturbation.
         */
        double       magnitude_T;

        /**
         * Magnitude of lithosphere thickness perturbation.
         */
        double       magnitude_lith;

        /**
         * Temperature difference from solidus, so the initial condition can
         * be super-solidus or sub-solidus.
         */
        double       deltaT;

        /**
         * The lateral wave number  of the harmonic perturbation in the first
         * dimension. This is the only lateral wave number in 2D and equals
         * the degree of the spherical harmonics in a 3D spherical shell.
         */
        int          lateral_wave_number_1;

        /**
         * The lateral wave number of the harmonic perturbation in the second
         * dimension. This is not used in 2D and equals the order of the
         * spherical harmonics in a 3D spherical shell.
         */
        int          lateral_wave_number_2;

        /**
         * The data file name for solidus data.
         */
        std::string  solidus_filename;

        /**
         * Data class for melting curve
         */
        MeltingCurve solidus_curve;
    };
  }
}


#endif
