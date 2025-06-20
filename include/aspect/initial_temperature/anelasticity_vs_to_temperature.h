/*
  Copyright (C) 2016 - 2019 by the authors of the ASPECT code.

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

#ifndef _aspect_initial_temperature_anelasticity_vs_to_temperature_h
#define _aspect_initial_temperature_anelasticity_vs_to_temperature_h

#include <aspect/initial_temperature/interface.h>
#include <aspect/utilities.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace InitialTemperature
  {
    using namespace dealii;

    /**
     * A class that describes an initial temperature field for a 2D or 3D shear wave velocity (Vs) model.
     * Vs values are read from ascii data files that follow the format of the AsciiDataLookup class and
     * are converted to temperature using the anelasticity parameterization of Yamauchi & Takei (2016).
     *
     * @ingroup InitialTemperatures
     */
    template <int dim>
    class AnelasticVsToTemperature : public Utilities::AsciiDataInitial<dim>, public Interface<dim>
    {
      public:
        /**
        * Constructor. Initialize variables.
        */
        AnelasticVsToTemperature ();

        /**
        * Initialization function. This function is called once at the
        * beginning of the program. Checks preconditions.
        */
        void
        initialize ();

        // Avoid -Woverloaded-virtual:
        using Utilities::AsciiDataInitial<dim>::initialize;

        /**
        * Return the boundary temperature as a function of position. For the
        * current class, this function returns temperatures that are computed
        * from seismic velocities read in from text files.
        */
        double
        initial_temperature (const Point<dim> &position) const;

        /**
        * Declare the parameters for the input files.
        */
        static
        void
        declare_parameters (ParameterHandler  &prm);

        /**
        * Read the parameters from the parameter file.
        */
        void
        parse_parameters (ParameterHandler &prm);

      private:

        /**
        * Function that assesses difference between input and calculated Vs for Brent minimization
        */
        double
        fVs(const double x, const double y, const double z, const double a, const double b,
            const double c, const double d, const double e, const double f,
            const double g, const bool h) const;

        /**
        * Function that assesses difference between input and calculated pressure for Brent minimization
        */
        double
        fdV(const double x, const double a, const double b, const double c) const;

        /**
          * Function to calculate Vs using Yamauchi & Takei 2016 anelasticity parameterization
        */
        double
        yamauchi_takei_Vs(const double x, const double y, const double a, const double b,
                          const double c, const double d, const double e,
                          const double f, const double g, const bool h) const;

        /**
        * Whether to remove temperature heterogeneity upper parts of model
        */
        double no_perturbation_depth;

        /**
        * Whether to use original parameters published in Yamauchi & Takei (2016)
        * or a more up-to-date version that accounts for non-linear pressure
        * dependence of thermal expansivity and depth dependence of shear wave period
        */
        bool use_original_model;

        /**
        * Constant temperature to set where variations have been removed
        */
        double upper_temperature;


        /**
        * Shear modulus parameters
        */
        double mu0;
        double dmudT;
        double dmudP;

        /**
        * Viscosity parameters
        */
        double viscosity_prefactor;
        double activation_energy;
        double activation_volume;
        double solidus_gradient;

    };
  }
}

#endif
