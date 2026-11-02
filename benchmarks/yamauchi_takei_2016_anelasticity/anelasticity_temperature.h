// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2016-2025 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#ifndef _aspect_initial_temperature_anelasticity_temperature_h
#define _aspect_initial_temperature_anelasticity_temperature_h

#include <aspect/initial_temperature/interface.h>
#include <aspect/utilities.h>
#include <aspect/simulator_access.h>
#include <cmath>
#include <algorithm>
#include <functional>

namespace aspect
{
  namespace InitialTemperature
  {
    /**
     * A class that describes an initial temperature field for a 2D or 3D shear wave velocity (Vs) model.
     * Vs values are converted to temperature using the anelasticity parameterization of Yamauchi & Takei (2016).
     *
     * @ingroup InitialTemperatures
     */
    template <int dim>
    class AnelasticVs2T : public Utilities::AsciiDataInitial<dim>, public Interface<dim>
    {
      public:
        /**
         * Constructor. Initialize variables.
         */
        AnelasticVs2T ();

        /**
         * Initialization function. This function is called once at the
         * beginning of the program. Checks preconditions.
         */
        void
        initialize () override;

        // Avoid -Woverloaded-virtual:
        using Utilities::AsciiDataInitial<dim>::initialize;

        /**
         * Return the boundary temperature as a function of position. For the
         * current class, this function returns value from the text files.
         */
        double
        initial_temperature (const Point<dim> &position) const override;

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
        parse_parameters (ParameterHandler &prm) override;

      private:

        /**
         * Function that assesses difference between input and calculated Vs for Brent minimization
         */
        double
        fVs(double x, double depth, double absolute_Vs, double mu0, double dmudT,
            double dmudP, double viscosity_prefactor, double activation_energy, double activation_volume,
            double solidus_gradient, bool density_model_flag) const;

        /**
         * Function that assesses difference between input and calculated pressure for Brent minimization
         */
        double
        fdV(double x, double bulk_modulus, double bulk_modulus_pressure_derivative, double pressure) const;

        /**
         * Function to calculate Vs using Yamauchi & Takei 2016 anelasticity parameterization
         */
        double
        yamauchi_takei_Vs(double temperature, double depth, double mu0, double dmudT,
                          double dmudP, double viscosity_prefactor, double activation_energy,
                          double activation_volume, double solidus_gradient, bool density_model_flag) const;

        /**
         * Whether to remove temperature heterogeneity upper parts of model
         */
        double no_perturbation_depth;

        /**
         * Constant temperature to set where variations have been removed
         */
        double reference_temperature;

        /**
         * Whether to use Yamauchi & Takei (2016) anelasticity
         * parameterization
         */
        bool use_yamauchi_takei;

        /**
         * Whether to use original parameters published in Yamauchi & Takei
         * (2016) or an updated version that accounts for non-linear pressure
         * dependence of thermal expansivity and depth dependence of shear
         * wave period
         */
        bool use_original_model;

    };
  }
}

#endif
