// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2025-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#ifndef _aspect_heating_model_tidal_heating_h
#define _aspect_heating_model_tidal_heating_h

#include <aspect/heating_model/interface.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace HeatingModel
  {
    /**
     * A class that implements tidal heating.
     *
     * @ingroup HeatingModels
     */
    template <int dim>
    class TidalHeating : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Initialize function, which sets the start time and
         * start timestep of tidal heating.
         */
        void initialize() override;

      public:
        /**
         * Return the tidal heating terms.
         */
        void
        evaluate (const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
                  const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
                  HeatingModel::HeatingModelOutputs &heating_model_outputs) const override;

        /**
         * Specify which material model outputs the heating model requires
         * for computing the heating terms.
         */
        MaterialModel::MaterialProperties::Property
        get_required_properties () const override;

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
         * Parameters used for tidal heating (H), which is defined using the following
         * Equation is from Tobie et al. (2003) (https://doi.org/10.1029/2003JE002099)
         * H = 2*(viscosity)*(time-averaged tidal strain rate)^2/(1+((viscosity)*(tidal frequency)/(elastic shear modulus))^2))
         * viscosity (Pa s) = viscosity calculated by the selected material in ASPECT
         * time-averaged strain rate = constant_tidal_strain_rate
         * tidal frequency = tidal_frequency
         * elastic shear modulus = elastic_shear_modulus
         *
         * strain_rate_distribution lets selection for distribution of tidal strain rate.
         * If 'constant', the tidal strain rate is fixed to 'Constant tidal strain rate'.
         * If 'latitudinal variation', 'Maximum tidal strain rate' and 'Minimum tidal strain rate' are used.
         * Then, the equation follows as (maximum_tidal_strain_rate - minimum_tidal_strain_rate)*cos(theta/2)+(maximum_tidal_strain_rate+minimum_tidal_strain_rate)/2.
         * This is shown in Fig.3 of Nimmo et al. (2007) (https://doi.org/10.1016/j.icarus.2007.04.021).
         *
         * This variable is read from the parameter file through the 'Tidal
         * frequency' parameter.
         */
        double tidal_frequency;
        /**
         * This variable is read from the parameter file through a parameter called 'Elastic shear modulus'.
         */
        double elastic_shear_modulus;
        /**
         * This variable is read from the parameter file through a parameter called 'Constant tidal strain rate'.
         */
        double constant_tidal_strain_rate;

        /**
         * Specify the distribution of time-averaged tidal strain rate.
         */
        enum StrainRateDistribution
        {
          constant,
          latitudinal_variation
        } strain_rate_distribution;

        /**
         * This variable is read from the parameter file through a parameter called 'Maximum tidal strain rate'.
         */
        double maximum_tidal_strain_rate;
        /**
         * This variable is read from the parameter file through a parameter called 'Minimum tidal strain rate'.
         */
        double minimum_tidal_strain_rate;
    };
  }
}


#endif
