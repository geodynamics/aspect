// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2016-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#ifndef _aspect_material_model_nondimensional_h
#define _aspect_material_model_nondimensional_h

#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace MaterialModel
  {
    /**
     * A material model for incompressible (using the Boussinesq approximation)
     * and compressible computations (with ALA or TALA) for a nondimensionalized
     * problem. The viscosity is (optionally) depth and temperature dependent.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class Nondimensional : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:

        /**
         * Initialize.
         */
        void initialize () override;

        /**
         * Evaluate the material model.
         */
        void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                      MaterialModel::MaterialModelOutputs<dim> &out) const override;

        /**
         * @name Qualitative properties one can ask a material model
         * @{
         */

        /**
         * Return whether the model is compressible or not.
         */
        bool is_compressible () const override;
        /**
         * @}
         */

        /**
         * @name Functions used in dealing with run-time parameters
         * @{
         */
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
        /**
         * @}
         */

      private:
        /**
         * Use truncated anelastic approximation?
         *
         * This variable is read from the parameter file through a parameter called 'Use TALA'.
         */
        bool tala;

        /**
         * Are we using ALA/TALA or BA? (determined by the choice of Di)
         */
        bool compressible;

        /**
         * Parameter describing the temperature prefactor of viscosity.
         * This variable is read from the parameter file through a parameter called 'Viscosity temperature prefactor'.
         */
        double exponential_viscosity_temperature_prefactor;

        /**
         * Parameter describing the depth prefactor of viscosity.
         * This variable is read from the parameter file through a parameter called 'Viscosity depth prefactor'.
         */
        double exponential_viscosity_depth_prefactor;

        /**
         * The surface density.
         * This variable is read from the parameter file through a parameter called 'Reference density'.
         */
        double reference_rho;

        /**
         * The nondimensional numbers (Dissipation number,
         * Rayleigh number, Grueneisen parameter).
         * This variable is read from the parameter file through parameters called
         * 'Di', 'Ra', and 'gamma'.
         */
        double Di, Ra, gamma;

        /**
         * The constant specific heat
         *
         * This variable is read from the parameter file through a parameter called 'Reference specific heat'.
         */
        double reference_specific_heat;

    };

  }
}

#endif
