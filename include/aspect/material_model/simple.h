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

#ifndef _aspect_material_model_simple_h
#define _aspect_material_model_simple_h

#include <aspect/material_model/interface.h>
#include <aspect/material_model/equation_of_state/linearized_incompressible.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace MaterialModel
  {
    /**
     * A material model that consists of globally constant values for all
     * material parameters except density and viscosity.
     *
     * The model is considered incompressible, following the definition
     * described in Interface::is_compressible. This is essentially the
     * material model used in the step-32 tutorial program.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class Simple : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                      MaterialModel::MaterialModelOutputs<dim> &out) const override;

        /**
         * @name Qualitative properties one can ask a material model
         * @{
         */

        /**
         * Return whether the model is compressible or not.  Incompressibility
         * does not necessarily imply that the density is constant; rather, it
         * may still depend on temperature or pressure. In the current
         * context, compressibility means whether we should solve the continuity
         * equation as $\nabla \cdot (\rho \mathbf u)=0$ (compressible Stokes)
         * or as $\nabla \cdot \mathbf{u}=0$ (incompressible Stokes).
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
         * This variable is read from the parameter file through a parameter called 'Reference temperature'.
         */
        double reference_T;
        /**
         * This variable is read from the parameter file through a parameter called 'Viscosity'.
         */
        double eta;
        /**
         * This variable is read from the parameter file through a parameter called 'Composition viscosity prefactor'.
         */
        double composition_viscosity_prefactor;
        /**
         * This variable is read from the parameter file through a parameter called 'Thermal viscosity exponent'.
         */
        double thermal_viscosity_exponent;
        /**
         * This variable is read from the parameter file through a parameter called 'Maximum thermal prefactor'.
         */
        double maximum_thermal_prefactor;
        /**
         * This variable is read from the parameter file through a parameter called 'Minimum thermal prefactor'.
         */
        double minimum_thermal_prefactor;

        /**
         * The thermal conductivity.
         *
         * This variable is read from the parameter file through a parameter called 'Thermal conductivity'.
         */
        double k_value;

        EquationOfState::LinearizedIncompressible<dim> equation_of_state;
    };

  }
}

#endif
