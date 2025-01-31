/*
  Copyright (C) 2013 - 2022 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_latent_heat_h
#define _aspect_material_model_latent_heat_h

#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace MaterialModel
  {
    /**
     * A material model that implements a standard approximation of the latent
     * heat terms following Christensen \& Yuen, 1986. The change of entropy
     * is calculated as $Delta S = \gamma \frac{\Delta\rho}{\rho^2}$ with the
     * Clapeyron slope $\gamma$ and the density change $\Delta\rho$ of the
     * phase transition being input parameters. This model employs an analytic
     * phase function in the form $X=\frac{1}{2} \left( 1 + \tanh \left( \frac{\Delta
     * p}{\Delta p_0} \right) \right)$ with $\Delta p = p - p_transition -
     * \gamma \left( T - T_transition \right)$ and $\Delta p_0$ being the
     * pressure difference over the width of the phase transition (specified
     * as input parameter).
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class LatentHeat : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:

        /**
         * Evaluate material properties.
         */
        void evaluate(const MaterialModelInputs<dim> &in,
                      MaterialModelOutputs<dim> &out) const override;

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
        double reference_rho;
        double reference_T;
        double eta;
        double composition_viscosity_prefactor;
        double thermal_viscosity_exponent;
        double thermal_alpha;
        double reference_specific_heat;
        double reference_compressibility;
        double maximum_viscosity;
        double minimum_viscosity;

        /**
         * The thermal conductivity.
         */
        double k_value;

        double compositional_delta_rho;

        // list of depth (or pressure), width and Clapeyron slopes
        // for the different phase transitions
        std::vector<double> density_jumps;
        std::vector<int> transition_phases;
        std::vector<double> phase_prefactors;

        MaterialUtilities::PhaseFunction<dim> phase_function;
    };

  }
}

#endif
