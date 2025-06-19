/*
  Copyright (C) 2011 - 2025 by the authors of the ASPECT code.

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

#ifndef _aspect_cookbooks_phase_transition_kinetics_phase_transition_kinetics_h
#define _aspect_cookbooks_phase_transition_kinetics_phase_transition_kinetics_h

#include <deal.II/base/patterns.h>
#include <deal.II/base/types.h>

#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/material_model/equation_of_state/interface.h>
#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>

#include <vector>

namespace aspect
{
  namespace MaterialModel
  {
    /**
     * A material model that reads in thermodynamic data from an ascii .txt file.
     * The thermodynamic data are evaluated along an adiabatic reference profile and
     * used to compute a reaction rate between two phases (or fixed mineral assemblages)
     * using the operator splitting technique. Average material properties are computed from
     * the mass fractions of the reacting phases. The model is considered compressible.
     *
     * The viscosity is computed as
     * \f[
     * \eta(z,T) = \eta_r(z) \eta_0 \exp\left(-A \frac{T - \bar{T}}{\bar{T}}\right)."
     * \f]
     *
     * where $\eta_r(z)$ is a depth-dependent viscosity prefactor, $\eta_0$ is the
     * reference viscosity, $A$ is the thermal viscosity exponent, $T$ is the full
     * temperature, and $\bar{T}$ is the reference adiabatic temperature.
     *
     * @ingroup MaterialModels
     */

    template <int dim>
    class PhaseTransitionKinetics : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor. Initialize variables.
         */
        PhaseTransitionKinetics();

        /**
         * Initialization function. This function is called once at the
         * beginning of the program after parse_parameters is run and after
         * the SimulatorAccess (if applicable) is initialized.
         */
        void
        initialize() override;

        /**
         * @name Qualitative properties one can ask a material model
         * @{
         */

        /**
         * Return whether the model is compressible or not.
         */
        bool
        is_compressible() const override;

        /**
         * @}
         */

        /**
         * Evaluate material properties.
         */
        void
        evaluate(const MaterialModel::MaterialModelInputs<dim> &in, MaterialModel::MaterialModelOutputs<dim> &out) const override;

        /**
         * @name Functions used in dealing with run-time parameters
         * @{
         */

        /**
         * Declare the parameters this class takes through input files.
         */
        static void
        declare_parameters(ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        void
        parse_parameters(ParameterHandler &prm) override;

        /**
         * @}
         */

        /**
         * Add the named outputs for reaction rates.
         */
        void
        create_additional_named_outputs(MaterialModel::MaterialModelOutputs<dim> &out) const override;

      private:
        /**
         * Object that stores the thermodynamic data used for computing the equation
         * of state and phase reaction. The thermodynamic data are evaluated along a
         * reference adiabatic profile.
         *
         * Data is read from a tab-separated .txt file with the following required columns:
         *
         * pressure density_a density_b thermal_expansivity_a thermal_expansivity_b
         * specific_heat_a specific_heat_b compressibility_a compressibility_b
         * delta_molar_gibbs delta_molar_entropy delta_molar_volume
         *
         * Note: "a" and "b" represent phases "a" and "b".
         */
        Utilities::AsciiDataProfile<dim> profile;

        /**
         * Column indices from the thermodynamic data along an adiabatic profile.
         * The pressure column is used for searching the table via interpolation.
         * All other columns store thermodynamic data required for calculating
         * the thermodynamic driving force in the PhaseTransitionKinetics material model.
         *
         * Note: "a" and "b" represent phases "a" and "b" and "dG", "dS", "dV"
         * represent the differences in the molar Gibbs free energy, molar entropy,
         * and molar volume of phases "a" and "b". For example, $\text{dG}$ =
         * $\text{G}_b$ - $\text{G}_a$
         */
        unsigned int rho_a_idx;
        unsigned int rho_b_idx;
        unsigned int alpha_a_idx;
        unsigned int alpha_b_idx;
        unsigned int beta_a_idx;
        unsigned int beta_b_idx;
        unsigned int cp_a_idx;
        unsigned int cp_b_idx;
        unsigned int dG_idx;
        unsigned int dS_idx;
        unsigned int dV_idx;

        /**
         * The reference viscosity.
         *
         * Units: Pa s
         */
        double viscosity;

        /**
         * The minimum viscosity.
         *
         * Units: Pa s
         */
        double minimum_viscosity;

        /**
         * The maximum viscosity.
         *
         * Units: Pa s
         */

        double maximum_viscosity;

        /**
         * The constant $A$ in the temperature dependence of viscosity
         * $\exp\left(-A \frac{T - \bar{T}}{\bar{T}}\right).$
         *
         * Units: none
         */
        double thermal_viscosity_exponent;

        /**
         * A list of depths that determine the locations of the jumps in
         * the piece-wise constant function $\eta_r(z)$, which describes the
         * depth dependence of viscosity.
         *
         * Units: m
         */
        std::vector<double> transition_depths;

        /**
         * A list of constants that make up the piece-wise constant function
         * $\eta_r(z)$, which determines the depth dependence of viscosity,
         * and is multiplied with the reference viscosity and the
         * temperature dependence to compute the viscosity $\eta(z,T)$.
         *
         * Units: none
         */
        std::vector<double> viscosity_prefactors;

        /**
         * The reference thermal conductivity
         *
         * Units: W/m/k
         */
        double k;

        /**
         * The kinetic prefactor Q used to calculate the reaction rate
         * $\frac{dX}{dt} = Q \Delta G (1 - X)$
         *
         * Units: J/mol/s
         */
        double Q_kinetic_prefactor;
    };
  } // namespace MaterialModel
} // namespace aspect

#endif
