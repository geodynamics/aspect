/*
  Copyright (C) 2011 - 2026 by the authors of the ASPECT code.

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

#ifndef _aspect_cookbooks_mtz_kinetics_mtz_kinetics_h
#define _aspect_cookbooks_mtz_kinetics_mtz_kinetics_h

#include <deal.II/base/patterns.h>
#include <deal.II/base/types.h>

#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/material_model/reaction_model/reaction_chain.h>
#include <aspect/material_model/equation_of_state/interface.h>
#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/introspection.h>
#include <aspect/utilities.h>

#include <vector>

namespace aspect
{
  namespace MaterialModel
  {
    /**
     * A compressible material model that implements a chain of phase transformations goverened by reaction kinetics.
     *
     * Each transformation step reads thermodynamic reference data from an ASCII .tsv file. The thermodynamic data are evaluated along an adiabatic
     * reference profile and used to compute reaction rates. Material properties are computed from the mass fractions of the reacting phases.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class MTZKinetics : public MaterialModel::Interface<dim>, public virtual ::aspect::SimulatorAccess<dim>
    {
      public:
        MTZKinetics();

        void initialize() override;

        bool is_compressible() const override;

        void evaluate(const MaterialModel::MaterialModelInputs<dim> &in, MaterialModel::MaterialModelOutputs<dim> &out) const override;

        void clamp_reaction_progress_fields(MaterialModel::MaterialModelInputs<dim> &in) const override;

        double correct_gibbs(const double dG_ref, const double dS, const double dV, const double P, const double P_ref, const double T, const double T_ref) const;

        static void declare_parameters(ParameterHandler &prm);

        void parse_parameters(ParameterHandler &prm) override;

        void create_additional_named_outputs(MaterialModel::MaterialModelOutputs<dim> &out) const override;

      private:
        /**
         * Reaction chain manager configuring kinetic models for all sequential phase transitions in the system.
         */
        ReactionModel::ReactionChain<dim> reaction_chain;

        /**
         * Array of ASCII data profiles storing the thermodynamic and physical data for each reaction.
         */
        std::vector<Utilities::AsciiDataProfile<dim>> profiles;

        /**
         * Data directory and file names for the ASCII profiles.
         */
        std::string data_directory;
        std::vector<std::string> data_file_names;

        /**
         * Per-phase material property column indices.
         */
        std::vector<unsigned int> rho_idx, alpha_idx, beta_idx, cp_idx;
        std::vector<unsigned int> Vp_idx, Vs_idx, dVp_dT_idx, dVs_dT_idx;

        /**
         * Per-reaction thermodynamic property column indices.
         */
        std::vector<unsigned int> dG_idx, dS_idx, dV_idx;

        /**
         * Cached indices corresponding to the cumulative reaction progress fields.
         */
        std::vector<unsigned int> reaction_field_indices;

        /**
         * Reference viscosity. Units: Pa s
         */
        double viscosity;

        /**
         * Minimum and maximum viscosity cutoffs. Units: Pa s
         */
        double minimum_viscosity;
        double maximum_viscosity;

        /**
         * Exponent A in the temperature dependence term exp(-A * (T - T_ref) / T_ref). Dimensionless
         */
        double thermal_viscosity_exponent;

        /**
         * Per-phase viscosity prefactors. Dimensionless
         */
        std::vector<double> phase_viscosity_prefactors;

        /**
         * Width of the tanh viscosity transition in Gibbs energy space. Smaller values give sharper transitions. Units: J/mol
         */
        double gibbs_viscosity_width;

        /**
         * Reference thermal conductivity. Units: W/m/K
         */
        double k;

        /**
         * Use dynamic or hydrostatic pressure to drive reaction kinetics.
         */
        bool use_dynamic_pressure_correction_for_density;
        bool use_dynamic_pressure_correction_for_gibbs;
    };



    /**
     * Additional output fields for the MTZKinetics material model.
     */
    template <int dim>
    struct MTZKineticsOutputs : public NamedAdditionalMaterialOutputs<dim>
    {
      MTZKineticsOutputs(const unsigned int n_points, const unsigned int n_phases);

      /**
       * Return the output vector corresponding to the requested output index.
       */
      std::vector<double>
      get_nth_output(const unsigned int idx) const override;

      // General diagnostics
      std::vector<double> visc_temperature_dependence;

      // Phase mass fraction array (size: [N_phases][n_q_points])
      std::vector<std::vector<double>> phase_mass_fractions;

      // Kinetics factor arrays (size: [N_reactions][n_q_points])
      std::vector<std::vector<double>> arrhenius_factors;
      std::vector<std::vector<double>> thermodynamic_factors;
    };
  } // namespace MaterialModel
} // namespace aspect

#endif
