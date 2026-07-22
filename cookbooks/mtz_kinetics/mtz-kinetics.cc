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

#include "mtz-kinetics.h"

#include <deal.II/base/numbers.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/patterns.h>
#include <deal.II/base/point.h>

#include <aspect/global.h>

namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    MTZKinetics<dim>::MTZKinetics() = default;



    template <int dim>
    void
    MTZKinetics<dim>::initialize()
    {
      reaction_chain.initialize_simulator(this->get_simulator());
      reaction_chain.initialize();

      const unsigned int n_reactions = reaction_chain.n_reactions();
      const unsigned int n_phases    = n_reactions + 1;

      const unsigned int projected_density_idx = this->introspection().compositional_name_exists("density_field")
                                                 ? this->introspection().compositional_index_for_name("density_field")
                                                 : numbers::invalid_unsigned_int;

      reaction_field_indices.clear();
      for (unsigned int c = 0; c < this->introspection().n_compositional_fields; ++c)
        if (c != projected_density_idx)
          reaction_field_indices.push_back(c);

      AssertThrow(reaction_field_indices.size() == n_reactions,
                  ExcMessage("The number of available compositional fields for reactions does not match the number of reactions in the chain."));

      // Size property index arrays
      rho_idx.assign(n_phases, numbers::invalid_unsigned_int);
      alpha_idx.assign(n_phases, numbers::invalid_unsigned_int);
      beta_idx.assign(n_phases, numbers::invalid_unsigned_int);
      cp_idx.assign(n_phases, numbers::invalid_unsigned_int);
      Vp_idx.assign(n_phases, numbers::invalid_unsigned_int);
      Vs_idx.assign(n_phases, numbers::invalid_unsigned_int);
      dVp_dT_idx.assign(n_phases, numbers::invalid_unsigned_int);
      dVs_dT_idx.assign(n_phases, numbers::invalid_unsigned_int);

      dG_idx.assign(n_reactions, numbers::invalid_unsigned_int);
      dS_idx.assign(n_reactions, numbers::invalid_unsigned_int);
      dV_idx.assign(n_reactions, numbers::invalid_unsigned_int);

      // Initialize profiles and assign indices
      for (unsigned int r = 0; r < n_reactions; ++r)
        {
          profiles[r].initialize(this->get_mpi_communicator());

          dG_idx[r] = profiles[r].get_column_index_from_name("delta_molar_gibbs");
          dS_idx[r] = profiles[r].get_column_index_from_name("delta_molar_entropy");
          dV_idx[r] = profiles[r].get_column_index_from_name("delta_molar_volume");
        }

      for (unsigned int p = 0; p < n_phases; ++p)
        {
          const bool from_a        = (p == 0);
          const unsigned int r     = from_a ? 0 : p - 1;
          const std::string suffix = from_a ? "_a" : "_b";

          rho_idx[p]    = profiles[r].get_column_index_from_name("density" + suffix);
          alpha_idx[p]  = profiles[r].get_column_index_from_name("thermal_expansivity" + suffix);
          beta_idx[p]   = profiles[r].get_column_index_from_name("compressibility" + suffix);
          cp_idx[p]     = profiles[r].get_column_index_from_name("specific_heat" + suffix);
          Vp_idx[p]     = profiles[r].maybe_get_column_index_from_name("pressure_wave_velocity" + suffix);
          dVp_dT_idx[p] = profiles[r].maybe_get_column_index_from_name("pressure_wave_velocity_T_derivative" + suffix);
          Vs_idx[p]     = profiles[r].maybe_get_column_index_from_name("shear_wave_velocity" + suffix);
          dVs_dT_idx[p] = profiles[r].maybe_get_column_index_from_name("shear_wave_velocity_T_derivative" + suffix);
        }
    }



    template <int dim>
    bool
    MTZKinetics<dim>::is_compressible() const
    {
      return true;
    }



    template <int dim>
    void
    MTZKinetics<dim>::evaluate(const MaterialModel::MaterialModelInputs<dim> &in, MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      // Set up additional output objects
      std::shared_ptr<PrescribedFieldOutputs<dim>>   prescribed_field_out = out.template get_additional_output_object<PrescribedFieldOutputs<dim>>();
      std::shared_ptr<ReactionRateOutputs<dim>>      reaction_rate_out    = out.template get_additional_output_object<ReactionRateOutputs<dim>>();
      std::shared_ptr<MTZKineticsOutputs<dim>>       mtz_kinetics_out     = out.template get_additional_output_object<MTZKineticsOutputs<dim>>();
      std::shared_ptr<SeismicAdditionalOutputs<dim>> seismic_out          = out.template get_additional_output_object<SeismicAdditionalOutputs<dim>>();

      const bool output_diagnostics  = (mtz_kinetics_out != nullptr);
      const unsigned int n_q_points  = in.n_evaluation_points();
      const unsigned int n_reactions = reaction_chain.n_reactions();
      const unsigned int n_phases    = n_reactions + 1;

      // Allocate space once outside the quadrature loop
      if (output_diagnostics)
        {
          mtz_kinetics_out->visc_temperature_dependence.resize(n_q_points);
          mtz_kinetics_out->phase_mass_fractions.assign(n_phases, std::vector<double>(n_q_points, 0.0));
        }

      // Get compositional field indices
      const unsigned int projected_density_idx = this->introspection().compositional_name_exists("density_field")
                                                 ? this->introspection().compositional_index_for_name("density_field")
                                                 : numbers::invalid_unsigned_int;

      for (unsigned int q = 0; q < in.n_evaluation_points(); ++q)
        {
          const Point<dim> position = in.position[q];

          const double P = in.pressure[q];
          const double T = in.temperature[q];

          const double P_ref = this->get_adiabatic_conditions().pressure(position);
          const double T_ref = this->get_adiabatic_conditions().temperature(position);
          const Point<1> profile_pos(P_ref);

          // Initialize reaction rates to zero for all compositional fields
          if (reaction_rate_out != nullptr)
            for (unsigned int c = 0; c < this->introspection().n_compositional_fields; ++c)
              reaction_rate_out->reaction_rates[q][c] = 0.0;

          // Get cumulative reaction progress from compositional fields
          std::vector<double> reaction_progress_values(n_reactions);
          for (unsigned int r = 0; r < n_reactions; ++r)
            reaction_progress_values[r] = in.composition[q][reaction_field_indices[r]];

          // Ensure physical range [0, 1] and monotonicity xi[i] <= xi[i-1]
          const std::vector<double> reaction_progress_clamped = reaction_chain.clamp_cumulative_progress(reaction_progress_values);

          // Storage for all phases
          // 0 = olivine (ol), 1 = wadsleyite (wd), 2 = ringwoodite (ri), 3 = post-spinel (ps)
          EquationOfStateOutputs<dim> eos_outputs(n_phases);

          // Seismic velocities are not part of EquationOfStateOutputs
          std::vector<double> Vps(n_phases, 0.0), Vss(n_phases, 0.0);

          // Get properties from profiles
          for (unsigned int p = 0; p < n_phases; ++p)
            {
              const bool from_a    = (p == 0);
              const unsigned int r = from_a ? 0 : p - 1;

              eos_outputs.densities[p]                      = profiles[r].get_data_component(profile_pos, rho_idx[p]);
              eos_outputs.thermal_expansion_coefficients[p] = profiles[r].get_data_component(profile_pos, alpha_idx[p]);
              eos_outputs.compressibilities[p]              = profiles[r].get_data_component(profile_pos, beta_idx[p]);
              eos_outputs.specific_heat_capacities[p]       = profiles[r].get_data_component(profile_pos, cp_idx[p]);

              if (seismic_out != nullptr && Vp_idx[p] != numbers::invalid_unsigned_int)
                {
                  Vps[p] = profiles[r].get_data_component(profile_pos, Vp_idx[p]);
                  if (dVp_dT_idx[p] != numbers::invalid_unsigned_int)
                    Vps[p] += profiles[r].get_data_component(profile_pos, dVp_dT_idx[p]) * (T - T_ref);
                }
              if (seismic_out != nullptr && Vs_idx[p] != numbers::invalid_unsigned_int)
                {
                  Vss[p] = profiles[r].get_data_component(profile_pos, Vs_idx[p]);
                  if (dVs_dT_idx[p] != numbers::invalid_unsigned_int)
                    Vss[p] += profiles[r].get_data_component(profile_pos, dVs_dT_idx[p]) * (T - T_ref);
                }
            }

          // Latent heat terms are handled separately via ReactionRateOutputs and the operator-splitting scheme
          std::fill(eos_outputs.entropy_derivative_pressure.begin(),    eos_outputs.entropy_derivative_pressure.end(),    0.0);
          std::fill(eos_outputs.entropy_derivative_temperature.begin(), eos_outputs.entropy_derivative_temperature.end(), 0.0);

          // Compute actual phase mass fractions from cumulative variables
          const std::vector<double> mass_fractions = reaction_chain.compute_phase_mass_fractions(reaction_progress_clamped, eos_outputs.densities);

          // Populate phase mass fractions in diagnostics
          if (output_diagnostics)
            for (unsigned int p = 0; p < n_phases; ++p)
              mtz_kinetics_out->phase_mass_fractions[p][q] = mass_fractions[p];

          // Compute volume fractions from mass fractions
          const std::vector<double> volume_fractions = MaterialUtilities::compute_volumes_from_masses(mass_fractions, eos_outputs.densities, true);

          // Update material model outputs (weighted averaging)
          MaterialUtilities::fill_averaged_equation_of_state_outputs(eos_outputs, mass_fractions, volume_fractions, q, out);

          // Seismic velocities are handled separately since EquationOfStateOutputs has no Vp/Vs fields
          if (seismic_out != nullptr && in.requests_property(MaterialProperties::additional_outputs))
            {
              if (Vps[0] > 0.0)
                seismic_out->vp[q] = MaterialUtilities::average_value(volume_fractions, Vps, MaterialUtilities::arithmetic);
              if (Vss[0] > 0.0)
                seismic_out->vs[q] = MaterialUtilities::average_value(volume_fractions, Vss, MaterialUtilities::arithmetic);
            }

          // Thermal conductivity held constant
          out.thermal_conductivities[q] = k;

          // Apply temperature and pressure corrections to density
          const double temperature_correction_density = (T - T_ref) * out.thermal_expansion_coefficients[q];
          const double pressure_correction_density    = use_dynamic_pressure_correction_for_density ? (P - P_ref) * out.compressibilities[q] : 0.0;
          const double density_factor                 = (1.0 - temperature_correction_density) * (1.0 + pressure_correction_density);
          out.densities[q]                             = out.densities[q] * density_factor;

          // Compute viscosity temperature dependence: exp(-A * (T - T_ref) / T_ref)
          double visc_temperature_dependence = std::max(std::min(std::exp(-thermal_viscosity_exponent * (T - T_ref) / T_ref), 1e3), 1e-3);
          if (std::isnan(visc_temperature_dependence))
            visc_temperature_dependence = 1.0;

          if (output_diagnostics)
            mtz_kinetics_out->visc_temperature_dependence[q] = visc_temperature_dependence;

          // Compute viscosity
          if (in.requests_property(MaterialProperties::viscosity))
            {
              // ---------------------------------------------------------------
              // Gibbs-based phase viscosity prefactor
              //
              // Viscosity jumps are anchored to the equilibrium phase boundaries (where dG = 0 along the reference adiabat) using tanh transitions
              // of width gibbs_viscosity_width (J/mol). This is fully decoupled from the kinetic reaction_progress_values fields, so metastable
              // nuclei in the slab do not produce spurious viscosity patches.
              //
              // Transition fractions f_r -> 1 when the high-pressure phase is stable (dG < 0).
              // The chained construction dynamically builds stability fractions (partition of unity) for N phases.
              // ---------------------------------------------------------------
              const double w = gibbs_viscosity_width;

              std::vector<double> f_transition(n_reactions, 0.0);
              for (unsigned int r = 0; r < n_reactions; ++r)
                {
                  const double dG_val = profiles[r].get_data_component(profile_pos, dG_idx[r]);
                  f_transition[r] = 0.5 * (1.0 - std::tanh(dG_val / w));
                }

              std::vector<double> phi(n_phases, 1.0);
              double accumulated_prod = 1.0;
              for (unsigned int p = 0; p < n_phases; ++p)
                {
                  if (p == 0)
                    phi[p] = 1.0 - f_transition[0];
                  else if (p < n_reactions)
                    {
                      accumulated_prod *= f_transition[p - 1];
                      phi[p] = accumulated_prod * (1.0 - f_transition[p]);
                    }
                  else
                    {
                      accumulated_prod *= f_transition[p - 1];
                      phi[p] = accumulated_prod;
                    }
                }

              // Log-linear (geometric) average of per-phase prefactors
              double log_prefactor = 0.0;
              for (unsigned int p = 0; p < n_phases; ++p)
                log_prefactor += phi[p] * std::log(phase_viscosity_prefactors[p]);

              const double visc_depth_dependence = std::exp(log_prefactor);
              const double eta                   = viscosity * visc_temperature_dependence * visc_depth_dependence;
              const double eta_effective         = std::min(std::max(eta, minimum_viscosity), maximum_viscosity);

              out.viscosities[q] = eta_effective;
            }

          // Evaluate reaction rates and individual kinetic factors
          if (in.requests_property(MaterialProperties::reaction_rates))
            {
              const double time_scale = this->convert_output_to_years() ? year_in_seconds : 1.0;

              for (unsigned int r = 0; r < n_reactions; ++r)
                {
                  const double dG_ref = profiles[r].get_data_component(profile_pos, dG_idx[r]);
                  const double dS     = profiles[r].get_data_component(profile_pos, dS_idx[r]);
                  const double dV     = profiles[r].get_data_component(profile_pos, dV_idx[r]);
                  const double dG     = correct_gibbs(dG_ref, dS, dV, P, P_ref, T, T_ref);

                  if (reaction_rate_out != nullptr)
                    reaction_rate_out->reaction_rates[q][reaction_field_indices[r]] =
                      reaction_chain.net_forward_reaction_rate(T, P, dG, reaction_progress_clamped[r], r) / time_scale;
                }
            }

          for (unsigned int c = 0; c < this->introspection().n_compositional_fields; ++c)
            out.reaction_terms[q][c] = 0.0;
        }

      // Calculate projected density reaction terms
      if (projected_density_idx != numbers::invalid_unsigned_int)
        for (unsigned int q = 0; q < in.n_evaluation_points(); ++q)
          out.reaction_terms[q][projected_density_idx] = out.densities[q] - in.composition[q][projected_density_idx];

      // Update projected density field
      if (prescribed_field_out != nullptr && projected_density_idx != numbers::invalid_unsigned_int)
        for (unsigned int i = 0; i < in.position.size(); ++i)
          prescribed_field_out->prescribed_field_outputs[i][projected_density_idx] = out.densities[i];
    }



    template <int dim>
    double
    MTZKinetics<dim>::correct_gibbs(const double dG_ref, const double dS, const double dV, const double P, const double P_ref, const double T, const double T_ref) const
    {
      const double dT = T - T_ref;
      const double dP = use_dynamic_pressure_correction_for_gibbs ? P - P_ref : 0.0;

      return dG_ref + (dP * dV) - (dT * dS);
    }



    template <int dim>
    void
    MTZKinetics<dim>::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("MTZ kinetics");
        {
          ReactionModel::ReactionChain<dim>::declare_parameters(prm);

          prm.declare_entry("Data directory", "data/", Patterns::DirectoryName(), "Directory containing thermodynamic data profiles.");
          prm.declare_entry("Data file names",
                            "profile.tsv",
                            Patterns::List(Patterns::Anything()),
                            "'|'-separated list of thermodynamic data files, one per reaction.");

          prm.enter_subsection("Ascii data model");
          {
            Utilities::AsciiDataProfile<dim>::declare_parameters(prm, "$ASPECT_SOURCE_DIR/cookbooks/mtz_kinetics/", "profile.tsv");
          }
          prm.leave_subsection();

          prm.declare_entry("Viscosity", "1e21", Patterns::Double(0.), "Reference viscosity. Units: Pa s");
          prm.declare_entry("Minimum viscosity", "1e19", Patterns::Double(0.), "Minimum viscosity cutoff. Units: Pa s");
          prm.declare_entry("Maximum viscosity", "1e24", Patterns::Double(0.), "Maximum viscosity cutoff. Units: Pa s");
          prm.declare_entry("Thermal viscosity exponent", "0.0", Patterns::Double(0.), "Exponent A in exp(-A*(T-T_ref)/T_ref). Dimensionless.");
          prm.declare_entry("Phase viscosity prefactors",
                            "1.0, 1.0, 1.0, 1.0",
                            Patterns::List(Patterns::Double(0.)),
                            "Comma-separated list of dimensionless viscosity prefactors for each phase. "
                            "Must be the same length as 'Phase names'.");
          prm.declare_entry("Gibbs viscosity width",
                            "1000.0",
                            Patterns::Double(0.),
                            "Width of the tanh viscosity transition in Gibbs free energy space. "
                            "Controls sharpness of viscosity jumps at equilibrium phase boundaries. "
                            "Smaller values give sharper transitions. ~1000 J/mol gives a transition "
                            "spread of a few km. Units: J/mol");

          prm.declare_entry("Thermal conductivity", "4.0", Patterns::Double(0.), "Reference thermal conductivity. Units: W/m/K");

          prm.declare_entry("Use dynamic pressure correction for density", "true", Patterns::Bool(), "Apply dynamic pressure correction to density.");
          prm.declare_entry("Use dynamic pressure correction for Gibbs",
                            "false",
                            Patterns::Bool(),
                            "Apply dynamic pressure correction to Gibbs free energy.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    MTZKinetics<dim>::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("MTZ kinetics");
        {
          data_directory  = prm.get("Data directory");
          data_file_names = Utilities::split_string_list(prm.get("Data file names"), "|");

          reaction_chain.parse_parameters(prm);

          const unsigned int n_reactions = reaction_chain.n_reactions();
          const unsigned int n_phases    = n_reactions + 1;

          AssertThrow(data_file_names.size() == n_reactions,
                      ExcMessage("The number of entries in 'Data file names' must match the number of reactions."));

          profiles.resize(n_reactions);
          for (unsigned int r = 0; r < n_reactions; ++r)
            {
              prm.enter_subsection("Ascii data model");
              profiles[r].parse_parameters(prm);
              prm.leave_subsection();

              profiles[r].data_directory = data_directory;
              profiles[r].data_file_name = data_file_names[r];
            }

          viscosity                  = prm.get_double("Viscosity");
          minimum_viscosity          = prm.get_double("Minimum viscosity");
          maximum_viscosity          = prm.get_double("Maximum viscosity");
          thermal_viscosity_exponent = prm.get_double("Thermal viscosity exponent");

          phase_viscosity_prefactors = Utilities::string_to_double(Utilities::split_string_list(prm.get("Phase viscosity prefactors")));
          AssertThrow(phase_viscosity_prefactors.size() == n_phases,
                      ExcMessage("The number of entries in 'Phase viscosity prefactors' must match the number of phases (n_reactions + 1)."));

          gibbs_viscosity_width      = prm.get_double("Gibbs viscosity width");

          k = prm.get_double("Thermal conductivity");

          use_dynamic_pressure_correction_for_density = prm.get_bool("Use dynamic pressure correction for density");
          use_dynamic_pressure_correction_for_gibbs   = prm.get_bool("Use dynamic pressure correction for Gibbs");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      AssertThrow(reaction_chain.n_reactions() > 0, ExcMessage("At least one phase transformation must be defined."));
      AssertThrow(gibbs_viscosity_width > 0.0,  ExcMessage("Gibbs viscosity width must be positive."));

      this->model_dependence.viscosity            = NonlinearDependence::none;
      this->model_dependence.density              = NonlinearDependence::pressure | NonlinearDependence::temperature;
      this->model_dependence.compressibility      = NonlinearDependence::none;
      this->model_dependence.specific_heat        = NonlinearDependence::none;
      this->model_dependence.thermal_conductivity = NonlinearDependence::none;
    }



    template <int dim>
    void
    MTZKinetics<dim>::create_additional_named_outputs(MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      if (out.template get_additional_output_object<ReactionRateOutputs<dim>>() == nullptr)
        {
          const unsigned int n_points = out.n_evaluation_points();
          out.additional_outputs.push_back(
            std::make_unique<MaterialModel::ReactionRateOutputs<dim>>(n_points, this->introspection().n_compositional_fields));
        }

      if (this->introspection().composition_type_exists(CompositionalFieldDescription::density) &&
          out.template get_additional_output_object<PrescribedFieldOutputs<dim>>() == nullptr)
        {
          const unsigned int n_points = out.n_evaluation_points();
          out.additional_outputs.push_back(
            std::make_unique<MaterialModel::PrescribedFieldOutputs<dim>>(n_points, this->introspection().n_compositional_fields));
        }

      if (out.template get_additional_output_object<MTZKineticsOutputs<dim>>() == nullptr)
        {
          const unsigned int n_points    = out.n_evaluation_points();
          const unsigned int n_reactions = reaction_chain.n_reactions();
          const unsigned int n_phases    = n_reactions + 1;
          out.additional_outputs.push_back(std::make_unique<MaterialModel::MTZKineticsOutputs<dim>>(n_points, n_phases));
        }

      if (out.template get_additional_output_object<SeismicAdditionalOutputs<dim>>() == nullptr)
        {
          const unsigned int n_points = out.n_evaluation_points();
          out.additional_outputs.push_back(std::make_unique<MaterialModel::SeismicAdditionalOutputs<dim>>(n_points));
        }
    }



    namespace
    {
      std::vector<std::string>
      make_additional_output_names(const unsigned int n_phases)
      {
        std::vector<std::string> names;
        names.reserve(1 + n_phases);

        names.emplace_back("visc_temperature_dependence");

        for (unsigned int p = 0; p < n_phases; ++p)
          names.push_back("X_phase_" + std::to_string(p));

        return names;
      }
    } // namespace

    template <int dim>
    MTZKineticsOutputs<dim>::MTZKineticsOutputs(const unsigned int n_points, const unsigned int n_phases)
      : NamedAdditionalMaterialOutputs<dim>(make_additional_output_names(n_phases))
      , visc_temperature_dependence(n_points, numbers::signaling_nan<double>())
      , phase_mass_fractions(n_phases, std::vector<double>(n_points, numbers::signaling_nan<double>()))
    {}



    template <int dim>
    std::vector<double>
    MTZKineticsOutputs<dim>::get_nth_output(const unsigned int idx) const
    {
      const unsigned int n_phases = phase_mass_fractions.size();

      AssertIndexRange(idx, 1 + n_phases);

      if (idx == 0)
        return visc_temperature_dependence;

      const unsigned int phase_idx = idx - 1;
      AssertIndexRange(phase_idx, phase_mass_fractions.size());
      return phase_mass_fractions[phase_idx];
    }
  } // namespace MaterialModel
} // namespace aspect

namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(
      MTZKinetics,
      "MTZ kinetics",
      "Models phase transformations for three reactions with stoichiometric coupling:\n"
      "(1: interface-controlled) Olivine <-> wadsleyite (Hosoya et al. 2005): dX/dt = Z * T * exp(-Ha + PVa/RT) * (1 - exp(dG/RT)) * X\n"
      "(2: interface-controlled) Wadsleyite <-> ringwoodite: dX/dt = Z * T * exp(-Ha + PVa/RT) * (1 - exp(dG/RT)) * X\n"
      "(3: diffusion-controlled) Ringwoodite <-> bridgmanite + periclase (Kubo et al. 2002): dX/dt = Z * -dG * |dG| * exp(-Ea/RT) * X\n\n"
      "Compositional fields xi_ol_wd, xi_wd_ri, xi_ri_ps are cumulative reaction progress variables.\n"
      "Actual phase fractions: X_ol = 1 - xi_ol_wd, X_wd = xi_ol_wd - xi_wd_ri, X_ri = xi_wd_ri - xi_ri_ps, X_ps = xi_ri_ps.\n\n"
      "Viscosity uses per-phase prefactors blended via tanh transitions centred on the equilibrium phase boundaries\n"
      "(where dG = 0 along the reference adiabat). This decouples the viscosity structure from the kinetic phase\n"
      "fields, preventing spurious viscosity patches from metastable nuclei in subducting slabs.\n"
      "Transition sharpness is controlled by 'Gibbs viscosity width' (J/mol).")
  } // namespace MaterialModel
} // namespace aspect
