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

#include "phase-transition-kinetics.h"

#include <deal.II/base/numbers.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>

#include <aspect/global.h>

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    PhaseTransitionKinetics<dim>::PhaseTransitionKinetics()
      : rho_a_idx(numbers::invalid_unsigned_int)
      , rho_b_idx(numbers::invalid_unsigned_int)
      , alpha_a_idx(numbers::invalid_unsigned_int)
      , alpha_b_idx(numbers::invalid_unsigned_int)
      , beta_a_idx(numbers::invalid_unsigned_int)
      , beta_b_idx(numbers::invalid_unsigned_int)
      , cp_a_idx(numbers::invalid_unsigned_int)
      , cp_b_idx(numbers::invalid_unsigned_int)
      , dG_idx(numbers::invalid_unsigned_int)
      , dS_idx(numbers::invalid_unsigned_int)
      , dV_idx(numbers::invalid_unsigned_int)
      , Vp_a_idx(numbers::invalid_unsigned_int)
      , Vp_b_idx(numbers::invalid_unsigned_int)
      , Vs_a_idx(numbers::invalid_unsigned_int)
      , Vs_b_idx(numbers::invalid_unsigned_int)
      , dVp_dT_a_idx(numbers::invalid_unsigned_int)
      , dVp_dT_b_idx(numbers::invalid_unsigned_int)
      , dVs_dT_a_idx(numbers::invalid_unsigned_int)
      , dVs_dT_b_idx(numbers::invalid_unsigned_int)
    {}



    template <int dim>
    void
    PhaseTransitionKinetics<dim>::initialize()
    {
      // Initialize the data reader
      profile.initialize(this->get_mpi_communicator());

      // Get column indices
      rho_a_idx = profile.get_column_index_from_name("density_a");
      rho_b_idx = profile.get_column_index_from_name("density_b");
      alpha_a_idx = profile.get_column_index_from_name("thermal_expansivity_a");
      alpha_b_idx = profile.get_column_index_from_name("thermal_expansivity_b");
      beta_a_idx = profile.get_column_index_from_name("compressibility_a");
      beta_b_idx = profile.get_column_index_from_name("compressibility_b");
      cp_a_idx = profile.get_column_index_from_name("specific_heat_a");
      cp_b_idx = profile.get_column_index_from_name("specific_heat_b");
      dG_idx = profile.get_column_index_from_name("delta_molar_gibbs");
      dS_idx = profile.get_column_index_from_name("delta_molar_entropy");
      dV_idx = profile.get_column_index_from_name("delta_molar_volume");

      // Get optional column indices
      Vp_a_idx = profile.maybe_get_column_index_from_name("pressure_wave_velocity_a");
      Vp_b_idx = profile.maybe_get_column_index_from_name("pressure_wave_velocity_b");
      dVp_dT_a_idx = profile.maybe_get_column_index_from_name("pressure_wave_velocity_T_derivative_a");
      dVp_dT_b_idx = profile.maybe_get_column_index_from_name("pressure_wave_velocity_T_derivative_b");
      Vs_a_idx = profile.maybe_get_column_index_from_name("shear_wave_velocity_a");
      Vs_b_idx = profile.maybe_get_column_index_from_name("shear_wave_velocity_b");
      dVs_dT_a_idx = profile.maybe_get_column_index_from_name("shear_wave_velocity_T_derivative_a");
      dVs_dT_b_idx = profile.maybe_get_column_index_from_name("shear_wave_velocity_T_derivative_b");
    }



    template <int dim>
    bool
    PhaseTransitionKinetics<dim>::is_compressible() const
    {
      return true;
    }



    template <int dim>
    void
    PhaseTransitionKinetics<dim>::evaluate(const MaterialModel::MaterialModelInputs<dim> &in, MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      // Get density field index
      const unsigned int projected_density_idx = this->introspection().compositional_index_for_name("density_field");

      // Get X field index
      const unsigned int X_idx = this->introspection().compositional_index_for_name("X_field");

      // Set up prescribed field outputs for PDA
      std::shared_ptr<PrescribedFieldOutputs<dim>> prescribed_field_out = out.template get_additional_output_object<PrescribedFieldOutputs<dim>>();

      // Set up reaction rate outputs
      std::shared_ptr<ReactionRateOutputs<dim>> reaction_rate_out = out.template get_additional_output_object<ReactionRateOutputs<dim>>();

      // Set up phase transition kinetics outputs
      std::shared_ptr<PhaseTransitionKineticsOutputs<dim>> phase_transition_kinetics_out = out.template get_additional_output_object<PhaseTransitionKineticsOutputs<dim>>();

      // Set up seismic velocity outputs
      std::shared_ptr<SeismicAdditionalOutputs<dim>> seismic_out = out.template get_additional_output_object<SeismicAdditionalOutputs<dim>>();

      for (unsigned int q = 0; q < in.n_evaluation_points(); ++q)
        {
          // Get mesh coordinate
          const Point<dim> position = in.position[q];

          // Get local PTX conditions at spatial mesh coordinate
          const double P = in.pressure[q];
          const double T = in.temperature[q];
          const double X = in.composition[q][X_idx];

          // Get material data from adiabatic profile at adiabatic pressure coordinate
          const double P_adiabatic = this->get_adiabatic_conditions().pressure(position);
          const double T_adiabatic = this->get_adiabatic_conditions().temperature(position);
          const Point<1> profile_pos(P_adiabatic);

          // Initialize material property arrays
          std::vector<double> densities(this->introspection().n_compositional_fields);
          std::vector<double> alphas(this->introspection().n_compositional_fields);
          std::vector<double> betas(this->introspection().n_compositional_fields);
          std::vector<double> cps(this->introspection().n_compositional_fields);
          std::vector<double> Vps(this->introspection().n_compositional_fields);
          std::vector<double> Vss(this->introspection().n_compositional_fields);
          std::vector<double> dVp_dTs(this->introspection().n_compositional_fields);
          std::vector<double> dVs_dTs(this->introspection().n_compositional_fields);

          // Fill material property arrays with material data from profile
          densities[0] = profile.get_data_component(profile_pos, rho_a_idx);
          densities[1] = profile.get_data_component(profile_pos, rho_b_idx);
          alphas[0] = profile.get_data_component(profile_pos, alpha_a_idx);
          alphas[1] = profile.get_data_component(profile_pos, alpha_b_idx);
          betas[0] = profile.get_data_component(profile_pos, beta_a_idx);
          betas[1] = profile.get_data_component(profile_pos, beta_b_idx);
          cps[0] = profile.get_data_component(profile_pos, cp_a_idx);
          cps[1] = profile.get_data_component(profile_pos, cp_b_idx);

          // Fill seismic velocity arrays with material data from profile (if they exist)
          if (seismic_out != nullptr && in.requests_property(MaterialProperties::additional_outputs))
            {
              if (Vp_a_idx != numbers::invalid_unsigned_int && Vp_b_idx != numbers::invalid_unsigned_int)
                {
                  Vps[0] = profile.get_data_component(profile_pos, Vp_a_idx);
                  Vps[1] = profile.get_data_component(profile_pos, Vp_b_idx);
                }
              if (Vs_a_idx != numbers::invalid_unsigned_int && Vs_b_idx != numbers::invalid_unsigned_int)
                {
                  Vss[0] = profile.get_data_component(profile_pos, Vs_a_idx);
                  Vss[1] = profile.get_data_component(profile_pos, Vs_b_idx);
                }
              if (dVp_dT_a_idx != numbers::invalid_unsigned_int && dVp_dT_b_idx != numbers::invalid_unsigned_int)
                {
                  Vps[0] += profile.get_data_component(profile_pos, dVp_dT_a_idx) * (T - T_adiabatic);
                  Vps[1] += profile.get_data_component(profile_pos, dVp_dT_b_idx) * (T - T_adiabatic);
                }
              if (dVs_dT_a_idx != numbers::invalid_unsigned_int && dVs_dT_b_idx != numbers::invalid_unsigned_int)
                {
                  Vss[0] += profile.get_data_component(profile_pos, dVs_dT_a_idx) * (T - T_adiabatic);
                  Vss[1] += profile.get_data_component(profile_pos, dVs_dT_b_idx) * (T - T_adiabatic);
                }
            }

          // Get thermodynamic terms from profile
          const double dG = profile.get_data_component(profile_pos, dG_idx);
          const double dS = profile.get_data_component(profile_pos, dS_idx);
          const double dV = profile.get_data_component(profile_pos, dV_idx);

          if (in.requests_property(MaterialProperties::viscosity))
            {
              // Calculate viscosity temperature factor
              double visc_temperature_dependence = std::max(std::min(std::exp(-thermal_viscosity_exponent * (T - T_adiabatic) / T_adiabatic), 1e3), 1e-3);
              if (std::isnan(visc_temperature_dependence))
                visc_temperature_dependence = 1.0;

              // Calculate viscosity depth factor
              double visc_depth_dependence = viscosity_prefactors[0];
              for (unsigned int z = 0; z < transition_depths.size(); ++z)
                {
                  const Point<dim> transition_point = this->get_geometry_model().representative_point(transition_depths[z]);
                  const double transition_pressure = this->get_adiabatic_conditions().pressure(transition_point);
                  if (P_adiabatic > transition_pressure)
                    visc_depth_dependence = viscosity_prefactors[z + 1];
                }

              // Calculate viscosity
              const double eta = viscosity * visc_temperature_dependence * visc_depth_dependence;
              const double eta_effective = std::min(std::max(eta, minimum_viscosity), maximum_viscosity);

              // Update viscosity
              out.viscosities[q] = eta_effective;
            }

          // Clamp composition to valid range
          const double X_clamped = std::max(0.0, std::min(1.0, X));

          // Fill mass fractions array using clamped mass fraction
          std::vector<double> mass_fractions(this->introspection().n_compositional_fields, 1.0);
          mass_fractions[0] = 1 - X_clamped;
          mass_fractions[1] = X_clamped;

          // Compute volume fractions from mass fraction and material densities
          const std::vector<double> volume_fractions = MaterialUtilities::compute_volumes_from_masses(mass_fractions, densities, true);

          // Average material properties using computed mass fraction X
          const double rho = MaterialUtilities::average_value(volume_fractions, densities, MaterialUtilities::arithmetic);
          const double beta = MaterialUtilities::average_value(volume_fractions, betas, MaterialUtilities::arithmetic);
          const double alpha = MaterialUtilities::average_value(volume_fractions, alphas, MaterialUtilities::arithmetic);
          const double cp = MaterialUtilities::average_value(mass_fractions, cps, MaterialUtilities::arithmetic);

          // Average seismic wave velocities (if they exist)
          if (seismic_out != nullptr && in.requests_property(MaterialProperties::additional_outputs))
            {
              if (Vp_a_idx != numbers::invalid_unsigned_int && Vp_b_idx != numbers::invalid_unsigned_int)
                {
                  const double Vp = MaterialUtilities::average_value(volume_fractions, Vps, MaterialUtilities::arithmetic);
                  seismic_out->vp[q] = Vp;
                }
              if (Vs_a_idx != numbers::invalid_unsigned_int && Vs_b_idx != numbers::invalid_unsigned_int)
                {
                  const double Vs = MaterialUtilities::average_value(volume_fractions, Vss, MaterialUtilities::arithmetic);
                  seismic_out->vs[q] = Vs;
                }
              if (dVp_dT_a_idx != numbers::invalid_unsigned_int && dVp_dT_b_idx != numbers::invalid_unsigned_int)
                {
                  const double dVp_dT = MaterialUtilities::average_value(volume_fractions, dVp_dTs, MaterialUtilities::arithmetic);
                  seismic_out->vp[q] += dVp_dT * (T - T_adiabatic);
                }
              if (dVs_dT_a_idx != numbers::invalid_unsigned_int && dVs_dT_b_idx != numbers::invalid_unsigned_int)
                {
                  const double dVs_dT = MaterialUtilities::average_value(volume_fractions, dVs_dTs, MaterialUtilities::arithmetic);
                  seismic_out->vs[q] += dVs_dT * (T - T_adiabatic);
                }
            }

          // Calculate density
          const double density_factor = (1.0 - alpha * (T - T_adiabatic)) * (1.0 + beta * (P - P_adiabatic));
          const double final_rho = rho * density_factor;

          // Update material model
          out.densities[q] = final_rho;
          out.thermal_expansion_coefficients[q] = alpha;
          out.compressibilities[q] = beta;
          out.thermal_conductivities[q] = k;
          out.specific_heat[q] = cp;
          out.entropy_derivative_pressure[q] = 0.0;
          out.entropy_derivative_temperature[q] = 0.0;

          // Set all reaction terms to zero
          for (unsigned int c = 0; c < this->introspection().n_compositional_fields; ++c)
            out.reaction_terms[q][c] = 0.0;

          // Compute reaction rates at each evaluation point
          if (reaction_rate_out != nullptr)
            {
              // Calculate driving force:
              // ΔG = ΔG₀ + (P - P_adiabatic)ΔV - (T - T_adiabatic)ΔS
              const double driving_force = dG + (P - P_adiabatic) * dV - (T - T_adiabatic) * dS;

              // Get time scale
              const double time_scale = this->convert_output_to_years() ? year_in_seconds : 1.0;

              // Calculate reaction rate: dX/dt = Q * ΔG * (1 - X)
              const double rate = (driving_force < 0) ? Q_kinetic_prefactor * std::abs(driving_force) * (1.0 - X_clamped) / time_scale : 0.0;

              // Update additional named outputs
              phase_transition_kinetics_out->driving_force[q] = driving_force;

              // Update reaction rates
              reaction_rate_out->reaction_rates[q][X_idx] = rate;

              // Set other compositional fields to zero reaction rate
              for (unsigned int c = 0; c < this->introspection().n_compositional_fields; ++c)
                {
                  if (c != X_idx)
                    reaction_rate_out->reaction_rates[q][c] = 0.0;
                }
            }
        }

      // Calculate projected density reaction terms
      if (projected_density_idx != numbers::invalid_unsigned_int)
        {
          for (unsigned int q = 0; q < in.n_evaluation_points(); ++q)
            {
              out.reaction_terms[q][projected_density_idx] = out.densities[q] - in.composition[q][projected_density_idx];
            }
        }

      // Update projected density field
      if (prescribed_field_out != nullptr && projected_density_idx != numbers::invalid_unsigned_int)
        {
          for (unsigned int i = 0; i < in.position.size(); ++i)
            prescribed_field_out->prescribed_field_outputs[i][projected_density_idx] = out.densities[i];
        }
    }



    template <int dim>
    void
    PhaseTransitionKinetics<dim>::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Phase transition kinetics");
        {
          Utilities::AsciiDataProfile<dim>::declare_parameters(prm,
                                                               "$ASPECT_SOURCE_DIR/cookbooks/phase_transition_kinetics/",
                                                               "thermodynamic-driving-force-profile-olivine-wadsleyite.txt");
          prm.declare_entry("Data directory",
                            "$ASPECT_SOURCE_DIR/cookbooks/phase_transition_kinetics/",
                            Patterns::DirectoryName(),
                            "The name of a directory where the model data is found. This path may either be "
                            "absolute (if starting with '/') or relative to the current directory.");
          prm.declare_entry("Data file name",
                            "thermodynamic-driving-force-profile-olivine-wadsleyite.txt",
                            Patterns::Anything(),
                            "File that stores the thermodynamic data used for computing the equation of state and phase reaction. "
                            "The thermodynamic data are evaluated along a reference adiabatic profile. Data is read from a "
                            "tab-separated .txt file with the following required columns: "
                            "pressure density_a density_b thermal_expansivity_a thermal_expansivity_b "
                            "specific_heat_a specific_heat_b compressibility_a compressibility_b "
                            "delta_molar_gibbs delta_molar_entropy delta_molar_volume"
                            "Note: 'a' and 'b' represent phases 'a' and 'b'");
          prm.declare_entry("Viscosity", "1e21", Patterns::Double(0.), "Viscosity. Units: Pa s");
          prm.declare_entry("Minimum viscosity", "1e19", Patterns::Double(0), "The minimum viscosity cutoff. Units: Pa s");
          prm.declare_entry("Maximum viscosity", "1e24", Patterns::Double(0), "The maximum viscosity cutoff. Units: Pa s");
          prm.declare_entry("Thermal viscosity exponent",
                            "0.0",
                            Patterns::Double(0.),
                            "The temperature dependence of viscosity. Dimensionless exponent.");
          prm.declare_entry("Transition depths",
                            "150e3, 410e3, 660e3",
                            Patterns::List(Patterns::Double(0.)),
                            "A list of depths where the viscosity changes. Values must monotonically increase. Units: m");
          prm.declare_entry("Viscosity prefactors",
                            "10.0, 0.1, 1.0, 10.0",
                            Patterns::List(Patterns::Double(0.)),
                            "A list of prefactors for the viscosity that determine the viscosity profile. Each prefactor "
                            "is applied in a depth range specified by the list of `Transition depths', i.e. the first "
                            "prefactor is applied above the first transition depth, the second one between the first and "
                            "second transition depth, and so on. To compute the viscosity profile, this prefactor is multiplied "
                            "by the reference viscosity specified through the parameter `Viscosity'. List must have one more "
                            "entry than Transition depths. Units: non-dimensional.");
          prm.declare_entry("Thermal conductivity", "4.0", Patterns::Double(0.), "Reference thermal conductivity. Units: W/m/K");
          prm.declare_entry("Kinetic prefactor Q",
                            "1e-5",
                            Patterns::Double(0.0),
                            "The scalar kinetic prefactor Q in the reaction rate equation dX/dt = Q * driving_force * (1-X). "
                            "Units: mol/J/s");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    PhaseTransitionKinetics<dim>::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Phase transition kinetics");
        {
          profile.parse_parameters(prm);
          viscosity = prm.get_double("Viscosity");
          minimum_viscosity = prm.get_double("Minimum viscosity");
          maximum_viscosity = prm.get_double("Maximum viscosity");
          thermal_viscosity_exponent = prm.get_double("Thermal viscosity exponent");
          transition_depths = Utilities::string_to_double(Utilities::split_string_list(prm.get("Transition depths")));
          viscosity_prefactors = Utilities::string_to_double(Utilities::split_string_list(prm.get("Viscosity prefactors")));
          k = prm.get_double("Thermal conductivity");
          Q_kinetic_prefactor = prm.get_double("Kinetic prefactor Q");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // Validate parameters
      AssertThrow(this->introspection().n_compositional_fields >= 1,
                  ExcMessage("Phase transition kinetics model requires at least one compositional field."));
      AssertThrow(Q_kinetic_prefactor > 0.0,
                  ExcMessage("Kinetic prefactor Q must be positive."));
      if (viscosity_prefactors.size() != transition_depths.size() + 1)
        AssertThrow(false,
                    ExcMessage("Error: The list of Viscosity prefactors needs to have exactly one more "
                               "entry than the list of Transition depths."));

      // Declare dependencies on solution variables
      this->model_dependence.viscosity = NonlinearDependence::none;
      this->model_dependence.density = NonlinearDependence::pressure | NonlinearDependence::temperature;
      this->model_dependence.compressibility = NonlinearDependence::none;
      this->model_dependence.specific_heat = NonlinearDependence::none;
      this->model_dependence.thermal_conductivity = NonlinearDependence::none;
    }



    template <int dim>
    void
    PhaseTransitionKinetics<dim>::create_additional_named_outputs(
      MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      // Create reaction rate outputs
      if (out.template get_additional_output_object<ReactionRateOutputs<dim>>() == nullptr)
        {
          const unsigned int n_points = out.n_evaluation_points();
          out.additional_outputs.push_back(std::make_unique<MaterialModel::ReactionRateOutputs<dim>>(n_points, this->introspection().n_compositional_fields));
        }

      // Create prescribed field outputs for PDA
      if (this->introspection().composition_type_exists(CompositionalFieldDescription::density) &&
          out.template get_additional_output_object<PrescribedFieldOutputs<dim>>() == nullptr)
        {
          const unsigned int n_points = out.n_evaluation_points();
          out.additional_outputs.push_back(std::make_unique<MaterialModel::PrescribedFieldOutputs<dim>>(n_points, this->introspection().n_compositional_fields));
        }

      // Create phase transition kinetics outputs
      if (out.template get_additional_output_object<PhaseTransitionKineticsOutputs<dim>>() == nullptr)
        {
          const unsigned int n_points = out.n_evaluation_points();
          out.additional_outputs.push_back(std::make_unique<MaterialModel::PhaseTransitionKineticsOutputs<dim>>(n_points));
        }

      // Create seismic velocity outputs
      if (out.template get_additional_output_object<SeismicAdditionalOutputs<dim>>() == nullptr)
        {
          const unsigned int n_points = out.n_evaluation_points();
          out.additional_outputs.push_back(std::make_unique<MaterialModel::SeismicAdditionalOutputs<dim>>(n_points));
        }
    }



    namespace
    {
      std::vector<std::string> make_additional_output_names()
      {
        std::vector<std::string> names;
        names.emplace_back("driving_force");
        return names;
      }
    }

    template <int dim>
    PhaseTransitionKineticsOutputs<dim>::PhaseTransitionKineticsOutputs(const unsigned int n_points)
      : NamedAdditionalMaterialOutputs<dim>(make_additional_output_names())
      , driving_force(n_points, numbers::signaling_nan<double>())
    {}



    template <int dim>
    std::vector<double>
    PhaseTransitionKineticsOutputs<dim>::get_nth_output(const unsigned int idx) const
    {
      AssertIndexRange (idx, 1);
      switch (idx)
        {
          case 0:
            return driving_force;

          default:
            AssertThrow(false, ExcInternalError());
        }
      return driving_force;
    }
  } // namespace MaterialModel
} // namespace aspect

namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(
      PhaseTransitionKinetics,
      "phase transition kinetics",
      "Models a phase transition using reaction kinetics that follow the equation: dX/dt = Q * ΔG * (1 - X), where Q is a "
      "user-defined kinetic prefactor and ΔG is the thermodynamic driving force calculated from internally-consistent "
      "thermodynamic data. The thermodynamic driving force is computed as ΔG = ΔG₀ + (P - P_adiabatic)ΔV - (T - T_adiabatic)ΔS, "
      "where all thermodynamic parameters are read from a user-specified ASCII data file. Requires at least one compositional "
      "field representing the reacting phase.")
  } // namespace MaterialModel
} // namespace aspect
