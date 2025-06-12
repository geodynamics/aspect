#include "kinetic-driving-force.h"

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
    KineticDrivingForce<dim>::KineticDrivingForce()
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
    {}

    template <int dim>
    void
    KineticDrivingForce<dim>::initialize()
    {
      // Initialize the data reader
      profile.initialize(this->get_mpi_communicator());

      // Get column indices
      rho_a_idx   = profile.get_column_index_from_name("density_a");
      rho_b_idx   = profile.get_column_index_from_name("density_b");
      alpha_a_idx = profile.get_column_index_from_name("thermal_expansivity_a");
      alpha_b_idx = profile.get_column_index_from_name("thermal_expansivity_b");
      beta_a_idx  = profile.get_column_index_from_name("compressibility_a");
      beta_b_idx  = profile.get_column_index_from_name("compressibility_b");
      cp_a_idx    = profile.get_column_index_from_name("specific_heat_a");
      cp_b_idx    = profile.get_column_index_from_name("specific_heat_b");
      dG_idx      = profile.get_column_index_from_name("delta_molar_gibbs");
      dS_idx      = profile.get_column_index_from_name("delta_molar_entropy");
      dV_idx      = profile.get_column_index_from_name("delta_molar_volume");
    }

    template <int dim>
    bool
    KineticDrivingForce<dim>::is_compressible() const
    {
      return true;
    }

    template <int dim>
    void
    KineticDrivingForce<dim>::evaluate(
      const MaterialModel::MaterialModelInputs<dim> &in,
      MaterialModel::MaterialModelOutputs<dim>      &out) const
    {
      // Validation checks
      AssertThrow(this->introspection().n_compositional_fields >= 1,
                  ExcMessage("Kinetic driving force model requires at least "
                             "one compositional field."));

      // Get density field index
      const unsigned int projected_density_idx =
        this->introspection().compositional_index_for_name("density_field");

      // Get X field index
      const unsigned int X_idx =
        this->introspection().compositional_index_for_name("X_field");

      // Set up prescribed field outputs
      PrescribedFieldOutputs<dim> *prescribed_field_out =
        out.template get_additional_output<PrescribedFieldOutputs<dim>>();

      // Set up reaction rate outputs
      ReactionRateOutputs<dim> *reaction_rate_out =
        out.template get_additional_output<ReactionRateOutputs<dim>>();

      for (unsigned int q = 0; q < in.n_evaluation_points(); ++q)
        {
          // Get mesh coordinate
          const Point<dim> position = in.position[q];

          // Get local PT conditions at spatial coordinate
          const double P = in.pressure[q];
          const double T = in.temperature[q];
          const double X = in.composition[q][X_idx];

          // Get data from profile at adiabatic pressure coordinate
          const double P_adiabatic =
            this->get_adiabatic_conditions().pressure(position);
          const double T_adiabatic =
            this->get_adiabatic_conditions().temperature(position);
          const Point<1> profile_pos(P_adiabatic);

          double rho_a, rho_b, alpha_a, alpha_b, beta_a, beta_b, cp_a, cp_b, dG,
                 dS, dV;
          rho_a   = profile.get_data_component(profile_pos, rho_a_idx);
          rho_b   = profile.get_data_component(profile_pos, rho_b_idx);
          alpha_a = profile.get_data_component(profile_pos, alpha_a_idx);
          alpha_b = profile.get_data_component(profile_pos, alpha_b_idx);
          beta_a  = profile.get_data_component(profile_pos, beta_a_idx);
          beta_b  = profile.get_data_component(profile_pos, beta_b_idx);
          cp_a    = profile.get_data_component(profile_pos, cp_a_idx);
          cp_b    = profile.get_data_component(profile_pos, cp_b_idx);
          dG      = profile.get_data_component(profile_pos, dG_idx);
          dS      = profile.get_data_component(profile_pos, dS_idx);
          dV      = profile.get_data_component(profile_pos, dV_idx);

          // Calculate viscosity temperature factor
          double visc_temperature_dependence =
            std::max(std::min(std::exp(-thermal_viscosity_exponent *
                                       (T - T_adiabatic) / T_adiabatic),
                              1e3),
                     1e-3);
          if (std::isnan(visc_temperature_dependence))
            visc_temperature_dependence = 1.0;

          // Calculate viscosity depth factor
          double visc_depth_dependence = viscosity_prefactors[0];
          for (unsigned int z = 0; z < transition_depths.size(); ++z)
            {
              const Point<dim, double> transition_point =
                this->get_geometry_model().representative_point(
                  transition_depths[z]);
              double transition_pressure =
                this->get_adiabatic_conditions().pressure(transition_point);
              if (P_adiabatic > transition_pressure)
                visc_depth_dependence = viscosity_prefactors[z + 1];
            }

          // Calculate viscosity
          const double eta =
            viscosity * visc_temperature_dependence * visc_depth_dependence;
          const double eta_effective =
            std::min(std::max(eta, minimum_viscosity), maximum_viscosity);

          // Clamp composition to valid range
          const double X_clamped = std::max(0.0, std::min(1.0, X));

          // Average material properties
          std::vector<double> mass_fractions(
            this->introspection().n_compositional_fields, 1.0);
          mass_fractions[0] = 1 - X_clamped;
          mass_fractions[1] = X_clamped;

          std::vector<double> densities(
            this->introspection().n_compositional_fields, rho_a);
          densities[0] = rho_a;
          densities[1] = rho_b;

          std::vector<double> alphas(
            this->introspection().n_compositional_fields, alpha_a);
          alphas[0] = alpha_a;
          alphas[1] = alpha_b;

          std::vector<double> betas(
            this->introspection().n_compositional_fields, beta_a);
          betas[0] = beta_a;
          betas[1] = beta_b;

          std::vector<double> cps(this->introspection().n_compositional_fields,
                                  cp_a);
          cps[0] = cp_a;
          cps[1] = cp_b;

          const std::vector<double> volume_fractions =
            MaterialUtilities::compute_volumes_from_masses(mass_fractions,
                                                           densities,
                                                           true);

          const double rho =
            MaterialUtilities::average_value(volume_fractions,
                                             densities,
                                             MaterialUtilities::arithmetic);
          const double beta =
            MaterialUtilities::average_value(volume_fractions,
                                             betas,
                                             MaterialUtilities::arithmetic);
          const double alpha =
            MaterialUtilities::average_value(volume_fractions,
                                             alphas,
                                             MaterialUtilities::arithmetic);
          const double cp =
            MaterialUtilities::average_value(mass_fractions,
                                             cps,
                                             MaterialUtilities::arithmetic);

          // Get pressure for density calculation
          double P_for_rho =
            use_adiabatic_pressure_for_density ? P_adiabatic : P;

          // Calculate density
          const double density_factor =
            (1.0 - alpha * (T - T_adiabatic)) * (1.0 + beta * (P - P_for_rho));
          const double final_rho = rho * density_factor;

          // Update material model
          out.densities[q]                      = final_rho;
          out.thermal_expansion_coefficients[q] = alpha;
          out.compressibilities[q]              = beta;
          out.thermal_conductivities[q]         = k;
          out.specific_heat[q]                  = cp;
          out.entropy_derivative_pressure[q]    = 0.0;
          out.entropy_derivative_temperature[q] = 0.0;
          out.viscosities[q]                    = eta_effective;

          // Set all reaction terms to zero
          for (unsigned int c = 0;
               c < this->introspection().n_compositional_fields;
               ++c)
            out.reaction_terms[q][c] = 0.0;

          // Compute reaction rates at each evaluation point
          if (reaction_rate_out != nullptr)
            {
              // Apply temperature cutoff if enabled
              double temperature_factor = 1.0;
              if (use_temperature_cutoff && T < reaction_cutoff_temperature)
                {
                  const double temp_diff = reaction_cutoff_temperature - T;
                  temperature_factor =
                    std::exp(-temp_diff / 100.0); // 100K characteristic scale
                }

              // Calculate driving force:
              // ΔG = ΔG₀ + (P - P_adiabatic)ΔV - (T - T_adiabatic)ΔS
              const double driving_force =
                dG + (P - P_adiabatic) * dV - (T - T_adiabatic) * dS;

              // Get time scale
              const double time_scale =
                this->convert_output_to_years() ? year_in_seconds : 1.0;

              // Calculate reaction rate: dX/dt = Q * ΔG * (1-X)
              const double rate =
                (driving_force < 0) ?
                Q_kinetic_prefactor * std::abs(driving_force) *
                (1.0 - X_clamped) * temperature_factor / time_scale :
                0.0;

              // Limit reaction rate
              const double epsilon      = 1e-15;
              double       rate_clamped = rate;

              if (X_clamped <= 0.0 + epsilon)
                rate_clamped = 0.0;
              if (X_clamped >= 1.0 - epsilon && rate > 0)
                rate_clamped = 0.0;

              reaction_rate_out->reaction_rates[q][X_idx] = rate_clamped;

              // Set other compositional fields to zero reaction rate
              for (unsigned int c = 0;
                   c < this->introspection().n_compositional_fields;
                   ++c)
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
              out.reaction_terms[q][projected_density_idx] =
                out.densities[q] - in.composition[q][projected_density_idx];
            }
        }

      // Calculate ranges of reaction terms and rates
      if (reaction_rate_out != nullptr && in.n_evaluation_points() > 0)
        {
          double min_in_X        = in.composition[0][X_idx];
          double max_in_X        = in.composition[0][X_idx];
          double min_in_density  = in.composition[0][projected_density_idx];
          double max_in_density  = in.composition[0][projected_density_idx];
          double min_out_density = out.densities[0];
          double max_out_density = out.densities[0];
          double min_term        = out.reaction_terms[0][projected_density_idx];
          double max_term        = out.reaction_terms[0][projected_density_idx];
          double min_rate        = reaction_rate_out->reaction_rates[0][X_idx];
          double max_rate        = reaction_rate_out->reaction_rates[0][X_idx];

          for (unsigned int q = 1; q < in.n_evaluation_points(); ++q)
            {
              if (in.composition[q][X_idx] < min_in_X)
                min_in_X = in.composition[q][X_idx];
              if (in.composition[q][X_idx] > max_in_X)
                max_in_X = in.composition[q][X_idx];
              if (in.composition[q][projected_density_idx] < min_in_density)
                min_in_density = in.composition[q][projected_density_idx];
              if (in.composition[q][projected_density_idx] > max_in_density)
                max_in_density = in.composition[q][projected_density_idx];
              if (out.densities[q] < min_out_density)
                min_out_density = out.densities[q];
              if (out.densities[q] > max_out_density)
                max_out_density = out.densities[q];
              if (out.reaction_terms[q][projected_density_idx] < min_term)
                min_term = out.reaction_terms[q][projected_density_idx];
              if (out.reaction_terms[q][projected_density_idx] > max_term)
                max_term = out.reaction_terms[q][projected_density_idx];
              if (reaction_rate_out->reaction_rates[q][X_idx] < min_rate)
                min_rate = reaction_rate_out->reaction_rates[q][X_idx];
              if (reaction_rate_out->reaction_rates[q][X_idx] > max_rate)
                max_rate = reaction_rate_out->reaction_rates[q][X_idx];
            }
        }

      // Update projected density field
      if (prescribed_field_out != nullptr &&
          projected_density_idx != numbers::invalid_unsigned_int)
        {
          for (unsigned int i = 0; i < in.position.size(); ++i)
            prescribed_field_out
            ->prescribed_field_outputs[i][projected_density_idx] =
              out.densities[i];
        }
    }

    template <int dim>
    void
    KineticDrivingForce<dim>::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Kinetic driving force");
        {
          // Data profile
          prm.declare_entry(
            "Data directory",
            "./",
            Patterns::DirectoryName(),
            "The name of a directory where the model data is found. This path "
            "may either be absolute (if starting with '/') or relative to "
            "the current directory.");
          prm.declare_entry("Data file name",
                            "thermo_data.txt",
                            Patterns::Anything(),
                            "The file name of the model data.");

          // Rheological parameters
          prm.declare_entry("Viscosity",
                            "1e21",
                            Patterns::Double(0.),
                            "Viscosity. Units: Pa s");
          prm.declare_entry(
            "Thermal viscosity exponent",
            "0.0",
            Patterns::Double(0.),
            "The temperature dependence of viscosity. Dimensionless exponent.");
          prm.declare_entry(
            "Transition depths",
            "150e3, 410e3, 660e3",
            Patterns::List(Patterns::Double(0.)),
            "A list of depths where the viscosity changes. Values must "
            "monotonically increase. "
            "Units: m");
          prm.declare_entry(
            "Viscosity prefactors",
            "10.0, 0.1, 1.0, 10.0",
            Patterns::List(Patterns::Double(0.)),
            "A list of prefactors for the viscosity that determine the "
            "viscosity profile. Each prefactor is applied in a depth range "
            "specified by the list of `Transition depths', i.e. the first "
            "prefactor is applied above the first transition depth, the "
            "second one between the first and second transition depth, and "
            "so on. To compute the viscosity profile, this prefactor is "
            "multiplied by the reference viscosity specified through the "
            "parameter `Viscosity'. List must have one more entry than "
            "Transition depths. Units: non-dimensional.");
          prm.declare_entry("Minimum viscosity",
                            "1e19",
                            Patterns::Double(0),
                            "The minimum viscosity cutoff. Units: Pa s");
          prm.declare_entry("Maximum viscosity",
                            "1e24",
                            Patterns::Double(0),
                            "The maximum viscosity cutoff. Units: Pa s");

          // Material parameters
          prm.declare_entry("Thermal conductivity",
                            "4.0",
                            Patterns::Double(0.),
                            "Reference conductivity");
          prm.declare_entry("Use adiabatic pressure for density",
                            "true",
                            Patterns::Bool(),
                            "Whether to use adiabatic pressure instead of full "
                            "pressure in density calculation.");
        }

        // Kinetic parameters
        prm.declare_entry(
          "Kinetic prefactor Q",
          "1e-5",
          Patterns::Double(0.0),
          "The scalar kinetic prefactor Q in the reaction rate "
          "equation dX/dt = Q * driving_force * (1-X). Units: mol/J/s");
        prm.declare_entry("Use temperature cutoff",
                          "false",
                          Patterns::Bool(),
                          "Whether to apply a temperature cutoff below which "
                          "reactions are suppressed.");
        prm.declare_entry(
          "Reaction cutoff temperature",
          "0.0",
          Patterns::Double(0.0),
          "Temperature below which reactions are suppressed. Units: K");
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    KineticDrivingForce<dim>::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Kinetic driving force");
        {
          // Parse data profile
          const std::string data_directory = prm.get("Data directory");
          const std::string data_file_name = prm.get("Data file name");
          ParameterHandler  profile_prm;
          Utilities::AsciiDataProfile<dim>::declare_parameters(profile_prm,
                                                               data_directory,
                                                               data_file_name);
          profile.parse_parameters(profile_prm);

          // Parse rheological parameters
          viscosity = prm.get_double("Viscosity");
          thermal_viscosity_exponent =
            prm.get_double("Thermal viscosity exponent");
          transition_depths = Utilities::string_to_double(
                                Utilities::split_string_list(prm.get("Transition depths")));
          viscosity_prefactors = Utilities::string_to_double(
                                   Utilities::split_string_list(prm.get("Viscosity prefactors")));
          minimum_viscosity = prm.get_double("Minimum viscosity");
          maximum_viscosity = prm.get_double("Maximum viscosity");

          // Parse material parameters
          k = prm.get_double("Thermal conductivity");
          use_adiabatic_pressure_for_density =
            prm.get_bool("Use adiabatic pressure for density");

          // Parse kinetic parameters
          Q_kinetic_prefactor    = prm.get_double("Kinetic prefactor Q");
          use_temperature_cutoff = prm.get_bool("Use temperature cutoff");
          reaction_cutoff_temperature =
            prm.get_double("Reaction cutoff temperature");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // Validate parameters
      if (viscosity_prefactors.size() != transition_depths.size() + 1)
        AssertThrow(
          false,
          ExcMessage(
            "Error: The list of Viscosity prefactors needs to have exactly "
            "one more entry than the list of Transition depths. "));

      AssertThrow(Q_kinetic_prefactor > 0.0,
                  ExcMessage("Kinetic prefactor Q must be positive."));

      if (use_temperature_cutoff)
        {
          AssertThrow(reaction_cutoff_temperature >= 0.0,
                      ExcMessage(
                        "Reaction cutoff temperature must be non-negative."));
        }

      // Declare dependencies on solution variables
      this->model_dependence.viscosity = NonlinearDependence::none;
      this->model_dependence.density =
        NonlinearDependence::pressure | NonlinearDependence::temperature;
      this->model_dependence.compressibility      = NonlinearDependence::none;
      this->model_dependence.specific_heat        = NonlinearDependence::none;
      this->model_dependence.thermal_conductivity = NonlinearDependence::none;
    }

    template <int dim>
    void
    KineticDrivingForce<dim>::create_additional_named_outputs(
      MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      if (out.template get_additional_output<ReactionRateOutputs<dim>>() ==
          nullptr)
        {
          const unsigned int n_points = out.n_evaluation_points();
          out.additional_outputs.push_back(
            std::make_unique<MaterialModel::ReactionRateOutputs<dim>>(
              n_points, this->introspection().n_compositional_fields));
        }

      const unsigned int projected_density_idx =
        this->introspection().find_composition_type(
          CompositionalFieldDescription::density);

      if (projected_density_idx != numbers::invalid_unsigned_int &&
          out.template get_additional_output<PrescribedFieldOutputs<dim>>() ==
          nullptr)
        {
          const unsigned int n_points = out.n_evaluation_points();
          out.additional_outputs.push_back(
            std::make_unique<MaterialModel::PrescribedFieldOutputs<dim>>(
              n_points, this->introspection().n_compositional_fields));
        }
    }
  } // namespace MaterialModel
} // namespace aspect

namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(
      KineticDrivingForce,
      "kinetic driving force",
      "A material model that computes kinetic driving force using the "
      "equation dX/dt = Q * ΔG * (1-X), where Q is a user-defined kinetic "
      "prefactor and ΔG is the thermodynamic driving force calculated from "
      "equilibrium data. The driving force is computed as "
      "ΔG = ΔG₀ + (P - P_adiabatic)ΔV - (T - T_adiabatic)ΔS, where all "
      "thermodynamic parameters are read from a user-specified ASCII data "
      "file. Requires at least one compositional field representing the "
      "reacting phase. This model acts as a decorator around a base material "
      "model and supports projected density formulation.")
  } // namespace MaterialModel
} // namespace aspect
