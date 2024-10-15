/*
  Copyright (C) 2019 - 2024 by the authors of the ASPECT code.

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


#include <aspect/material_model/rheology/diffusion_creep.h>
#include <aspect/material_model/utilities.h>
#include <aspect/utilities.h>

#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/parameter_handler.h>


namespace aspect
{
  namespace MaterialModel
  {
    namespace Rheology
    {
      DiffusionCreepParameters::DiffusionCreepParameters()
        : prefactor (numbers::signaling_nan<double>()),
          activation_energy (numbers::signaling_nan<double>()),
          activation_volume (numbers::signaling_nan<double>()),
          stress_exponent (numbers::signaling_nan<double>()),
          grain_size_exponent (numbers::signaling_nan<double>())
      {}



      template <int dim>
      DiffusionCreep<dim>::DiffusionCreep ()
        = default;



      template <int dim>
      const DiffusionCreepParameters
      DiffusionCreep<dim>::compute_creep_parameters (const unsigned int composition,
                                                     const std::vector<double> &phase_function_values,
                                                     const std::vector<unsigned int> &n_phase_transitions_per_composition) const
      {
        DiffusionCreepParameters creep_parameters;
        if (phase_function_values == std::vector<double>())
          {
            // no phases
            creep_parameters.prefactor = prefactors_diffusion[composition];
            creep_parameters.activation_energy = activation_energies_diffusion[composition];
            creep_parameters.activation_volume = activation_volumes_diffusion[composition];
            creep_parameters.stress_exponent = stress_exponents_diffusion[composition];
            creep_parameters.grain_size_exponent = grain_size_exponents_diffusion[composition];
          }
        else
          {
            // Average among phases
            creep_parameters.prefactor = MaterialModel::MaterialUtilities::phase_average_value(phase_function_values, n_phase_transitions_per_composition,
                                         prefactors_diffusion, composition,  MaterialModel::MaterialUtilities::PhaseUtilities::logarithmic);
            creep_parameters.activation_energy = MaterialModel::MaterialUtilities::phase_average_value(phase_function_values, n_phase_transitions_per_composition,
                                                 activation_energies_diffusion, composition);
            creep_parameters.activation_volume = MaterialModel::MaterialUtilities::phase_average_value(phase_function_values, n_phase_transitions_per_composition,
                                                 activation_volumes_diffusion, composition);
            creep_parameters.stress_exponent = MaterialModel::MaterialUtilities::phase_average_value(phase_function_values, n_phase_transitions_per_composition,
                                               stress_exponents_diffusion, composition);
            creep_parameters.grain_size_exponent = MaterialModel::MaterialUtilities::phase_average_value(phase_function_values, n_phase_transitions_per_composition,
                                                   grain_size_exponents_diffusion, composition);
          }
        return creep_parameters;
      }



      template <int dim>
      double
      DiffusionCreep<dim>::compute_viscosity (const double pressure,
                                              const double temperature,
                                              const unsigned int composition,
                                              const std::vector<double> &phase_function_values,
                                              const std::vector<unsigned int> &n_phase_transitions_per_composition) const
      {
        return compute_viscosity(pressure, temperature, fixed_grain_size, composition, phase_function_values, n_phase_transitions_per_composition);
      }



      template <int dim>
      double
      DiffusionCreep<dim>::compute_viscosity (const double pressure,
                                              const double temperature,
                                              const double grain_size,
                                              const unsigned int composition,
                                              const std::vector<double> &phase_function_values,
                                              const std::vector<unsigned int> &n_phase_transitions_per_composition) const
      {
        const DiffusionCreepParameters p = compute_creep_parameters(composition,
                                                                    phase_function_values,
                                                                    n_phase_transitions_per_composition);

        // Power law creep equation
        //    viscosity = 0.5 * A^(-1) * d^(m) * exp((E + P*V)/(RT))
        // A: prefactor,
        // d: grain size, m: grain size exponent, E: activation energy, P: pressure,
        // V; activation volume, R: gas constant, T: temperature.
        double viscosity_diffusion = 0.5 / p.prefactor *
                                     std::exp((p.activation_energy +
                                               pressure*p.activation_volume)/
                                              (constants::gas_constant*temperature)) *
                                     std::pow(grain_size, p.grain_size_exponent);

        Assert (viscosity_diffusion > 0.0,
                ExcMessage ("Negative diffusion viscosity detected. This is unphysical and should not happen. "
                            "Check for negative parameters. Temperature and pressure are "
                            + Utilities::to_string(temperature) + " K, " + Utilities::to_string(pressure) + " Pa. "));

        // Creep viscosities become extremely large at low
        // temperatures and can therefore provoke floating-point overflow errors. In
        // real rocks, other deformation mechanisms become dominant at low temperatures,
        // so these high viscosities are never achieved. It is therefore both reasonable
        // and desirable to require the single-mechanism viscosity to be smaller than
        // std::sqrt(max_double).
        viscosity_diffusion = std::min(viscosity_diffusion, std::sqrt(std::numeric_limits<double>::max()));

        return viscosity_diffusion;
      }



      template <int dim>
      std::pair<double, double>
      DiffusionCreep<dim>::compute_strain_rate_and_derivative (const double stress,
                                                               const double pressure,
                                                               const double temperature,
                                                               const DiffusionCreepParameters creep_parameters) const
      {
        return compute_strain_rate_and_derivative(stress, pressure, temperature, fixed_grain_size, creep_parameters);
      }



      template <int dim>
      std::pair<double, double>
      DiffusionCreep<dim>::compute_strain_rate_and_derivative(const double stress,
                                                              const double pressure,
                                                              const double temperature,
                                                              const double grain_size,
                                                              const DiffusionCreepParameters creep_parameters) const
      {
        // Power law creep equation
        //   edot_ii_partial = A * stress^n * d^-m * exp(-(E + P*V)/(RT))
        //   d(edot_ii_partial)/d(stress) = A * n * stress^(n-1) * d^-m * exp(-(E + P*V)/(RT))
        // A: prefactor, edot_ii_partial: square root of second invariant of deviatoric strain rate tensor attributable to the creep mechanism,
        // d: grain size, m: grain size exponent, E: activation energy, P: pressure,
        // V; activation volume, R: gas constant, T: temperature.
        // For diffusion creep, n = 1 (strain rate is linearly dependent on stress).
        const double dstrain_rate_dstress_diffusion = creep_parameters.prefactor *
                                                      std::pow(grain_size, -creep_parameters.grain_size_exponent) *
                                                      std::exp(-(creep_parameters.activation_energy + pressure * creep_parameters.activation_volume) /
                                                               (constants::gas_constant * temperature));

        const double strain_rate_diffusion = stress * dstrain_rate_dstress_diffusion;

        return std::make_pair(strain_rate_diffusion, dstrain_rate_dstress_diffusion);
      }



      template <int dim>
      std::pair<double, double>
      DiffusionCreep<dim>::compute_log_strain_rate_and_derivative (const double log_stress,
                                                                   const double pressure,
                                                                   const double temperature,
                                                                   const DiffusionCreepParameters creep_parameters) const
      {
        return compute_log_strain_rate_and_derivative (log_stress, pressure, temperature, fixed_grain_size, creep_parameters);
      }



      template <int dim>
      std::pair<double, double>
      DiffusionCreep<dim>::compute_log_strain_rate_and_derivative (const double log_stress,
                                                                   const double pressure,
                                                                   const double temperature,
                                                                   const double grain_size,
                                                                   const DiffusionCreepParameters creep_parameters) const
      {
        // Power law creep equation
        // log(edot_ii_partial) = std::log(A) + n*std::log(stress) - m*std::log(d) - (E + P*V)/(RT)
        //   d(log_edot_ii_partial)/d(log_stress) = n
        // A: prefactor, edot_ii_partial: square root of second invariant of deviatoric strain rate tensor attributable to the creep mechanism,
        // d: grain size, m: grain size exponent, E: activation energy, P: pressure,
        // V; activation volume, R: gas constant, T: temperature.
        // For diffusion creep, n = 1 (strain rate is linearly dependent on stress).
        const double log_strain_rate_diffusion = std::log(creep_parameters.prefactor) +
                                                 log_stress -
                                                 creep_parameters.grain_size_exponent * std::log(grain_size) -
                                                 (creep_parameters.activation_energy + pressure*creep_parameters.activation_volume)/
                                                 (constants::gas_constant*temperature);

        const double dlog_strain_rate_dlog_stress_diffusion = 1.0;

        return std::make_pair(log_strain_rate_diffusion, dlog_strain_rate_dlog_stress_diffusion);
      }



      template <int dim>
      void
      DiffusionCreep<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.declare_entry ("Prefactors for diffusion creep", "1.5e-15",
                           Patterns::Anything(),
                           "List of viscosity prefactors, $A$, for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of all compositional fields or only "
                           "those corresponding to chemical compositions. "
                           "If only one value is given, then all use the same value. "
                           "Units: \\si{\\per\\pascal\\meter}$^{m_{\\text{diffusion}}}$\\si{\\per\\second}.");
        prm.declare_entry ("Stress exponents for diffusion creep", "1.",
                           Patterns::List(Patterns::Double(0.)),
                           "List of stress exponents, $n_{\\text{diffusion}}$, for background mantle and compositional fields, "
                           "for a total of N+1 values, where N is the number of all compositional fields or only "
                           "those corresponding to chemical compositions. "
                           "The stress exponent for diffusion creep is almost always equal to one. "
                           "If only one value is given, then all use the same value.  Units: None.");
        prm.declare_entry ("Grain size exponents for diffusion creep", "3.",
                           Patterns::Anything(),
                           "List of grain size exponents, $m_{\\text{diffusion}}$, for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of all compositional fields or only "
                           "those corresponding to chemical compositions. "
                           "If only one value is given, then all use the same value. Units: None.");
        prm.declare_entry ("Activation energies for diffusion creep", "375e3",
                           Patterns::Anything(),
                           "List of activation energies, $E_a$, for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of all compositional fields or only "
                           "those corresponding to chemical compositions. "
                           "If only one value is given, then all use the same value. "
                           "Units: \\si{\\joule\\per\\mole}.");
        prm.declare_entry ("Activation volumes for diffusion creep", "6e-6",
                           Patterns::Anything(),
                           "List of activation volumes, $V_a$, for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of all compositional fields or only "
                           "those corresponding to chemical compositions. "
                           "If only one value is given, then all use the same value. "
                           "Units: \\si{\\meter\\cubed\\per\\mole}.");
        prm.declare_entry ("Grain size", "1e-3", Patterns::Double (0.),
                           "The fixed grain size of the material. "
                           "This grain size is only used if the parent "
                           "material model does not provide its own "
                           "(possibly variable) grain size when "
                           "calling this rheology."
                           "Units: \\si{\\meter}.");
      }



      template <int dim>
      void
      DiffusionCreep<dim>::parse_parameters (ParameterHandler &prm,
                                             const std::unique_ptr<std::vector<unsigned int>> &expected_n_phases_per_composition)
      {
        // Retrieve the list of composition names
        std::vector<std::string> compositional_field_names = this->introspection().get_composition_names();

        // Retrieve the list of names of fields that represent chemical compositions, and not, e.g.,
        // plastic strain
        std::vector<std::string> chemical_field_names = this->introspection().chemical_composition_field_names();

        // Establish that a background field is required here
        compositional_field_names.insert(compositional_field_names.begin(), "background");
        chemical_field_names.insert(chemical_field_names.begin(), "background");

        // Make options file for parsing maps to double arrays
        Utilities::MapParsing::Options options(chemical_field_names, "Prefactors for diffusion creep");
        options.list_of_allowed_keys = compositional_field_names;
        options.allow_multiple_values_per_key = true;
        if (expected_n_phases_per_composition)
          {
            options.n_values_per_key = *expected_n_phases_per_composition;

            // check_values_per_key is required to be true to duplicate single values
            // if they are to be used for all phases associated with a given key.
            options.check_values_per_key = true;
          }

        // Read parameters, each of size of number of composition + number of phases + 1
        prefactors_diffusion = Utilities::MapParsing::parse_map_to_double_array(prm.get("Prefactors for diffusion creep"),
                                                                                options);

        options.property_name = "Stress exponents for diffusion creep";
        stress_exponents_diffusion = Utilities::MapParsing::parse_map_to_double_array(prm.get("Stress exponents for diffusion creep"),
                                                                                      options);

        options.property_name = "Grain size exponents for diffusion creep";
        grain_size_exponents_diffusion = Utilities::MapParsing::parse_map_to_double_array(prm.get("Grain size exponents for diffusion creep"),
                                         options);

        options.property_name = "Activation energies for diffusion creep";
        activation_energies_diffusion = Utilities::MapParsing::parse_map_to_double_array(prm.get("Activation energies for diffusion creep"),
                                        options);

        options.property_name = "Activation volumes for diffusion creep";
        activation_volumes_diffusion = Utilities::MapParsing::parse_map_to_double_array(prm.get("Activation volumes for diffusion creep"),
                                                                                        options);

        fixed_grain_size = prm.get_double("Grain size");

        // Check that there are no entries set to zero,
        // for example because the entry is for a field
        // that is masked anyway, like strain. Despite
        // these compositions being masked, their viscosities
        // are computed anyway and this will lead to division by zero.
        for (const double prefactor : prefactors_diffusion)
          AssertThrow(prefactor > 0.,
                      ExcMessage("The diffusion prefactor should be larger than zero."));
      }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
#define INSTANTIATE(dim) \
  namespace Rheology \
  { \
    template class DiffusionCreep<dim>; \
  }

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
