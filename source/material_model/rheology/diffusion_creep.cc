/*
  Copyright (C) 2019 - 2020 by the authors of the ASPECT code.

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
      template <int dim>
      DiffusionCreep<dim>::DiffusionCreep ()
      {}



      template <int dim>
      const DiffusionCreepParameters
      DiffusionCreep<dim>::compute_creep_parameters (const unsigned int composition,
                                                     const std::vector<double> &phase_function_values,
                                                     const std::vector<unsigned int> &n_phases_per_composition) const
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
            creep_parameters.prefactor = MaterialModel::MaterialUtilities::phase_average_value(phase_function_values, n_phases_per_composition,
                                         prefactors_diffusion, composition,  MaterialModel::MaterialUtilities::PhaseUtilities::logarithmic);
            creep_parameters.activation_energy = MaterialModel::MaterialUtilities::phase_average_value(phase_function_values, n_phases_per_composition,
                                                 activation_energies_diffusion, composition);
            creep_parameters.activation_volume = MaterialModel::MaterialUtilities::phase_average_value(phase_function_values, n_phases_per_composition,
                                                 activation_volumes_diffusion, composition);
            creep_parameters.stress_exponent = MaterialModel::MaterialUtilities::phase_average_value(phase_function_values, n_phases_per_composition,
                                               stress_exponents_diffusion, composition);
            creep_parameters.grain_size_exponent = MaterialModel::MaterialUtilities::phase_average_value(phase_function_values, n_phases_per_composition,
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
                                              const std::vector<unsigned int> &n_phases_per_composition) const
      {
        const DiffusionCreepParameters p = compute_creep_parameters(composition,
                                                                    phase_function_values,
                                                                    n_phases_per_composition);

        // Power law creep equation
        //    viscosity = 0.5 * A^(-1/n) * d^(m/n) * exp((E + P*V)/(nRT))
        // A: prefactor,
        // d: grain size, m: grain size exponent, E: activation energy, P: pressure,
        // V; activation volume, R: gas constant, T: temperature.
        const double viscosity_diffusion = 0.5 / p.prefactor *
                                           std::exp((p.activation_energy +
                                                     pressure*p.activation_volume)/
                                                    (constants::gas_constant*temperature)) *
                                           std::pow(grain_size, p.grain_size_exponent);

        return viscosity_diffusion;
      }



      template <int dim>
      std::pair<double, double>
      DiffusionCreep<dim>::compute_strain_rate_and_derivative (const double stress,
                                                               const double pressure,
                                                               const double temperature,
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
                                                      std::exp(-(creep_parameters.activation_energy + pressure*creep_parameters.activation_volume)/
                                                               (constants::gas_constant*temperature));

        const double strain_rate_diffusion = stress * dstrain_rate_dstress_diffusion;

        return std::make_pair(strain_rate_diffusion, dstrain_rate_dstress_diffusion);
      }



      template <int dim>
      void
      DiffusionCreep<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.declare_entry ("Prefactors for diffusion creep", "1.5e-15",
                           Patterns::Anything(),
                           "List of viscosity prefactors, $A$, for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of compositional fields. "
                           "If only one value is given, then all use the same value. "
                           "Units: \\si{\\per\\pascal\\meter}$^{m_{\\text{diffusion}}}$\\si{\\per\\second}.");
        prm.declare_entry ("Stress exponents for diffusion creep", "1.",
                           Patterns::List(Patterns::Double(0.)),
                           "List of stress exponents, $n_{\\text{diffusion}}$, for background mantle and compositional fields, "
                           "for a total of N+1 values, where N is the number of compositional fields. "
                           "The stress exponent for diffusion creep is almost always equal to one. "
                           "If only one value is given, then all use the same value.  Units: None.");
        prm.declare_entry ("Grain size exponents for diffusion creep", "3.",
                           Patterns::Anything(),
                           "List of grain size exponents, $m_{\\text{diffusion}}$, for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of compositional fields. "
                           "If only one value is given, then all use the same value. Units: None.");
        prm.declare_entry ("Activation energies for diffusion creep", "375e3",
                           Patterns::Anything(),
                           "List of activation energies, $E_a$, for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of compositional fields. "
                           "If only one value is given, then all use the same value. "
                           "Units: \\si{\\joule\\per\\mole}.");
        prm.declare_entry ("Activation volumes for diffusion creep", "6e-6",
                           Patterns::Anything(),
                           "List of activation volumes, $V_a$, for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of compositional fields. "
                           "If only one value is given, then all use the same value. "
                           "Units: \\si{\\meter\\cubed\\per\\mole}.");
        prm.declare_entry ("Grain size", "1e-3", Patterns::Double (0.),
                           "Units: \\si{\\meter}.");
      }



      template <int dim>
      void
      DiffusionCreep<dim>::parse_parameters (ParameterHandler &prm,
                                             const std::shared_ptr<std::vector<unsigned int>> &expected_n_phases_per_composition)
      {
        // Retrieve the list of composition names
        const std::vector<std::string> list_of_composition_names = this->introspection().get_composition_names();

        // Establish that a background field is required here
        const bool has_background_field = true;

        // Read parameters, each of size of number of composition + number of phases + 1
        prefactors_diffusion = Utilities::parse_map_to_double_array(prm.get("Prefactors for diffusion creep"),
                                                                    list_of_composition_names,
                                                                    has_background_field,
                                                                    "Prefactors for diffusion creep",
                                                                    true,
                                                                    expected_n_phases_per_composition);
        stress_exponents_diffusion = Utilities::parse_map_to_double_array(prm.get("Stress exponents for diffusion creep"),
                                                                          list_of_composition_names,
                                                                          has_background_field,
                                                                          "Prefactors for diffusion creep",
                                                                          true,
                                                                          expected_n_phases_per_composition);
        grain_size_exponents_diffusion = Utilities::parse_map_to_double_array(prm.get("Grain size exponents for diffusion creep"),
                                                                              list_of_composition_names,
                                                                              has_background_field,
                                                                              "Grain size exponents for diffusion creep",
                                                                              true,
                                                                              expected_n_phases_per_composition);
        activation_energies_diffusion = Utilities::parse_map_to_double_array(prm.get("Activation energies for diffusion creep"),
                                                                             list_of_composition_names,
                                                                             has_background_field,
                                                                             "Activation energies for diffusion creep",
                                                                             true,
                                                                             expected_n_phases_per_composition);
        activation_volumes_diffusion = Utilities::parse_map_to_double_array(prm.get("Activation volumes for diffusion creep"),
                                                                            list_of_composition_names,
                                                                            has_background_field,
                                                                            "Activation volumes for diffusion creep",
                                                                            true,
                                                                            expected_n_phases_per_composition);
        grain_size = prm.get_double("Grain size");
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
