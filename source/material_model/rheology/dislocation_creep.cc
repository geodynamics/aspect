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


#include <aspect/material_model/rheology/dislocation_creep.h>
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
      DislocationCreep<dim>::DislocationCreep ()
      {}



      template <int dim>
      const DislocationCreepParameters
      DislocationCreep<dim>::compute_creep_parameters (const unsigned int composition,
                                                       const std::vector<double> &phase_function_values,
                                                       const std::vector<unsigned int> &n_phases_per_composition) const
      {
        DislocationCreepParameters creep_parameters;
        if (phase_function_values == std::vector<double>())
          {
            // no phases
            creep_parameters.prefactor = prefactors_dislocation[composition];
            creep_parameters.activation_energy = activation_energies_dislocation[composition];
            creep_parameters.activation_volume = activation_volumes_dislocation[composition];
            creep_parameters.stress_exponent = stress_exponents_dislocation[composition];
          }
        else
          {
            // Average among phases
            creep_parameters.prefactor = MaterialModel::MaterialUtilities::phase_average_value(phase_function_values, n_phases_per_composition,
                                         prefactors_dislocation, composition,
                                         MaterialModel::MaterialUtilities::PhaseUtilities::logarithmic);
            creep_parameters.activation_energy = MaterialModel::MaterialUtilities::phase_average_value(phase_function_values, n_phases_per_composition,
                                                 activation_energies_dislocation, composition);
            creep_parameters.activation_volume = MaterialModel::MaterialUtilities::phase_average_value(phase_function_values, n_phases_per_composition,
                                                 activation_volumes_dislocation , composition);
            creep_parameters.stress_exponent = MaterialModel::MaterialUtilities::phase_average_value(phase_function_values, n_phases_per_composition,
                                               stress_exponents_dislocation, composition);
          }

        return creep_parameters;
      }



      template <int dim>
      double
      DislocationCreep<dim>::compute_viscosity (const double strain_rate,
                                                const double pressure,
                                                const double temperature,
                                                const unsigned int composition,
                                                const std::vector<double> &phase_function_values,
                                                const std::vector<unsigned int> &n_phases_per_composition) const
      {
        const DislocationCreepParameters p = compute_creep_parameters(composition,
                                                                      phase_function_values,
                                                                      n_phases_per_composition);

        // Power law creep equation:
        //    viscosity = 0.5 * A^(-1/n) * edot_ii^((1-n)/n) * exp((E + P*V)/(nRT))
        // A: prefactor, edot_ii: square root of second invariant of deviatoric strain rate tensor,
        // E: activation energy, P: pressure,
        // V; activation volume, n: stress exponent, R: gas constant, T: temperature.
        const double viscosity_dislocation = 0.5 * std::pow(p.prefactor,-1/p.stress_exponent) *
                                             std::exp((p.activation_energy + pressure*p.activation_volume)/
                                                      (constants::gas_constant*temperature*p.stress_exponent)) *
                                             std::pow(strain_rate,((1. - p.stress_exponent)/p.stress_exponent));
        return viscosity_dislocation;
      }



      template <int dim>
      std::pair<double, double>
      DislocationCreep<dim>::compute_strain_rate_and_derivative (const double stress,
                                                                 const double pressure,
                                                                 const double temperature,
                                                                 const DislocationCreepParameters creep_parameters) const
      {
        // Power law creep equation:
        //   edot_ii_partial = A * stress^n * exp(-(E + P*V)/(RT))
        //   d(edot_ii_partial)/d(stress) = A * n * stress^(n-1) * exp(-(E + P*V)/(RT))
        // A: prefactor, edot_ii_partial: square root of second invariant of deviatoric strain rate tensor attributable to the creep mechanism,
        // stress: deviatoric stress, E: activation energy, P: pressure,
        // V; activation volume, n: stress exponent, R: gas constant, T: temperature.
        const double strain_rate_dislocation = creep_parameters.prefactor *
                                               std::pow(stress,creep_parameters.stress_exponent) *
                                               std::exp(-(creep_parameters.activation_energy + pressure*creep_parameters.activation_volume)/
                                                        (constants::gas_constant*temperature));

        const double dstrain_rate_dstress_dislocation = creep_parameters.prefactor *
                                                        creep_parameters.stress_exponent *
                                                        std::pow(stress,creep_parameters.stress_exponent-1.) *
                                                        std::exp(-(creep_parameters.activation_energy + pressure*creep_parameters.activation_volume)/
                                                                 (constants::gas_constant*temperature));

        return std::make_pair(strain_rate_dislocation, dstrain_rate_dstress_dislocation);
      }



      template <int dim>
      void
      DislocationCreep<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.declare_entry ("Prefactors for dislocation creep", "1.1e-16",
                           Patterns::Anything(),
                           "List of viscosity prefactors, $A$, for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of compositional fields. "
                           "If only one value is given, then all use the same value. "
                           "Units: \\si{\\pascal}$^{-n_{\\text{dislocation}}}$ \\si{\\per\\second}.");
        prm.declare_entry ("Stress exponents for dislocation creep", "3.5",
                           Patterns::Anything(),
                           "List of stress exponents, $n_{\\text{dislocation}}$, for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of compositional fields. "
                           "If only one value is given, then all use the same value.  Units: None.");
        prm.declare_entry ("Activation energies for dislocation creep", "530e3",
                           Patterns::Anything(),
                           "List of activation energies, $E_a$, for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of compositional fields. "
                           "If only one value is given, then all use the same value. "
                           "Units: \\si{\\joule\\per\\mole}.");
        prm.declare_entry ("Activation volumes for dislocation creep", "1.4e-5",
                           Patterns::Anything(),
                           "List of activation volumes, $V_a$, for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of compositional fields. "
                           "If only one value is given, then all use the same value. "
                           "Units: \\si{\\meter\\cubed\\per\\mole}.");
      }



      template <int dim>
      void
      DislocationCreep<dim>::parse_parameters (ParameterHandler &prm,
                                               const std::shared_ptr<std::vector<unsigned int>> &expected_n_phases_per_composition)
      {
        // Retrieve the list of composition names
        const std::vector<std::string> list_of_composition_names = this->introspection().get_composition_names();

        const bool has_background_field = true;

        // Read parameters, each of size of number of composition + number of phases + 1
        prefactors_dislocation = Utilities::parse_map_to_double_array(prm.get("Prefactors for dislocation creep"),
                                                                      list_of_composition_names,
                                                                      has_background_field,
                                                                      "Prefactors for dislocation creep",
                                                                      true,
                                                                      expected_n_phases_per_composition);

        stress_exponents_dislocation = Utilities::parse_map_to_double_array(prm.get("Stress exponents for dislocation creep"),
                                                                            list_of_composition_names,
                                                                            has_background_field,
                                                                            "Stress exponents for dislocation creep",
                                                                            true,
                                                                            expected_n_phases_per_composition);

        activation_energies_dislocation = Utilities::parse_map_to_double_array(prm.get("Activation energies for dislocation creep"),
                                                                               list_of_composition_names,
                                                                               has_background_field,
                                                                               "Activation energies for dislocation creep",
                                                                               true,
                                                                               expected_n_phases_per_composition);
        activation_volumes_dislocation  = Utilities::parse_map_to_double_array(prm.get("Activation volumes for dislocation creep"),
                                                                               list_of_composition_names,
                                                                               has_background_field,
                                                                               "Activation volumes for dislocation creep",
                                                                               true,
                                                                               expected_n_phases_per_composition);
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
    template class DislocationCreep<dim>; \
  }

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
