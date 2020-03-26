/*
  Copyright (C) 2019 by the authors of the ASPECT code.

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
      double
      DislocationCreep<dim>::compute_viscosity (const double strain_rate,
                                                const double pressure,
                                                const double temperature,
                                                const unsigned int composition,
                                                const std::pair<std::vector<double>, const std::vector<unsigned int>> &gamma_inputs) const
      {
        // If there is phase chage
        //  Phase change: average flow lay parameters from phase, where gamma_inputs is not nullptr
        //  No phase change: take flow lay parameters from composition
        double prefactors_dislocation_c;
        double activation_energies_dislocation_c;
        double activation_volumes_dislocation_c;
        double stress_exponents_dislocation_c;
        if (gamma_inputs != std::pair<std::vector<double>, const std::vector<unsigned int>>())
          {
            // Log values of prefactors are weighted by gamma values and then converted to a value by taking the exponetial
            // in order to change values of viscosity by magnitude with gamma values.
            prefactors_dislocation_c = MaterialModel::MaterialUtilities::phase_average_value(gamma_inputs, prefactors_dislocation, composition, true);
            activation_energies_dislocation_c = MaterialModel::MaterialUtilities::phase_average_value(gamma_inputs, activation_energies_dislocation, composition);
            activation_volumes_dislocation_c = MaterialModel::MaterialUtilities::phase_average_value(gamma_inputs, activation_volumes_dislocation , composition);
            stress_exponents_dislocation_c = MaterialModel::MaterialUtilities::phase_average_value(gamma_inputs, stress_exponents_dislocation, composition);
          }
        else
          {
            prefactors_dislocation_c = prefactors_dislocation[composition];
            activation_energies_dislocation_c = activation_energies_dislocation[composition];
            activation_volumes_dislocation_c = activation_volumes_dislocation[composition];
            stress_exponents_dislocation_c = stress_exponents_dislocation[composition];
          }

        // Power law creep equation:
        //    viscosity = 0.5 * A^(-1/n) * edot_ii^((1-n)/n) * exp((E + P*V)/(nRT))
        // A: prefactor, edot_ii: square root of second invariant of deviatoric strain rate tensor,
        // E: activation energy, P: pressure,
        // V; activation volume, n: stress exponent, R: gas constant, T: temperature.
        const double viscosity_dislocation = 0.5 * std::pow(prefactors_dislocation_c,-1/stress_exponents_dislocation_c) *
                                             std::exp((activation_energies_dislocation_c + pressure*activation_volumes_dislocation_c)/
                                                      (constants::gas_constant*temperature*stress_exponents_dislocation_c)) *
                                             std::pow(strain_rate,((1. - stress_exponents_dislocation_c)/stress_exponents_dislocation_c));
        return viscosity_dislocation;
      }



      template <int dim>
      void
      DislocationCreep<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.declare_entry ("Prefactors for dislocation creep", "1.1e-16",
                           Patterns::Anything(),
                           "List of viscosity prefactors, $A$, for background material and compositional fields, "
                           "for a total of N+M+1 values, where N is the number of compositional fields and M is the number of phases. "
                           "If only one value is given, then all use the same value. "
                           "Units: $Pa^{-n_{\\text{dislocation}}} s^{-1}$");
        prm.declare_entry ("Stress exponents for dislocation creep", "3.5",
                           Patterns::Anything(),
                           "List of stress exponents, $n_{\\text{dislocation}}$, for background material and compositional fields, "
                           "for a total of N+M+1 values, where N is the number of compositional fields and M is the number of phases. "
                           "If only one value is given, then all use the same value.  Units: None");
        prm.declare_entry ("Activation energies for dislocation creep", "530e3",
                           Patterns::Anything(),
                           "List of activation energies, $E_a$, for background material and compositional fields, "
                           "for a total of N+M+1 values, where N is the number of compositional fields and M is the number of phases. "
                           "If only one value is given, then all use the same value.  Units: $J / mol$");
        prm.declare_entry ("Activation volumes for dislocation creep", "1.4e-5",
                           Patterns::Anything(),
                           "List of activation volumes, $V_a$, for background material and compositional fields, "
                           "for a total of N+M+1 values, where N is the number of compositional fields and M is the number of phases. "
                           "If only one value is given, then all use the same value.  Units: $m^3 / mol$");
      }



      template <int dim>
      void
      DislocationCreep<dim>::parse_parameters (ParameterHandler &prm, const std::shared_ptr<std::vector<unsigned int>> expected_n_phases_per_composition)
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
  }
}
