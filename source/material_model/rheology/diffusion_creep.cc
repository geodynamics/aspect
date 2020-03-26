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
      double
      DiffusionCreep<dim>::compute_viscosity (const double pressure,
                                              const double temperature,
                                              const unsigned int composition,
                                              const std::pair<std::vector<double>, const std::vector<unsigned int>> &gamma_inputs) const
      {
        // If there is phase chage
        //  Phase change: average flow lay parameters from phase, where gamma_inputs is not nullptr
        //  No phase change: take flow lay parameters from composition
        double prefactors_diffusion_c;
        double activation_energies_diffusion_c;
        double activation_volumes_diffusion_c;
        double grain_size_exponents_diffusion_c;
        double grain_size_c;
        if (gamma_inputs != std::pair<std::vector<double>, const std::vector<unsigned int>>())
          {
            // Log values of prefactors are weighted by gamma values and then converted to a value by taking the exponetial
            // in order to change values of viscosity by magnitude with gamma values.
            prefactors_diffusion_c = MaterialModel::MaterialUtilities::phase_average_value(gamma_inputs, prefactors_diffusion, composition, true);
            activation_energies_diffusion_c = MaterialModel::MaterialUtilities::phase_average_value(gamma_inputs, activation_energies_diffusion, composition);
            activation_volumes_diffusion_c = MaterialModel::MaterialUtilities::phase_average_value(gamma_inputs, activation_volumes_diffusion , composition);
            grain_size_exponents_diffusion_c = MaterialModel::MaterialUtilities::phase_average_value(gamma_inputs, grain_size_exponents_diffusion, composition);
            grain_size_c = MaterialModel::MaterialUtilities::phase_average_value(gamma_inputs, grain_size, composition);
          }
        else
          {
            prefactors_diffusion_c = prefactors_diffusion[composition];
            activation_energies_diffusion_c = activation_energies_diffusion[composition];
            activation_volumes_diffusion_c = activation_volumes_diffusion[composition];
            grain_size_exponents_diffusion_c = grain_size_exponents_diffusion[composition];
            grain_size_c = grain_size[composition];
          }
        // Power law creep equation
        //    viscosity = 0.5 * A^(-1/n) * d^(m/n) * exp((E + P*V)/(nRT))
        // A: prefactor,
        // d: grain size, m: grain size exponent, E: activation energy, P: pressure,
        // V; activation volume, R: gas constant, T: temperature.
        const double viscosity_diffusion = 0.5 / prefactors_diffusion_c *
                                           std::exp((activation_energies_diffusion_c +
                                                     pressure*activation_volumes_diffusion_c)/
                                                    (constants::gas_constant*temperature)) *
                                           std::pow(grain_size_c, grain_size_exponents_diffusion_c);

        return viscosity_diffusion;
      }



      template <int dim>
      void
      DiffusionCreep<dim>::declare_parameters (ParameterHandler &prm)
      {
        // Here, a default value is assigned to every phase and composition if no entry is given
        prm.declare_entry ("Prefactors for diffusion creep", "1.5e-15",
                           Patterns::Anything(),
                           "List of viscosity prefactors, $A$, for background material and compositional fields, "
                           "for a total of N+M+1 values, where N is the number of compositional fields and M is the number of phases. "
                           "If only one value is given, then all use the same value. "
                           "Units: $Pa^{-1} m^{m_{\\text{diffusion}}} s^{-1}$");
        prm.declare_entry ("Grain size exponents for diffusion creep", "3",
                           Patterns::Anything(),
                           "List of grain size exponents, $m_{\\text{diffusion}}$, for background material and compositional fields, "
                           "for a total of N+M+1 values, where N is the number of compositional fields and M is the number of phases. "
                           "If only one value is given, then all use the same value. Units: None");
        prm.declare_entry ("Activation energies for diffusion creep", "375e3",
                           Patterns::Anything(),
                           "List of activation energies, $E_a$, for background material and compositional fields, "
                           "for a total of N+M+1 values, where N is the number of compositional fields and M is the number of phases. "
                           "If only one value is given, then all use the same value.  Units: $J / mol$");
        prm.declare_entry ("Activation volumes for diffusion creep", "6e-6",
                           Patterns::Anything(),
                           "List of activation volumes, $V_a$, for background material and compositional fields, "
                           "for a total of N+M+1 values, where N is the number of compositional fields and M is the number of phases. "
                           "If only one value is given, then all use the same value.  Units: $m^3 / mol$");
        prm.declare_entry ("Grain size", "1e-3",
                           Patterns::Anything(),
                           "List of grain sizes, $V_a$, for background material and compositional fields, "
                           "for a total of N+M+1 values, where N is the number of compositional fields and M is the number of phases. "
                           "If only one value is given, then all use the same value.  Units: $m$");
      }



      template <int dim>
      void
      DiffusionCreep<dim>::parse_parameters (ParameterHandler &prm,
                                             const std::shared_ptr<std::vector<unsigned int>> expected_n_phases_per_composition)
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
        grain_size = Utilities::parse_map_to_double_array(prm.get("Grain size"),
                                                          list_of_composition_names,
                                                          has_background_field,
                                                          "Grain_size",
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
    template class DiffusionCreep<dim>; \
  }

    ASPECT_INSTANTIATE(INSTANTIATE)
  }
}
