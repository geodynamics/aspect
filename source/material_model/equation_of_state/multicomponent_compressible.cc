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


#include <aspect/material_model/equation_of_state/multicomponent_compressible.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/utilities.h>


namespace aspect
{
  namespace MaterialModel
  {
    namespace EquationOfState
    {
      template <int dim>
      void
      MulticomponentCompressible<dim>::
      evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
               const unsigned int q,
               MaterialModel::EquationOfStateOutputs<dim> &eos_outputs) const
      {
        const double temperature = std::max(in.temperature[q], 1.); // temperature can't be zero for correct evaluation

        // If we are using the projected density approximation then we need to use the adiabatic pressure
        double pressure_for_density;
        if (this->introspection().composition_type_exists(Parameters<dim>::CompositionalFieldDescription::density) &&
            this->get_parameters().formulation_mass_conservation==Parameters<dim>::Formulation::MassConservation::projected_density_field)
          pressure_for_density = this->get_adiabatic_conditions().pressure(in.position[q]);
        else
          pressure_for_density = in.pressure[q];

        for (unsigned int c=0; c < eos_outputs.densities.size(); ++c)
          {
            const double ak = reference_thermal_expansivities[c]/reference_isothermal_compressibilities[c];
            const double f = (1. + (pressure_for_density - ak*(temperature - reference_temperatures[c])) *
                              isothermal_bulk_modulus_pressure_derivatives[c] *
                              reference_isothermal_compressibilities[c]);

            eos_outputs.densities[c] = reference_densities[c]*std::pow(f, 1./isothermal_bulk_modulus_pressure_derivatives[c]);
            eos_outputs.thermal_expansion_coefficients[c] = reference_thermal_expansivities[c] / f;
            eos_outputs.specific_heat_capacities[c] = (isochoric_specific_heats[c] +
                                                       (temperature*reference_thermal_expansivities[c] *
                                                        ak * std::pow(f, -1.-(1./isothermal_bulk_modulus_pressure_derivatives[c]))
                                                        / reference_densities[c]));
            eos_outputs.compressibilities[c] = reference_isothermal_compressibilities[c]/f;
            eos_outputs.entropy_derivative_pressure[c] = 0.;
            eos_outputs.entropy_derivative_temperature[c] = 0.;
          }
      }



      template <int dim>
      bool
      MulticomponentCompressible<dim>::
      is_compressible () const
      {
        return true;
      }



      template <int dim>
      void
      MulticomponentCompressible<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.declare_entry ("Reference temperatures", "298.15",
                           Patterns::Anything(),
                           "List of reference temperatures $T_0$ for background and compositional fields (N), "
                           "for a total of N+1 values for models with no phase transitions (or models where the "
                           "value does not change across any of the phase transitions). For models with phase "
                           "transitions, the list needs to contain each field name, including the background, for "
                           "a total of N+1 names, and for each of these names, specify the value for each phase "
                           "(giving P_c+1 values for each field, with P_c being the number of phase transitions "
                           "for field c). Therefore, the total number of values given is N+P+1, with P = sum(P_c) "
                           "the total number of phase transitions, summed over all phases. The format is "
                           "background: value1|value2|...|valueP_1+1, field1:value1|...|valueP_2+1, ..., fieldN: value1|...|valueP_N+1. "
                           "If only one value is given, then all fields/phases use the same value. "
                           "Units: \\si{\\kelvin}.");
        prm.declare_entry ("Reference densities", "3300.",
                           Patterns::Anything(),
                           "List of reference densities $T_0$ for background and compositional fields (N), "
                           "for a total of N+1 values for models with no phase transitions (or models where the "
                           "value does not change across any of the phase transitions). For models with phase "
                           "transitions, the list needs to contain each field name, including the background, for "
                           "a total of N+1 names, and for each of these names, specify the value for each phase "
                           "(giving P_c+1 values for each field, with P_c being the number of phase transitions "
                           "for field c). Therefore, the total number of values given is N+P+1, with P = sum(P_c) "
                           "the total number of phase transitions, summed over all phases. The format is "
                           "background: value1|value2|...|valueP_1+1, field1:value1|...|valueP_2+1, ..., fieldN: value1|...|valueP_N+1. "
                           "If only one value is given, then all fields/phases use the same value. "
                           "Units: \\si{\\kilogram\\per\\meter\\cubed}.");
        prm.declare_entry ("Reference isothermal compressibilities", "4e-12",
                           Patterns::Anything(),
                           "List of isothermal compressibilities for background and compositional fields (N), "
                           "for a total of N+1 values for models with no phase transitions (or models where the "
                           "value does not change across any of the phase transitions). For models with phase "
                           "transitions, the list needs to contain each field name, including the background, for "
                           "a total of N+1 names, and for each of these names, specify the value for each phase "
                           "(giving P_c+1 values for each field, with P_c being the number of phase transitions "
                           "for field c). Therefore, the total number of values given is N+P+1, with P = sum(P_c) "
                           "the total number of phase transitions, summed over all phases. The format is "
                           "background: value1|value2|...|valueP_1+1, field1:value1|...|valueP_2+1, ..., fieldN: value1|...|valueP_N+1. "
                           "If only one value is given, then all fields/phases use the same value. "
                           "Units: \\si{\\per\\pascal}.");
        prm.declare_entry ("Isothermal bulk modulus pressure derivatives", "4.",
                           Patterns::Anything(),
                           "List of isothermal bulk modulus pressure derivatives for background and compositional fields (N), "
                           "for a total of N+1 values for models with no phase transitions (or models where the "
                           "value does not change across any of the phase transitions). For models with phase "
                           "transitions, the list needs to contain each field name, including the background, for "
                           "a total of N+1 names, and for each of these names, specify the value for each phase "
                           "(giving P_c+1 values for each field, with P_c being the number of phase transitions "
                           "for field c). Therefore, the total number of values given is N+P+1, with P = sum(P_c) "
                           "the total number of phase transitions, summed over all phases. The format is "
                           "background: value1|value2|...|valueP_1+1, field1:value1|...|valueP_2+1, ..., fieldN: value1|...|valueP_N+1. "
                           "If only one value is given, then all fields/phases use the same value. "
                           "Units: [].");
        prm.declare_entry ("Reference thermal expansivities", "4.e-5",
                           Patterns::Anything(),
                           "List of reference thermal expansivities for background and compositional fields (N), "
                           "for a total of N+1 values for models with no phase transitions (or models where the "
                           "value does not change across any of the phase transitions). For models with phase "
                           "transitions, the list needs to contain each field name, including the background, for "
                           "a total of N+1 names, and for each of these names, specify the value for each phase "
                           "(giving P_c+1 values for each field, with P_c being the number of phase transitions "
                           "for field c). Therefore, the total number of values given is N+P+1, with P = sum(P_c) "
                           "the total number of phase transitions, summed over all phases. The format is "
                           "background: value1|value2|...|valueP_1+1, field1:value1|...|valueP_2+1, ..., fieldN: value1|...|valueP_N+1. "
                           "If only one value is given, then all fields/phases use the same value. "
                           "Units: \\si{\\per\\kelvin}.");
        prm.declare_entry ("Isochoric specific heats", "1250.",
                           Patterns::Anything(),
                           "List of isochoric specific heats for background and compositional fields (N), "
                           "for a total of N+1 values for models with no phase transitions (or models where the "
                           "value does not change across any of the phase transitions). For models with phase "
                           "transitions, the list needs to contain each field name, including the background, for "
                           "a total of N+1 names, and for each of these names, specify the value for each phase "
                           "(giving P_c+1 values for each field, with P_c being the number of phase transitions "
                           "for field c). Therefore, the total number of values given is N+P+1, with P = sum(P_c) "
                           "the total number of phase transitions, summed over all phases. The format is "
                           "background: value1|value2|...|valueP_1+1, field1:value1|...|valueP_2+1, ..., fieldN: value1|...|valueP_N+1. "
                           "If only one value is given, then all fields/phases use the same value. "
                           "Units: \\si{\\joule\\per\\kelvin\\per\\kilogram}.");
        prm.declare_entry ("Enable phase transitions", "false",
                           Patterns::Bool (),
                           "Whether to enable the use of phase transitions, which break the thermodynamic "
                           "consistency of the equation of state for properties (heat capacity, thermal "
                           "expansivity and compressibility) that are affected by the P-T-X dependence "
                           "of the phase transition.");
      }



      template <int dim>
      void
      MulticomponentCompressible<dim>::parse_parameters (ParameterHandler &prm,
                                                         const std::unique_ptr<std::vector<unsigned int>> &expected_n_phases_per_composition)
      {
        std::vector<std::string> compositional_field_names = this->introspection().get_composition_names();

        // Establish that a background field is required here
        compositional_field_names.insert(compositional_field_names.begin(),"background");

        std::vector<std::string> chemical_field_names = this->introspection().chemical_composition_field_names();
        chemical_field_names.insert(chemical_field_names.begin(),"background");

        Utilities::MapParsing::Options options(compositional_field_names, "");
        options.list_of_allowed_keys = compositional_field_names;
        options.allow_multiple_values_per_key = true;
        if (expected_n_phases_per_composition)
          {
            options.n_values_per_key = *expected_n_phases_per_composition;

            // check_values_per_key is required to be true to duplicate single values
            // if they are to be used for all phases associated with a given key.
            options.check_values_per_key = true;
          }
        else
          {
            // If the material model does not tell us how many phases per composition to expect,
            // at least check that the parameters parsed below have the same number of values.
            options.store_values_per_key = true;
          }

        // Parse multicomponent properties
        options.property_name = "Reference temperatures";
        reference_temperatures = Utilities::MapParsing::parse_map_to_double_array (prm.get(options.property_name), options);

        // Now that we know we have stored the number of phases per composition, we can check them for subsequent properties.
        options.store_values_per_key = false;
        options.check_values_per_key = true;

        options.property_name = "Reference densities";
        reference_densities = Utilities::MapParsing::parse_map_to_double_array (prm.get(options.property_name), options);

        options.property_name = "Reference isothermal compressibilities";
        reference_isothermal_compressibilities = Utilities::MapParsing::parse_map_to_double_array (prm.get(options.property_name), options);

        options.property_name = "Isothermal bulk modulus pressure derivatives";
        isothermal_bulk_modulus_pressure_derivatives = Utilities::MapParsing::parse_map_to_double_array (prm.get(options.property_name), options);

        options.property_name = "Reference thermal expansivities";
        reference_thermal_expansivities = Utilities::MapParsing::parse_map_to_double_array (prm.get(options.property_name), options);

        options.property_name = "Isochoric specific heats";
        isochoric_specific_heats = Utilities::MapParsing::parse_map_to_double_array (prm.get(options.property_name), options);

        enable_phase_transitions = prm.get_bool ("Enable phase transitions");

        // Check that phase transitions have explicitly been enabled if the number of phases
        // for any compositional field is greater than 0
        if (!enable_phase_transitions)
          {
            for (unsigned int n=0; n < (*expected_n_phases_per_composition).size(); ++n)
              {
                AssertThrow((*expected_n_phases_per_composition)[n] == 1, ExcMessage(
                              "The number of expected phases per composition for the compositional "
                              "field named (" + compositional_field_names[n] + ") is equal to "
                              "(" + Utilities::to_string((*expected_n_phases_per_composition)[n]) + "). "
                              "Use of phase transitions with the multicomponent compressible equation of "
                              "state model is currently implemented in a thermodynamically inconsistent manner "
                              "and must be explicitly enabled with 'set Enable phase transitions = true'."));
              }
          }
      }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    namespace EquationOfState
    {
#define INSTANTIATE(dim) \
  template class MulticomponentCompressible<dim>;

      ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
    }
  }
}
