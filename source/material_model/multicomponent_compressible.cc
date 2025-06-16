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


#include <aspect/material_model/multicomponent_compressible.h>
#include <aspect/utilities.h>
#include <deal.II/base/parameter_handler.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/gravity_model/interface.h>

#include <numeric>


namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    void
    MulticomponentCompressible<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      EquationOfStateOutputs<dim> eos_outputs (this->n_compositional_fields()+1);
      EquationOfStateOutputs<dim> eos_outputs_all_phases (n_phases);

      // Store the phase function value for each phase and composition
      // While the number of phases is fixed, the value of the phase function is updated for every point
      std::vector<double> phase_function_values(phase_function.n_phase_transitions(), 0.0);

      unsigned int density_field_index = numbers::invalid_unsigned_int;

      if (this->introspection().composition_type_exists(Parameters<dim>::CompositionalFieldDescription::density))
        density_field_index = this->introspection().find_composition_type(Parameters<dim>::CompositionalFieldDescription::density);

      for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
        {
          equation_of_state.evaluate(in, i, eos_outputs_all_phases);

          const double gravity_norm = this->get_gravity_model().gravity_vector(in.position[i]).norm();
          const double reference_density = (this->get_adiabatic_conditions().is_initialized())
                                           ?
                                           this->get_adiabatic_conditions().density(in.position[i])
                                           :
                                           eos_outputs_all_phases.densities[0];

          // The phase index is set to invalid_unsigned_int, because it is only used internally
          // in phase_average_equation_of_state_outputs to loop over all existing phases
          MaterialUtilities::PhaseFunctionInputs<dim> phase_inputs(in.temperature[i],
                                                                   in.pressure[i],
                                                                   this->get_geometry_model().depth(in.position[i]),
                                                                   gravity_norm*reference_density,
                                                                   numbers::invalid_unsigned_int);

          // Compute value of phase functions
          for (unsigned int j=0; j < phase_function.n_phase_transitions(); ++j)
            {
              phase_inputs.phase_transition_index = j;
              phase_function_values[j] = phase_function.compute_value(phase_inputs);
            }

          // Average by value of the phase function to get value of compositions
          phase_average_equation_of_state_outputs(eos_outputs_all_phases,
                                                  phase_function_values,
                                                  n_phase_transitions_for_each_chemical_composition,
                                                  eos_outputs);

          // Calculate volume fractions from mass fractions
          const std::vector<double> mass_fractions = MaterialUtilities::compute_only_composition_fractions(in.composition[i],
                                                     this->introspection().chemical_composition_field_indices());
          const std::vector<double> volume_fractions = MaterialUtilities::compute_volumes_from_masses(mass_fractions,
                                                       eos_outputs.densities,
                                                       true);

          // Specific quantities are mass averaged
          out.specific_heat[i] = MaterialUtilities::average_value (mass_fractions, eos_outputs.specific_heat_capacities, MaterialUtilities::arithmetic);
          out.entropy_derivative_pressure[i] = MaterialUtilities::average_value (mass_fractions, eos_outputs.entropy_derivative_pressure, MaterialUtilities::arithmetic);
          out.entropy_derivative_temperature[i] = MaterialUtilities::average_value (mass_fractions, eos_outputs.entropy_derivative_temperature, MaterialUtilities::arithmetic);

          // Density, thermal expansivity and compressibility are all volume averaged
          out.densities[i] = MaterialUtilities::average_value (volume_fractions, eos_outputs.densities, MaterialUtilities::arithmetic);
          out.thermal_expansion_coefficients[i] = MaterialUtilities::average_value (volume_fractions, eos_outputs.thermal_expansion_coefficients, MaterialUtilities::arithmetic);
          out.compressibilities[i] = MaterialUtilities::average_value (volume_fractions, eos_outputs.compressibilities, MaterialUtilities::arithmetic);

          // Arithmetic volume fraction averaging of thermal conductivities
          // This may not be the most reasonable thing, but for most Earth materials we hope
          // that they do not vary so much that it is a big problem.
          std::vector<double> phase_averaged_thermal_conductivities(volume_fractions.size());
          for (unsigned int c=0; c<volume_fractions.size(); ++c)
            phase_averaged_thermal_conductivities[c] = MaterialUtilities::phase_average_value(phase_function_values,
                                                       n_phase_transitions_for_each_chemical_composition,
                                                       thermal_conductivities,
                                                       c);
          out.thermal_conductivities[i] = MaterialUtilities::average_value (volume_fractions, phase_averaged_thermal_conductivities, MaterialUtilities::arithmetic);


          // User-defined volume fraction averaging of viscosities
          if (in.requests_property(MaterialProperties::viscosity) || in.requests_property(MaterialProperties::additional_outputs) ||
              (this->get_parameters().enable_elasticity && in.requests_property(MaterialProperties::reaction_rates) ))
            {
              std::vector<double> phase_averaged_viscosities(volume_fractions.size());
              for (unsigned int c=0; c<volume_fractions.size(); ++c)
                phase_averaged_viscosities[c] = MaterialUtilities::phase_average_value(phase_function_values,
                                                                                       n_phase_transitions_for_each_chemical_composition,
                                                                                       viscosities,
                                                                                       c);

              out.viscosities[i] = MaterialUtilities::average_value (volume_fractions, phase_averaged_viscosities, viscosity_averaging);
            }
          else
            {
              out.viscosities[i] = numbers::signaling_nan<double>();
            }

          for (unsigned int c=0; c<in.composition[i].size(); ++c)
            out.reaction_terms[i][c] = 0.0;

          // set up variable to interpolate prescribed field outputs onto compositional fields
          if (const std::shared_ptr<PrescribedFieldOutputs<dim>> prescribed_field_out
              = out.template get_additional_output_object<PrescribedFieldOutputs<dim>>())
            {
              prescribed_field_out->prescribed_field_outputs[i][density_field_index] = out.densities[i];
            }
        }
    }



    template <int dim>
    bool
    MulticomponentCompressible<dim>::
    is_compressible () const
    {
      return equation_of_state.is_compressible();
    }



    template <int dim>
    void
    MulticomponentCompressible<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Multicomponent compressible");
        {
          MaterialUtilities::PhaseFunction<dim>::declare_parameters(prm);

          EquationOfState::MulticomponentCompressible<dim>::declare_parameters (prm);

          prm.declare_entry ("Viscosities", "1.e21",
                             Patterns::Anything(),
                             "List of viscosities for background and compositional fields (N), "
                             "for a total of N+1 values for models with no phase transitions (or models where the "
                             "value does not change across any of the phase transitions). For models with phase "
                             "transitions, the list needs to contain each field name, including the background, for "
                             "a total of N+1 names, and for each of these names, specify the value for each phase "
                             "(giving P_c+1 values for each field, with P_c being the number of phase transitions "
                             "for field c). Therefore, the total number of values given is N+P+1, with P = sum(P_c) "
                             "the total number of phase transitions, summed over all phases. The format is "
                             "background: value1|value2|...|valueP_1+1, field1:value1|...|valueP_2+1, ..., fieldN: value1|...|valueP_N+1. "
                             "If only one value is given, then all fields/phases use the same value. "
                             "Units: \\si{\\pascal\\second}.");
          prm.declare_entry ("Thermal conductivities", "4.7",
                             Patterns::Anything(),
                             "List of thermal conductivities for background and compositional fields (N), "
                             "for a total of N+1 values for models with no phase transitions (or models where the "
                             "value does not change across any of the phase transitions). For models with phase "
                             "transitions, the list needs to contain each field name, including the background, for "
                             "a total of N+1 names, and for each of these names, specify the value for each phase "
                             "(giving P_c+1 values for each field, with P_c being the number of phase transitions "
                             "for field c). Therefore, the total number of values given is N+P+1, with P = sum(P_c) "
                             "the total number of phase transitions, summed over all phases. The format is "
                             "background: value1|value2|...|valueP_1+1, field1:value1|...|valueP_2+1, ..., fieldN: value1|...|valueP_N+1. "
                             "If only one value is given, then all fields/phases use the same value. "
                             "Units: \\si{\\watt\\per\\meter\\per\\kelvin}.");
          prm.declare_entry ("Viscosity averaging scheme", "harmonic",
                             Patterns::Selection("arithmetic|harmonic|geometric|maximum composition"),
                             "When more than one compositional field is present at a point "
                             "with different viscosities, we need to come up with an average "
                             "viscosity at that point.  Select a weighted harmonic, arithmetic, "
                             "geometric, or maximum composition.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    MulticomponentCompressible<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Multicomponent compressible");
        {

          // Phase transition parameters
          phase_function.initialize_simulator (this->get_simulator());
          phase_function.parse_parameters (prm);

          const std::vector<unsigned int> n_phases_for_each_chemical_composition = phase_function.n_phases_for_each_chemical_composition();
          n_phase_transitions_for_each_chemical_composition = phase_function.n_phase_transitions_for_each_chemical_composition();
          n_phases = phase_function.n_phases_over_all_chemical_compositions();

          equation_of_state.initialize_simulator (this->get_simulator());
          equation_of_state.parse_parameters (prm,
                                              std::make_unique<std::vector<unsigned int>>(n_phases_for_each_chemical_composition));

          viscosity_averaging = MaterialUtilities::parse_compositional_averaging_operation ("Viscosity averaging scheme",
                                prm);

          // Make options file for parsing maps to double arrays
          std::vector<std::string> chemical_field_names = this->introspection().chemical_composition_field_names();
          chemical_field_names.insert(chemical_field_names.begin(),"background");

          std::vector<std::string> compositional_field_names = this->introspection().get_composition_names();
          compositional_field_names.insert(compositional_field_names.begin(),"background");

          Utilities::MapParsing::Options options(chemical_field_names, "");
          options.list_of_allowed_keys = compositional_field_names;
          options.allow_multiple_values_per_key = true;
          options.check_values_per_key = true;
          options.n_values_per_key = n_phases_for_each_chemical_composition;
          options.check_values_per_key = (options.n_values_per_key.size() != 0);
          options.store_values_per_key = (options.n_values_per_key.size() == 0);

          options.property_name = "Viscosities";
          viscosities = Utilities::MapParsing::parse_map_to_double_array (prm.get(options.property_name), options);

          options.property_name = "Thermal conductivities";
          thermal_conductivities = Utilities::MapParsing::parse_map_to_double_array (prm.get(options.property_name), options);
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // Declare dependencies on solution variables
      this->model_dependence.thermal_conductivity = NonlinearDependence::compositional_fields;
      this->model_dependence.viscosity = NonlinearDependence::compositional_fields;
      this->model_dependence.density |= NonlinearDependence::temperature
                                        | NonlinearDependence::pressure
                                        | NonlinearDependence::compositional_fields;
      this->model_dependence.compressibility = NonlinearDependence::temperature
                                               | NonlinearDependence::pressure
                                               | NonlinearDependence::compositional_fields;
      this->model_dependence.specific_heat = NonlinearDependence::temperature
                                             | NonlinearDependence::pressure
                                             | NonlinearDependence::compositional_fields;
    }



    template <int dim>
    void
    MulticomponentCompressible<dim>::create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      if (this->introspection().composition_type_exists(Parameters<dim>::CompositionalFieldDescription::density))
        {
          if (out.template has_additional_output_object<PrescribedFieldOutputs<dim>>() == false)
            {
              const unsigned int n_points = out.n_evaluation_points();
              out.additional_outputs.push_back(
                std::make_unique<MaterialModel::PrescribedFieldOutputs<dim>>
                (n_points, this->n_compositional_fields()));
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
    ASPECT_REGISTER_MATERIAL_MODEL(MulticomponentCompressible,
                                   "multicomponent compressible",
                                   "This model is for use with an arbitrary number of compositional fields, where each field"
                                   " represents a rock type which can have completely different properties from the others."
                                   " Each rock type is described by a self-consistent equation of state.  The value of the "
                                   " compositional field is interpreted as a mass fraction. If the sum of the fields is"
                                   " greater than one, they are renormalized.  If it is less than one, material properties "
                                   " for ``background mantle'' make up the rest. When more than one field is present, the "
                                   "material properties are averaged arithmetically by mass fraction (for specific heat), "
                                   "or volume fraction (for density, thermal expansivity and compressibility). "
                                   "The thermal conductivity is also arithmetically averaged by volume fraction. "
                                   "Finally, the viscosity is averaged by volume fraction, but the user can choose "
                                   "between arithmetic, harmonic, geometric or maximum composition averaging.")
  }
}
