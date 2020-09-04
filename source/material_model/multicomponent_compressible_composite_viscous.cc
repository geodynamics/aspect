/*
  Copyright (C) 2011 - 2020 by the authors of the ASPECT code.

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


#include <aspect/material_model/multicomponent_compressible_composite_viscous.h>
#include <aspect/utilities.h>

#include <numeric>


namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    void
    MulticomponentCompressibleCompositeViscous<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      EquationOfStateOutputs<dim> eos_outputs (this->n_compositional_fields()+1);

      // Store value of phase function for each phase and composition
      // While the number of phases is fixed, the value of the phase function is updated for every point
      std::vector<double> phase_function_values(phase_function.n_phase_transitions(), 0.0);

      for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
        {
          equation_of_state.evaluate(in, i, eos_outputs);

          const double gravity_norm = this->get_gravity_model().gravity_vector(in.position[i]).norm();
          const double reference_density = (this->get_adiabatic_conditions().is_initialized())
                                           ?
                                           this->get_adiabatic_conditions().density(in.position[i])
                                           :
                                           eos_outputs.densities[0];

          // The phase index is set to invalid_unsigned_int, because it is only used internally
          // in phase_average_equation_of_state_outputs to loop over all existing phases
          MaterialUtilities::PhaseFunctionInputs<dim> phase_inputs(in.temperature[i],
                                                                   in.pressure[i],
                                                                   this->get_geometry_model().depth(in.position[i]),
                                                                   gravity_norm*reference_density,
                                                                   numbers::invalid_unsigned_int);

          // Compute value of phase functions
          for (unsigned int j=0; j < phase_function.n_phase_transitions(); j++)
            {
              phase_inputs.phase_index = j;
              phase_function_values[j] = phase_function.compute_value(phase_inputs);
            }

          // Calculate volume fractions from mass fractions
          const std::vector<double> mass_fractions = MaterialUtilities::compute_field_fractions(in.composition[i]);
          const unsigned int n_fields = mass_fractions.size();
          std::vector<double> volume_fractions(n_fields);
          std::vector<double> viscosities(n_fields);
          double sum_volume_fractions = 0.0;
          for (unsigned int j=0; j < n_fields; ++j)
            {
              volume_fractions[j] = mass_fractions[j] / eos_outputs.densities[j];
              sum_volume_fractions += volume_fractions[j];

              viscosities[j] = rheology.compute_viscosity(in.pressure[i], in.temperature[i],
                                                          j, in.strain_rate[i],
                                                          phase_function_values,
                                                          phase_function.n_phase_transitions_for_each_composition());
            }

          for (unsigned int j=0; j < n_fields; ++j)
            volume_fractions[j] /= sum_volume_fractions;

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
          out.thermal_conductivities[i] = MaterialUtilities::average_value (volume_fractions, thermal_conductivities, MaterialUtilities::arithmetic);

          // User-defined volume fraction averaging of viscosities
          out.viscosities[i] = MaterialUtilities::average_value (volume_fractions, viscosities, viscosity_averaging);

          for (unsigned int c=0; c<in.composition[i].size(); ++c)
            out.reaction_terms[i][c] = 0.0;

        }
    }

    template <int dim>
    double
    MulticomponentCompressibleCompositeViscous<dim>::
    reference_viscosity () const
    {
      return viscosities[0]; // background
    }

    template <int dim>
    bool
    MulticomponentCompressibleCompositeViscous<dim>::
    is_compressible () const
    {
      return equation_of_state.is_compressible();
    }

    template <int dim>
    void
    MulticomponentCompressibleCompositeViscous<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Multicomponent compressible composite viscous");
        {
          EquationOfState::MulticomponentCompressible<dim>::declare_parameters (prm);
          Rheology::CompositeViscousCreep<dim>::declare_parameters (prm);

          prm.declare_entry ("Thermal conductivities", "4.7",
                             Patterns::Anything(),
                             "List of thermal conductivities for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one value is given, then all use the same value. "
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
    MulticomponentCompressibleCompositeViscous<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Multicomponent compressible composite viscous");
        {
          // Phase transition parameters
          phase_function.initialize_simulator (this->get_simulator());
          phase_function.parse_parameters (prm);

          std::vector<unsigned int> n_phase_transitions_for_each_composition
          (phase_function.n_phase_transitions_for_each_composition());

          // We require one more entry for density, etc as there are phase transitions
          // (for the low-pressure phase before any transition).
          for (unsigned int i=0; i<n_phase_transitions_for_each_composition.size(); ++i)
            n_phase_transitions_for_each_composition[i] += 1;

          // Equation of state parameters
          equation_of_state.initialize_simulator (this->get_simulator());
          equation_of_state.parse_parameters (prm);

          // Rheology parameters
          rheology.initialize_simulator (this->get_simulator());
          rheology.parse_parameters (prm, std::make_shared<std::vector<unsigned int>>(n_phase_transitions_for_each_composition));

          viscosity_averaging = MaterialUtilities::parse_compositional_averaging_operation ("Viscosity averaging scheme",
                                prm);

          // Establish that a background field is required here
          const bool has_background_field = true;

          // Retrieve the list of composition names
          const std::vector<std::string> list_of_composition_names = this->introspection().get_composition_names();

          thermal_conductivities = Utilities::parse_map_to_double_array (prm.get("Thermal conductivities"),
                                                                         list_of_composition_names,
                                                                         has_background_field,
                                                                         "Thermal conductivities");
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
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(MulticomponentCompressibleCompositeViscous,
                                   "multicomponent compressible composite viscous",
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
