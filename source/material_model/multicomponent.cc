/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

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


#include <aspect/material_model/multicomponent.h>
#include <aspect/utilities.h>

#include <numeric>


namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    void
    Multicomponent<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      EquationOfStateOutputs<dim> eos_outputs (this->introspection().n_chemical_composition_fields()+1);

      for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
        {
          // The (incompressible) Boussinesq approximation treats the
          // buoyancy term as Delta rho[i] * C[i], which implies that
          // compositional fields are given as volume fractions.
          const std::vector<double> volume_fractions = MaterialUtilities::compute_only_composition_fractions(in.composition[i],
                                                       this->introspection().chemical_composition_field_indices());

          equation_of_state.evaluate(in, i, eos_outputs);
          out.viscosities[i] = MaterialUtilities::average_value(volume_fractions, viscosities, viscosity_averaging);
          out.specific_heat[i] = MaterialUtilities::average_value (volume_fractions, eos_outputs.specific_heat_capacities, MaterialUtilities::arithmetic);
          out.thermal_expansion_coefficients[i] = MaterialUtilities::average_value (volume_fractions, eos_outputs.thermal_expansion_coefficients, MaterialUtilities::arithmetic);
          out.compressibilities[i] = MaterialUtilities::average_value (volume_fractions, eos_outputs.compressibilities, MaterialUtilities::arithmetic);
          out.entropy_derivative_pressure[i] = MaterialUtilities::average_value (volume_fractions, eos_outputs.entropy_derivative_pressure, MaterialUtilities::arithmetic);
          out.entropy_derivative_temperature[i] = MaterialUtilities::average_value (volume_fractions, eos_outputs.entropy_derivative_temperature, MaterialUtilities::arithmetic);

          // Arithmetic averaging of thermal conductivities
          // This may not be strictly the most reasonable thing, but for most Earth materials we hope
          // that they do not vary so much that it is a big problem.
          out.thermal_conductivities[i] = MaterialUtilities::average_value (volume_fractions, thermal_conductivities, MaterialUtilities::arithmetic);

          // not strictly correct if thermal expansivities are different, since we are interpreting
          // these compositions as volume fractions, but the error introduced should not be too bad.
          out.densities[i] = MaterialUtilities::average_value (volume_fractions, eos_outputs.densities, MaterialUtilities::arithmetic);

          for (unsigned int c=0; c<in.composition[i].size(); ++c)
            out.reaction_terms[i][c] = 0.0;

        }
    }

    template <int dim>
    bool
    Multicomponent<dim>::
    is_compressible () const
    {
      return equation_of_state.is_compressible();
    }

    template <int dim>
    void
    Multicomponent<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Multicomponent");
        {
          EquationOfState::MulticomponentIncompressible<dim>::declare_parameters (prm, 4.e-5);

          prm.declare_entry ("Reference temperature", "293.",
                             Patterns::Double (0.),
                             "The reference temperature $T_0$. Units: $\\text{K}$.");
          prm.declare_entry ("Viscosities", "1.e21",
                             Patterns::Anything(),
                             "List of viscosities for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one value is given, then all use the same value. Units: $\\text{Pa}\\text{s}$.");
          prm.declare_entry ("Thermal conductivities", "4.7",
                             Patterns::Anything(),
                             "List of thermal conductivities for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one value is given, then all use the same value. "
                             "Units: $\\frac{\\text{W}{\\text{m}\\text{K}}$.");
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
    Multicomponent<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Multicomponent");
        {
          // Make options file for parsing maps to double arrays
          std::vector<std::string> chemical_field_names = this->introspection().chemical_composition_field_names();
          chemical_field_names.insert(chemical_field_names.begin(),"background");

          std::vector<std::string> compositional_field_names = this->introspection().get_composition_names();
          compositional_field_names.insert(compositional_field_names.begin(),"background");

          Utilities::MapParsing::Options options(chemical_field_names, "");
          options.list_of_allowed_keys = compositional_field_names;
          options.allow_multiple_values_per_key = true;
          options.check_values_per_key = true;

          // We do not have phase transitions. Make sure there is only one entry for each field.
          const std::vector<unsigned int> n_phases_for_each_chemical_composition (chemical_field_names.size(), 1);
          options.n_values_per_key = n_phases_for_each_chemical_composition;

          equation_of_state.initialize_simulator (this->get_simulator());
          equation_of_state.parse_parameters (prm,
                                              std::make_unique<std::vector<unsigned int>>(n_phases_for_each_chemical_composition));

          reference_T = prm.get_double ("Reference temperature");

          viscosity_averaging = MaterialUtilities::parse_compositional_averaging_operation ("Viscosity averaging scheme",
                                prm);

          // Parse multicomponent properties
          options.property_name = "Viscosities";
          viscosities = Utilities::MapParsing::parse_map_to_double_array (prm.get("Viscosities"), options);
          options.property_name = "Thermal conductivities";
          thermal_conductivities = Utilities::MapParsing::parse_map_to_double_array (prm.get("Thermal conductivities"), options);
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // Declare dependencies on solution variables
      this->model_dependence.viscosity = NonlinearDependence::compositional_fields;
      this->model_dependence.density = NonlinearDependence::temperature | NonlinearDependence::compositional_fields;
      this->model_dependence.compressibility = NonlinearDependence::none;
      this->model_dependence.specific_heat = NonlinearDependence::compositional_fields;
      this->model_dependence.thermal_conductivity = NonlinearDependence::compositional_fields;
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(Multicomponent,
                                   "multicomponent",
                                   "This incompressible model is for use with an arbitrary number of"
                                   " compositional fields, where each field represents a rock type which"
                                   " can have completely different properties from the others."
                                   " However, each rock type itself has constant material properties.  The value of the "
                                   " compositional field is interpreted as a volume fraction. If the sum of the fields is"
                                   " greater than one, they are renormalized.  If it is less than one, material properties "
                                   " for ``background mantle'' make up the rest. When more than one field is present, the"
                                   " material properties are averaged arithmetically.  An exception is the viscosity,"
                                   " where the averaging should make more of a difference.  For this, the user selects"
                                   " between arithmetic, harmonic, geometric, or maximum composition averaging.")
  }
}
