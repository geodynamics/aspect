/*
  Copyright (C) 2011 - 2022 by the authors of the ASPECT code.

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

      for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
        {
          equation_of_state.evaluate(in, i, eos_outputs);

          // Calculate volume fractions from mass fractions
          const std::vector<double> mass_fractions = MaterialUtilities::compute_composition_fractions(in.composition[i]);
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
          out.thermal_conductivities[i] = MaterialUtilities::average_value (volume_fractions, thermal_conductivities, MaterialUtilities::arithmetic);

          // User-defined volume fraction averaging of viscosities
          out.viscosities[i] = MaterialUtilities::average_value (volume_fractions, viscosities, viscosity_averaging);

          for (unsigned int c=0; c<in.composition[i].size(); ++c)
            out.reaction_terms[i][c] = 0.0;

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
          EquationOfState::MulticomponentCompressible<dim>::declare_parameters (prm);

          prm.declare_entry ("Viscosities", "1.e21",
                             Patterns::Anything(),
                             "List of viscosities for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one value is given, then all use the same value. Units: \\si{\\pascal\\second}.");
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
    MulticomponentCompressible<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Multicomponent compressible");
        {
          equation_of_state.initialize_simulator (this->get_simulator());
          equation_of_state.parse_parameters (prm);

          viscosity_averaging = MaterialUtilities::parse_compositional_averaging_operation ("Viscosity averaging scheme",
                                prm);

          std::vector<std::string> compositional_field_names = this->introspection().get_composition_names();
          // Establish that a background field is required here
          compositional_field_names.insert(compositional_field_names.begin(),"background");
          Utilities::MapParsing::Options options(compositional_field_names, "");

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
