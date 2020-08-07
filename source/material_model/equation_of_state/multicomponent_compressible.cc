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


#include <aspect/material_model/equation_of_state/multicomponent_compressible.h>
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
               MaterialModel::EquationOfStateOutputs<dim> &out) const
      {

        const double pressure = in.pressure[q];
        const double temperature = std::max(in.temperature[q], 1.); // temperature can't be zero for correct evaluation

        for (unsigned int c=0; c < out.densities.size(); ++c)
          {
            const double ak = reference_thermal_expansivities[c]/reference_isothermal_compressibilities[c];
            const double f = (1. + (pressure - ak*(temperature - reference_temperatures[c])) *
                              isothermal_bulk_modulus_pressure_derivatives[c] *
                              reference_isothermal_compressibilities[c]);

            out.densities[c] = reference_densities[c]*std::pow(f, 1./isothermal_bulk_modulus_pressure_derivatives[c]);
            out.thermal_expansion_coefficients[c] = reference_thermal_expansivities[c] / f;
            out.specific_heat_capacities[c] = (isochoric_specific_heats[c] +
                                               (temperature*reference_thermal_expansivities[c] *
                                                ak * std::pow(f, -1.-(1./isothermal_bulk_modulus_pressure_derivatives[c]))
                                                / reference_densities[c]));
            out.compressibilities[c] = reference_isothermal_compressibilities[c]/f;
            out.entropy_derivative_pressure[c] = 0.;
            out.entropy_derivative_temperature[c] = 0.;
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
                           Patterns::List(Patterns::Double (0.)),
                           "List of reference temperatures $T_0$ for background mantle and compositional fields,"
                           "for a total of N+1 values, where N is the number of compositional fields."
                           "If only one value is given, then all use the same value. Units: \\si{\\kelvin}.");
        prm.declare_entry ("Reference densities", "3300.",
                           Patterns::List(Patterns::Double (0.)),
                           "List of densities for background mantle and compositional fields,"
                           "for a total of N+1 values, where N is the number of compositional fields."
                           "If only one value is given, then all use the same value. "
                           "Units: \\si{\\kilogram\\per\\meter\\cubed}.");
        prm.declare_entry ("Reference isothermal compressibilities", "4e-12",
                           Patterns::List(Patterns::Double (0.)),
                           "List of isothermal compressibilities for background mantle and compositional fields,"
                           "for a total of N+1 values, where N is the number of compositional fields."
                           "If only one value is given, then all use the same value. "
                           "Units: \\si{\\per\\pascal}.");
        prm.declare_entry ("Isothermal bulk modulus pressure derivatives", "4.",
                           Patterns::List(Patterns::Double (0.)),
                           "List of isothermal pressure derivatives of the bulk moduli for background mantle and compositional fields,"
                           "for a total of N+1 values, where N is the number of compositional fields."
                           "If only one value is given, then all use the same value. "
                           "Units: [].");
        prm.declare_entry ("Reference thermal expansivities", "4.e-5",
                           Patterns::List(Patterns::Double (0.)),
                           "List of thermal expansivities for background mantle and compositional fields,"
                           "for a total of N+1 values, where N is the number of compositional fields."
                           "If only one value is given, then all use the same value. Units: \\si{\\per\\kelvin}.");
        prm.declare_entry ("Isochoric specific heats", "1250.",
                           Patterns::List(Patterns::Double (0.)),
                           "List of isochoric specific heats $C_v$ for background mantle and compositional fields,"
                           "for a total of N+1 values, where N is the number of compositional fields."
                           "If only one value is given, then all use the same value. "
                           "Units: \\si{\\joule\\per\\kelvin\\per\\kilogram}.");
      }



      template <int dim>
      void
      MulticomponentCompressible<dim>::parse_parameters (ParameterHandler &prm)
      {
        // Establish that a background field is required here
        const bool has_background_field = true;

        // Retrieve the list of composition names
        const std::vector<std::string> list_of_composition_names = this->introspection().get_composition_names();

        // Parse multicomponent properties
        reference_temperatures = Utilities::parse_map_to_double_array (prm.get("Reference temperatures"),
                                                                       list_of_composition_names,
                                                                       has_background_field,
                                                                       "Reference temperatures");

        reference_densities = Utilities::parse_map_to_double_array (prm.get("Reference densities"),
                                                                    list_of_composition_names,
                                                                    has_background_field,
                                                                    "Reference densities");

        reference_isothermal_compressibilities = Utilities::parse_map_to_double_array (prm.get("Reference isothermal compressibilities"),
                                                                                       list_of_composition_names,
                                                                                       has_background_field,
                                                                                       "Reference isothermal compressibilities");

        isothermal_bulk_modulus_pressure_derivatives = Utilities::parse_map_to_double_array (prm.get("Isothermal bulk modulus pressure derivatives"),
                                                       list_of_composition_names,
                                                       has_background_field,
                                                       "Isothermal bulk modulus pressure derivatives");

        reference_thermal_expansivities = Utilities::parse_map_to_double_array (prm.get("Reference thermal expansivities"),
                                                                                list_of_composition_names,
                                                                                has_background_field,
                                                                                "Reference thermal expansivities");

        isochoric_specific_heats = Utilities::parse_map_to_double_array (prm.get("Isochoric specific heats"),
                                                                         list_of_composition_names,
                                                                         has_background_field,
                                                                         "Isochoric specific heats");
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
