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
        std::vector<std::string> compositional_field_names = this->introspection().get_composition_names();
        // Establish that a background field is required here
        compositional_field_names.insert(compositional_field_names.begin(),"background");
        Utilities::MapParsing::Options options(compositional_field_names, "");

        // Parse multicomponent properties
        options.property_name = "Reference temperatures";
        reference_temperatures = Utilities::MapParsing::parse_map_to_double_array (prm.get(options.property_name), options);

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
