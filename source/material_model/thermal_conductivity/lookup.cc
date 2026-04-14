/*
  Copyright (C) 2025 - by the authors of the ASPECT code.

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

#include <aspect/material_model/thermal_conductivity/lookup.h>

namespace aspect
{
  namespace MaterialModel
  {
    namespace ThermalConductivity
    {
      template <int dim>
      void
      Lookup<dim>::evaluate (const MaterialModel::MaterialModelInputs<dim> &in,
                               MaterialModel::MaterialModelOutputs<dim> &out) const
      {
        // Fetch the names of the chemical composition fields
        const std::vector<std::string> chemical_field_names = this->introspection().chemical_composition_field_names();

        const unsigned int n_evaluation_points = in.n_evaluation_points();
        for (unsigned int i=0; i<n_evaluation_points; ++i)
          {
            const std::vector<double> volume_fractions = MaterialUtilities::compute_only_composition_fractions(in.composition[i],
                                                         this->introspection().chemical_composition_field_indices());

            // Reset the thermal diffusivity value before entering a new evaluation point
            double thermal_diffusivity = 0.0;
            const double &temperature = in.temperature[i];

            for (unsigned int j = 0; j < volume_fractions.size(); ++j)
              {
                // By default do not scale the diffusivity values
                double scale_factor = 1.0;
                
                if (j > 0)
                {
                  // The volume_fractions vector has N+1 values as it has background included, so take that into account 
                  const std::string &field_name = chemical_field_names[j - 1];
                  // Check if this composition has a specific name and apply scaling
                  if (field_name == "1_pyrolite_M" || field_name == "2_harzb_LM" || field_name == "7_wcont_LM")
                    scale_factor = 7673.33 * std::pow(temperature, -1.366) + 2.967e-4 * temperature;
                  else if (field_name == "5_cont_UC" || field_name == "6_cont_LC" || field_name == "8_wcont_UC" || field_name == "9_wcont_LC")
                    scale_factor = (temperature < 846.0) ? 567.3 / temperature - 0.062 : 0.732 - 1.35e-4 * temperature;
                  else if (field_name == "3_oc_UC")
                    scale_factor = 0.515 + 1.74 * std::exp(-temperature / 326.5);
                  else if (field_name == "4_oc_LC")
                    scale_factor = 0.2586 + 253.37 / temperature;
                }  
                thermal_diffusivity += volume_fractions[j] * 1e-6 * scale_factor;
              }
            // Compute the thermal conductivity  
            out.thermal_conductivities[i] = thermal_diffusivity * out.specific_heat[i] * out.densities[i];
            // out.thermal_conductivities[i] = thermal_diffusivity;
          }
      }



      template <int dim>
      void
      Lookup<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.declare_entry ("Thermal diffusivities", "0.8e-6",
                           Patterns::List(Patterns::Double (0.)),
                           "List of thermal diffusivities, for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of all compositional fields or only "
                           "those corresponding to chemical compositions. "
                           "If only one value is given, then all use the same value.  "
                           "Units: \\si{\\meter\\squared\\per\\second}.");
      }



      template <int dim>
      void
      Lookup<dim>::parse_parameters (ParameterHandler &prm)
      {
        // Make options file for parsing maps to double arrays
        std::vector<std::string> chemical_field_names = this->introspection().chemical_composition_field_names();
        chemical_field_names.insert(chemical_field_names.begin(),"background");

        std::vector<std::string> compositional_field_names = this->introspection().get_composition_names();
        compositional_field_names.insert(compositional_field_names.begin(),"background");

        Utilities::MapParsing::Options options(chemical_field_names, "Thermal diffusivities");
        options.list_of_allowed_keys = compositional_field_names;

        thermal_diffusivities = Utilities::MapParsing::parse_map_to_double_array(prm.get("Thermal diffusivities"), options);
      }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    namespace ThermalConductivity
    {
#define INSTANTIATE(dim) \
  template class Lookup<dim>;

      ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
    }
  }
}
