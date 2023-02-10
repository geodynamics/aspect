/*
  Copyright (C) 2011 - 2021 by the authors of the ASPECT code.
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


#include <aspect/material_model/equation_of_state/multicomponent_incompressible.h>
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
      MulticomponentIncompressible<dim>::
      evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
               const unsigned int q,
               MaterialModel::EquationOfStateOutputs<dim> &out) const
      {

        const double pressure = in.pressure[q];
        const double temperature = std::max(in.temperature[q], 1.); // temperature can't be zero for correct evaluation
        
        if(use_Birch_Murnaghan)
        {
            for (unsigned int c=0; c < out.densities.size(); ++c)
            {
              const double ak = reference_thermal_expansivities[c]/reference_isothermal_compressibilities[c];
              const double n = isothermal_bulk_modulus_pressure_derivatives[c];
              const double V0 = reference_volumes[c];
              const double B0 = reference_bulk_moduli[c];
              const double P0 = reference_pressures[c];
              const double third_order_correction = (3.*n-1.)/(2.*(n-1.)*(n-1.)) * (pressure - P0) * (pressure - P0) / B0;
              const double f = (1. + (pressure - ak*(temperature - reference_temperatures[c])) * n * reference_isothermal_compressibilities[c]);

              const double V = V0 * (1. + (pressure - P0) / (B0 * (n - 1.))) / pow(1. + (n - 1.) * (pressure - P0) / (2. * B0) + third_order_correction, 1./(n - 1.));
              out.densities[c] = reference_densities[c] * pow(V0 / V, 1.);
              out.thermal_expansion_coefficients[c] = reference_thermal_expansivities[c] / f;
              out.specific_heat_capacities[c] = (isochoric_specific_heats[c] +
                                              (temperature*reference_thermal_expansivities[c] *
                                                  ak * std::pow(f, -1.-(1./n))
                                                  / reference_densities[c]));
              out.compressibilities[c] = reference_isothermal_compressibilities[c]/f;
              out.entropy_derivative_pressure[c] = 0.;
              out.entropy_derivative_temperature[c] = 0.;
            }
        }else if(use_Murnaghan)
        {
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
        }else if(use_compressible_density_only)
        {
            for (unsigned int c=0; c < out.densities.size(); ++c)
            {
                const double ak = reference_thermal_expansivities[c]/reference_isothermal_compressibilities[c];
                const double f = (1. + (pressure - ak*(temperature - reference_temperatures[c])) *
                                isothermal_bulk_modulus_pressure_derivatives[c] *
                                reference_isothermal_compressibilities[c]);
                
                out.thermal_expansion_coefficients[c] = reference_thermal_expansivities[c] / f;

//                 if(c==composition_number_affected && in.temperature[q]>temperature_threshold)
//                 {
//                 out.densities[c] = 4000;
//                 
//                 }else{
                out.densities[c] = reference_densities[c]*std::pow(f, 1./isothermal_bulk_modulus_pressure_derivatives[c]);    
//                 }
                
                if(c==composition_number_affected && in.temperature[q]>temperature_threshold)
                {
                out.specific_heat_capacities[c] = ((isochoric_specific_heats[c] +
                                                (temperature*reference_thermal_expansivities[c] *
                                                    ak * std::pow(f, -1.-(1./isothermal_bulk_modulus_pressure_derivatives[c]))
                                                    / reference_densities[c]))*conductivity_increase_factor);
                }else
                {
                out.specific_heat_capacities[c] = (isochoric_specific_heats[c] +
                                                (temperature*reference_thermal_expansivities[c] *
                                                    ak * std::pow(f, -1.-(1./isothermal_bulk_modulus_pressure_derivatives[c]))
                                                    / reference_densities[c]));                    
                }
                
                //we send back that the model is incompressible, the mass will be balanced with the lithostatic pressure
                out.compressibilities[c] = 0;
                out.entropy_derivative_pressure[c] = 0.;
                out.entropy_derivative_temperature[c] = 0.;
            }
        
        }else if(use_incompressibility){
            for (unsigned int c=0; c < out.densities.size(); ++c)
            {
                out.densities[c] = reference_densities[c] * (1 - reference_thermal_expansivities[c] * (temperature - 274));
                out.thermal_expansion_coefficients[c] = reference_thermal_expansivities[c];
                out.specific_heat_capacities[c] = isochoric_specific_heats[c];
                out.compressibilities[c] = 0.0;
                out.entropy_derivative_pressure[c] = 0.0;
                out.entropy_derivative_temperature[c] = 0.0;
            }
        }else if(use_incompressibility_adiabat){    
        for (unsigned int c=0; c < out.densities.size(); ++c)
          {
                // If adiabatic heating is used, the reference temperature used to calculate density should be the adiabatic
        // temperature at the current position. This definition is consistent with the Extended Boussinesq Approximation.
        const double reference_temperature = (this->include_adiabatic_heating()
                                              ?
                                              this->get_adiabatic_conditions().temperature(in.position[q])
                                              :
                                              reference_temperatures[c]);              
            out.densities[c] = reference_densities[c] * (1 - reference_thermal_expansivities[c] * (in.temperature[q] - reference_temperature));
            out.thermal_expansion_coefficients[c] = reference_thermal_expansivities[c];
            out.specific_heat_capacities[c] = isochoric_specific_heats[c];
            out.compressibilities[c] = 0.0;
            out.entropy_derivative_pressure[c] = 0.0;
            out.entropy_derivative_temperature[c] = 0.0;
          }            
        }else if(use_gms_densities){
            for (unsigned int c=0; c < out.densities.size(); ++c)
            {
                out.densities[c] = reference_densities[c];
                out.thermal_expansion_coefficients[c] = reference_thermal_expansivities[c];
                out.specific_heat_capacities[c] = isochoric_specific_heats[c];
                out.compressibilities[c] = 0.0;
                out.entropy_derivative_pressure[c] = 0.0;
                out.entropy_derivative_temperature[c] = 0.0;
            }            
        }
      }



      template <int dim>
      bool
      MulticomponentIncompressible<dim>::
      is_compressible () const
      {
        if(use_Murnaghan || use_Birch_Murnaghan)
        {
            return true;
        }else{
            return false;
        }
      }



      template <int dim>
      void
      MulticomponentIncompressible<dim>::declare_parameters (ParameterHandler &prm,
                                                             const double default_thermal_expansion)
      {
        prm.declare_entry ("Reference temperatures", "298.15",
                           Patterns::Anything(),
                           "List of reference temperatures $T_0$ for background mantle and compositional fields,"
                           "for a total of N+1 values, where N is the number of compositional fields."
                           "If only one value is given, then all use the same value. Units: \\si{\\kelvin}.");
        prm.declare_entry ("Densities", "3300.",
                           Patterns::Anything(),
                           "List of densities for background mantle and compositional fields,"
                           "for a total of N+1 values, where N is the number of compositional fields."
                           "If only one value is given, then all use the same value. "
                           "Units: \\si{\\kilogram\\per\\meter\\cubed}.");
        prm.declare_entry ("Reference isothermal compressibilities", "4e-12",
                           Patterns::Anything(),
                           "List of isothermal compressibilities for background mantle and compositional fields,"
                           "for a total of N+1 values, where N is the number of compositional fields."
                           "If only one value is given, then all use the same value. "
                           "Units: \\si{\\per\\pascal}.");
        prm.declare_entry ("Isothermal bulk modulus pressure derivatives", "4.",
                           Patterns::Anything(),
                           "List of isothermal pressure derivatives of the bulk moduli for background mantle and compositional fields,"
                           "for a total of N+1 values, where N is the number of compositional fields."
                           "If only one value is given, then all use the same value. "
                           "Units: [].");
        prm.declare_entry ("Thermal expansivities", std::to_string(default_thermal_expansion),
                           Patterns::Anything(),
                           "List of thermal expansivities for background mantle and compositional fields,"
                           "for a total of N+1 values, where N is the number of compositional fields."
                           "If only one value is given, then all use the same value. Units: \\si{\\per\\kelvin}.");
        prm.declare_entry ("Heat capacities", "1250.",
                           Patterns::Anything(),
                           "List of isochoric specific heats $C_v$ for background mantle and compositional fields,"
                           "for a total of N+1 values, where N is the number of compositional fields."
                           "If only one value is given, then all use the same value. "
                           "Units: \\si{\\joule\\per\\kelvin\\per\\kilogram}.");

        prm.declare_entry ("Reference pressures", "1e5",
                           Patterns::Anything(),
                           "List of reference volumes $C_v$ for background mantle and compositional fields,"
                           "for a total of N+1 values, where N is the number of compositional fields."
                           "If only one value is given, then all use the same value. "
                           "Units: \\si{\\joule\\per\\Pascal}.");
        prm.declare_entry ("Reference volumes", "1e-5",
                           Patterns::Anything(),
                           "List of reference volumes $C_v$ for background mantle and compositional fields,"
                           "for a total of N+1 values, where N is the number of compositional fields."
                           "If only one value is given, then all use the same value. "
                           "Units: \\si{\\joule\\per\\Pascal}.");
        prm.declare_entry ("Reference bulk moduli", "163e9",
                           Patterns::Anything(),
                           "List of bulk moduli $C_v$ for background mantle and compositional fields,"
                           "for a total of N+1 values, where N is the number of compositional fields."
                           "If only one value is given, then all use the same value. "
                           "Units: \\si{\\Pa}.");                           
        prm.declare_entry("Use Birch-Murnaghan", "false",
                            Patterns::Bool(),
                             "Define if the density is incompressible of compressible and how it should be handled for the mass conservation");                            
        prm.declare_entry("Use Murnaghan", "false",
                            Patterns::Bool(),
                             "Define if the density is incompressible of compressible and how it should be handled for the mass conservation"); 
        prm.declare_entry("Use compressible density only", "false",
                            Patterns::Bool(),
                             "Define if the density is incompressible of compressible and how it should be handled for the mass conservation");  
        prm.declare_entry("Use incompressibility", "false",
                            Patterns::Bool(),
                             "Define if the density is incompressible of compressible and how it should be handled for the mass conservation");
        prm.declare_entry("Use incompressibility with adiabat", "false",
                            Patterns::Bool(),
                             "Define if the density is incompressible of compressible and how it should be handled for the mass conservation");
        prm.declare_entry("Use GMS densities", "false",
                            Patterns::Bool(),
                             "Define if the density is incompressible of compressible and how it should be handled for the mass conservation");         

        prm.declare_entry("Use conductivity temperature dependent", "false",
                            Patterns::Bool(),
                             "change the  spcific heat by a factor so it applies later to the conductivity");
        prm.declare_entry ("Composition number affected", "5", Patterns::Double (0.),
                             "Composition affected by the change of conductivity");
        prm.declare_entry ("Temperature of activation", "1000", Patterns::Double (0.),
                             "Temperature at which conductivity factor applies in Kelvin"); 
        prm.declare_entry ("Conductivity increase factor", "10", Patterns::Double (0.),
                             "Factor of increase of conductivity");        
        
      }



      template <int dim>
      void
      MulticomponentIncompressible<dim>::parse_parameters (ParameterHandler &prm,
                                                           const std::unique_ptr<std::vector<unsigned int>> &expected_n_phases_per_composition)
      {
       // Establish that a background field is required here
        const bool has_background_field = true;

        // Retrieve the list of composition names
        const std::vector<std::string> list_of_composition_names = this->introspection().get_composition_names();

        reference_volumes = Utilities::parse_map_to_double_array (prm.get("Reference volumes"),
                                                                       list_of_composition_names,
                                                                       has_background_field,
                                                                       "Reference volumes",
                                                                       true,
                                                                       expected_n_phases_per_composition);
        // Parse multicomponent properties
        reference_bulk_moduli = Utilities::parse_map_to_double_array (prm.get("Reference bulk moduli"),
                                                                       list_of_composition_names,
                                                                       has_background_field,
                                                                       "Reference bulk moduli",
                                                                       true,
                                                                       expected_n_phases_per_composition);
        reference_pressures = Utilities::parse_map_to_double_array (prm.get("Reference pressures"),
                                                                       list_of_composition_names,
                                                                       has_background_field,
                                                                       "Reference pressures",
                                                                       true,
                                                                       expected_n_phases_per_composition);
        reference_temperatures = Utilities::parse_map_to_double_array (prm.get("Reference temperatures"),
                                                                       list_of_composition_names,
                                                                       has_background_field,
                                                                       "Reference temperatures",
                                                                       true,
                                                                       expected_n_phases_per_composition);

        reference_densities = Utilities::parse_map_to_double_array (prm.get("Densities"),
                                                                    list_of_composition_names,
                                                                    has_background_field,
                                                                    "Densities",
                                                                    true,
                                                                    expected_n_phases_per_composition);

        reference_isothermal_compressibilities = Utilities::parse_map_to_double_array (prm.get("Reference isothermal compressibilities"),
                                                                                       list_of_composition_names,
                                                                                       has_background_field,
                                                                                       "Reference isothermal compressibilities",
                                                                                       true,
                                                                                       expected_n_phases_per_composition);


        isothermal_bulk_modulus_pressure_derivatives = Utilities::parse_map_to_double_array (prm.get("Isothermal bulk modulus pressure derivatives"),
                                                       list_of_composition_names,
                                                       has_background_field,
                                                       "Isothermal bulk modulus pressure derivatives",
                                                       true,
                                                       expected_n_phases_per_composition);


        reference_thermal_expansivities = Utilities::parse_map_to_double_array (prm.get("Thermal expansivities"),
                                                                                list_of_composition_names,
                                                                                has_background_field,
                                                                                "Thermal expansivities",
                                                                                true,
                                                                                expected_n_phases_per_composition);


        isochoric_specific_heats = Utilities::parse_map_to_double_array (prm.get("Heat capacities"),
                                                                         list_of_composition_names,
                                                                         has_background_field,
                                                                         "Heat capacities",
                                                                         true,
                                                                         expected_n_phases_per_composition);
        use_Birch_Murnaghan = prm.get_bool ("Use Birch-Murnaghan");
        use_Murnaghan = prm.get_bool ("Use Murnaghan");
        use_compressible_density_only = prm.get_bool ("Use compressible density only");
        use_incompressibility = prm.get_bool ("Use incompressibility");
        use_incompressibility_adiabat=prm.get_bool ("Use incompressibility with adiabat");
        use_gms_densities==prm.get_bool ("Use GMS densities");
        
        use_conductivity_temperature_dependent = prm.get_bool("Use conductivity temperature dependent");
        composition_number_affected = prm.get_double("Composition number affected");
        temperature_threshold = prm.get_double("Temperature of activation"); 
        conductivity_increase_factor = prm.get_double("Conductivity increase factor");          

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
  template class MulticomponentIncompressible<dim>;

      ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
    }
  }
}