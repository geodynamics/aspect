/*
  Copyright (C) 2020 - 2023 by the authors of the ASPECT code.

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


#include <aspect/material_model/rheology/frank_kamenetskii.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/utilities.h>

#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/parameter_handler.h>


namespace aspect
{
  namespace MaterialModel
  {
    namespace Rheology
    {
      template <unsigned int>
      FrankKamenetskii<dim>::FrankKamenetskii ()
        = default;



      template <unsigned int>
      double
      FrankKamenetskii<dim>::compute_viscosity (const double temperature,
                                                const unsigned int composition,
                                                const double pressure,
                                                const double density,
                                                const double gravity) const
      {
        const double reference_temperature = this->get_adiabatic_surface_temperature();
        const double reference_pressure = this->get_surface_pressure();
        const double max_depth = this->get_geometry_model().maximal_depth();

        //Frank-Kamenetskii equation with added pressure dependence terms
        const double viscosity_frank_kamenetskii = prefactors_frank_kamenetskii[composition] * std::exp(viscosity_ratios_frank_kamenetskii[composition] * 0.5 * (1.0-temperature/reference_temperature)
                                                   + pressure_prefactors_frank_kamenetskii[composition] * (pressure-reference_pressure)/(density*gravity*max_depth));


        return viscosity_frank_kamenetskii;
      }





      template <unsigned int>
      void
      FrankKamenetskii<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.declare_entry ("Viscosity ratios for Frank Kamenetskii", "15.",
                           Patterns::List(Patterns::Double (0.)),
                           "An adjusted viscosity ratio, $E$, for the viscosity approximation, "
                           "for a total of N+1 values, where N is the number of all compositional fields or only "
                           "those corresponding to chemical compositions. "
                           "If only one value is given, then all use the same value. "
                           "Units: None");
        prm.declare_entry ("Prefactors for Frank Kamenetskii", "1.e21",
                           Patterns::List(Patterns::Double (0.)),
                           "A viscosity prefactor for the viscosity approximation, "
                           "for a total of N+1 values, where N is the number of all compositional fields or only "
                           "those corresponding to chemical compositions. "
                           "If only one value is given, then all use the same value.  Units: None");

        prm.declare_entry ("Pressure prefactors for Frank Kamenetskii", "0.0",
                           Patterns::List(Patterns::Double (0.)),
                           "A prefactor for the pressure term in the viscosity approximation, "
                           "for a total of N+1 values, where N is the number of all compositional fields or only "
                           "those corresponding to chemical compositions. "
                           "If only one value is given, then all use the same value. "
                           "Units: None");
      }



      template <unsigned int>
      void
      FrankKamenetskii<dim>::parse_parameters (ParameterHandler &prm)
      {

        AssertThrow (this->include_adiabatic_heating() == false,
                     ExcMessage("The Frank-Kamenetskii rheology is currently only implemented for "
                                "models without adiabatic heating. Please implement the necessary "
                                "temperature adjustment if you need this feature."));

        AssertThrow (this->get_adiabatic_surface_temperature() > 0.0,
                     ExcMessage("The Frank-Kamenetskii rheology can only be used when the adiabatic "
                                "surface temperature (reference_temperature in equation for viscosity) "
                                "is non-zero."));

        // Retrieve the list of composition names
        std::vector<std::string> compositional_field_names = this->introspection().get_composition_names();

        // Retrieve the list of names of fields that represent chemical compositions, and not, e.g.,
        // plastic strain
        std::vector<std::string> chemical_field_names = this->introspection().chemical_composition_field_names();

        // Establish that a background field is required here
        compositional_field_names.insert(compositional_field_names.begin(), "background");
        chemical_field_names.insert(chemical_field_names.begin(),"background");

        // Make options file for parsing maps to double arrays
        Utilities::MapParsing::Options options(chemical_field_names, "Viscosity ratios for Frank Kamenetskii");
        options.list_of_allowed_keys = compositional_field_names;

        viscosity_ratios_frank_kamenetskii = Utilities::MapParsing::parse_map_to_double_array (prm.get("Viscosity ratios for Frank Kamenetskii"),
                                             options);

        options.property_name = "Prefactors for Frank Kamenetskii";
        prefactors_frank_kamenetskii = Utilities::MapParsing::parse_map_to_double_array(prm.get("Prefactors for Frank Kamenetskii"),
                                                                                        options);

        options.property_name = "Pressure prefactors for Frank Kamenetskii";
        pressure_prefactors_frank_kamenetskii = Utilities::MapParsing::parse_map_to_double_array(prm.get("Pressure prefactors for Frank Kamenetskii"),
                                                options);
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
    template class FrankKamenetskii<dim>; \
  }

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
