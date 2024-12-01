/*
  Copyright (C) 2024 by the authors of the ASPECT code.

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


#include <aspect/material_model/rheology/compositional_viscosity_prefactors.h>
#include <aspect/utilities.h>
#include <aspect/global.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/adiabatic_conditions/interface.h>


#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/parameter_handler.h>
#include <aspect/simulator_signals.h>
#include <aspect/parameters.h>

namespace aspect
{
  namespace MaterialModel
  {
    namespace Rheology
    {
      template <int dim>
      CompositionalViscosityPrefactors<dim>::CompositionalViscosityPrefactors ()
        = default;

      template <int dim>
      double
      CompositionalViscosityPrefactors<dim>::compute_viscosity (const MaterialModel::MaterialModelInputs<dim> &in,
                                                                const double base_viscosity,
                                                                const unsigned int composition_index,
                                                                const unsigned int q,
                                                                const ModifiedFlowLaws &modified_flow_laws) const
      {
        double factored_viscosities = base_viscosity;
        switch (viscosity_prefactor_scheme)
          {
            case none:
            {
              factored_viscosities = base_viscosity;
              break;
            }
            case hk04_olivine_hydration:
            {
              // We calculate the atomic H/Si ppm (C_OH) at each point to compute the water fugacity of
              // olivine assuming a composition of 90 mol% Forsterite and 10 mol% Fayalite from Hirth
              // and Kohlstaedt 2004 10.1029/138GM06.
              const double temperature_for_fugacity = (this->simulator_is_past_initialization())
                                                      ?
                                                      in.temperature[q]
                                                      :
                                                      this->get_adiabatic_conditions().temperature(in.position[q]);
              AssertThrow(temperature_for_fugacity != 0, ExcMessage(
                            "The temperature used in the calculation for determining the water fugacity is zero. "
                            "This is not allowed, because this value is used to divide through. It is probably "
                            "being caused by the temperature being zero somewhere in the model. The relevant "
                            "values for debugging are: temperature (" + Utilities::to_string(in.temperature[q]) +
                            "), and adiabatic temperature ("
                            + Utilities::to_string(this->get_adiabatic_conditions().temperature(in.position[q])) +
                            "). If the adiabatic temperature is 0, double check that you are correctly defining an "
                            " 'Adiabatic conditions model'."));

              const unsigned int bound_fluid_idx = this->introspection().compositional_index_for_name("bound_fluid");
              const double mass_fraction_H2O = std::max(minimum_mass_fraction_water_for_dry_creep[composition_index], in.composition[q][bound_fluid_idx]); // mass fraction of bound water
              const double mass_fraction_olivine = 1 - mass_fraction_H2O; // mass fraction of olivine
              const double COH = (mass_fraction_H2O/molar_mass_H2O) / (mass_fraction_olivine/molar_mass_olivine) * 1e6; // COH in H / Si ppm
              const double point_water_fugacity = COH / A_H2O *
                                                  std::exp((activation_energy_H2O + in.pressure[q]*activation_volume_H2O)/
                                                           (constants::gas_constant * temperature_for_fugacity));
              const double r = modified_flow_laws == diffusion
                               ?
                               -diffusion_water_fugacity_exponents[composition_index]
                               :
                               -dislocation_water_fugacity_exponents[composition_index];
              factored_viscosities = base_viscosity*std::pow(point_water_fugacity, r);
              break;
            }
          }
        return factored_viscosities;
      }


      template <int dim>
      void
      CompositionalViscosityPrefactors<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.declare_entry ("Minimum mass fraction bound water content for fugacity", "6.15e-6",
                           Patterns::List(Patterns::Double(0)),
                           "The minimum water content for the HK04 olivine hydration viscosity "
                           "prefactor scheme. This acts as the cutoff between 'dry' creep and 'wet' creep "
                           "for olivine, and the default value is chosen based on the value reported by "
                           "Hirth & Kohlstaedt 2004. For a mass fraction of bound water beneath this value, "
                           "this value is used instead to compute the water fugacity. Units: \\si{\\kg} / \\si{\\kg} %.");

        prm.declare_entry ("Water fugacity exponents for diffusion creep", "0.0",
                           Patterns::List(Patterns::Double(0)),
                           "List of water fugacity exponents for diffusion creep for "
                           "background material and compositional fields, for a total of N+1 "
                           "where N is the number of all compositional fields or only those "
                           "corresponding to chemical compositions. This is only applied when using the "
                           "Viscosity prefactor scheme 'HK04 olivine hydration'. Note, the water fugacity exponent "
                           "required by ASPECT for diffusion creep is r/n, where n is the stress exponent "
                           "for diffusion creep, which typically is 1. Units: none.");

        prm.declare_entry ("Water fugacity exponents for dislocation creep", "0.0",
                           Patterns::List(Patterns::Double(0)),
                           "List of water fugacity exponents for dislocation creep for "
                           "background material and compositional fields, for a total of N+1 "
                           "where N is the number of all compositional fields or only those "
                           "corresponding to chemical compositions. This is only applied when using the "
                           "Viscosity prefactor scheme 'HK04 olivine hydration'. Note, the water fugacity exponent "
                           "required by ASPECT for dislocation creep is r/n, where n is the stress exponent "
                           "for dislocation creep, which typically is 3.5. Units: none.");

        prm.declare_entry ("Viscosity prefactor scheme", "none",
                           Patterns::Selection("none|HK04 olivine hydration"),
                           "Select what type of viscosity multiplicative prefactor scheme to apply. "
                           "Allowed entries are 'none', and 'HK04 olivine hydration'. HK04 olivine "
                           "hydration calculates the viscosity change due to hydrogen incorporation "
                           "into olivine following Hirth & Kohlstaedt 2004 (10.1029/138GM06). none "
                           "does not modify the viscosity. Units: none.");
      }



      template <int dim>
      void
      CompositionalViscosityPrefactors<dim>::parse_parameters (ParameterHandler &prm)
      {
        if (prm.get ("Viscosity prefactor scheme") == "none")
          viscosity_prefactor_scheme = none;
        if (prm.get ("Viscosity prefactor scheme") == "HK04 olivine hydration")
          {
            // Retrieve the list of compositional names
            std::vector<std::string> compositional_field_names = this->introspection().get_composition_names();
            AssertThrow(this->introspection().compositional_name_exists("bound_fluid"),
                        ExcMessage("The HK04 olivine hydration pre-exponential factor only works if "
                                   "there is a compositional field called bound_fluid."));
            viscosity_prefactor_scheme = hk04_olivine_hydration;

            // Retrieve the list of chemical names
            std::vector<std::string> chemical_field_names = this->introspection().chemical_composition_field_names();

            // Establish that a background field is required here
            compositional_field_names.insert(compositional_field_names.begin(), "background");
            chemical_field_names.insert(chemical_field_names.begin(),"background");

            Utilities::MapParsing::Options options(chemical_field_names, "Water fugacity exponents for diffusion creep");

            options.list_of_allowed_keys = compositional_field_names;
            diffusion_water_fugacity_exponents = Utilities::MapParsing::parse_map_to_double_array (prm.get("Water fugacity exponents for diffusion creep"),
                                                 options);
            dislocation_water_fugacity_exponents = Utilities::MapParsing::parse_map_to_double_array (prm.get("Water fugacity exponents for dislocation creep"),
                                                   options);
            minimum_mass_fraction_water_for_dry_creep = Utilities::MapParsing::parse_map_to_double_array (prm.get("Minimum mass fraction bound water content for fugacity"),
                                                        options);
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
#define INSTANTIATE(dim) \
  namespace Rheology \
  { \
    template class CompositionalViscosityPrefactors<dim>; \
  }

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
