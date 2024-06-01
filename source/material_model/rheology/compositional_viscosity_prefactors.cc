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
      std::vector<double>
      CompositionalViscosityPrefactors<dim>::compute_viscosities (const MaterialModel::MaterialModelInputs<dim> &in,
                                                                  const double base_viscosity,
                                                                  const unsigned int composition_index,
                                                                  const unsigned int q) const
      {
        std::vector<double> factored_viscosities(number_of_prefactors, base_viscosity);
        switch (viscosity_prefactor_scheme)
          {
            case none:
            {
              break;
            }
            case hk04_olivine_hydration:
            {
              // We calculate the atomic H/Si ppm (C_OH) at each point to compute the water fugacity of
              // olivine assuming a composition of 90 wt% Forsterite and 10 wt% Fayalite from Hirth
              // and Kohlstaedt 2004 10.1029/138GM06.
              const unsigned int bound_fluid_idx = this->introspection().compositional_index_for_name("bound_fluid");
              const double weight_fraction_H2O = in.composition[q][bound_fluid_idx]; // mass fraction of bound water
              const double weight_fraction_olivine = 1 - weight_fraction_H2O; // mass fraction of olivine
              const double COH = (weight_fraction_H2O/molar_mass_H2O) / (weight_fraction_olivine/molar_mass_olivine) * 1e6; // COH in H / Si ppm
              const double point_water_fugacity = COH / A_H2O * \
                                                  std::exp((activation_energy_H2O + in.pressure[q]*activation_volume_H2O)/ \
                                                           (constants::gas_constant * in.temperature[q]));

              factored_viscosities[0] = base_viscosity*std::pow(point_water_fugacity, -diffusion_water_fugacity_exponents[composition_index]);
              factored_viscosities[1] = base_viscosity*std::pow(point_water_fugacity, -dislocation_water_fugacity_exponents[composition_index]);
              break;
            }
          }
        return factored_viscosities;
      }


      template <int dim>
      void
      CompositionalViscosityPrefactors<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.declare_entry("Number of prefactors", "2",
                          Patterns::List(Patterns::Integer(1)),
                          "The number of flow laws to apply viscosity multiplicative prefactors. "
                          "Default value is 2, for diffusion creep and dislocation creep. ");

        prm.declare_entry ("Water fugacity exponents for diffusion creep", "0.0",
                           Patterns::List(Patterns::Double(0)),
                           "List of water fugacity exponents for diffusion creep for "
                           "background material and compositional fields, for a total of N+1 "
                           "where N is the number of all compositional fields or only those "
                           "corresponding to chemical compositions. Note, the water fugacity exponent "
                           "required by ASPECT for diffusion creep is r/n, where n is the stress exponent "
                           "for diffusion creep, which typically is 1. Units: none.");

        prm.declare_entry ("Water fugacity exponents for dislocation creep", "0.0",
                           Patterns::List(Patterns::Double(0)),
                           "List of water fugacity exponents for dislocation creep for "
                           "background material and compositional fields, for a total of N+1 "
                           "where N is the number of all compositional fields or only those "
                           "corresponding to chemical compositions. Note, the water fugacity exponent "
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
        number_of_prefactors = prm.get_integer ("Number of prefactors");

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

            Utilities::MapParsing::Options options1(chemical_field_names, "Water fugacity exponents for diffusion creep");
            Utilities::MapParsing::Options options2(chemical_field_names, "Water fugacity exponents for dislocation creep");

            options1.list_of_allowed_keys = compositional_field_names;
            options2.list_of_allowed_keys = compositional_field_names;
            diffusion_water_fugacity_exponents = Utilities::MapParsing::parse_map_to_double_array (prm.get("Water fugacity exponents for diffusion creep"),
                                                 options1);
            dislocation_water_fugacity_exponents = Utilities::MapParsing::parse_map_to_double_array (prm.get("Water fugacity exponents for dislocation creep"),
                                                   options2);
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
