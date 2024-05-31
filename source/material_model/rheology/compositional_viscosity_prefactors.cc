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
      double
      CompositionalViscosityPrefactors<dim>::compute_viscosity (const MaterialModel::MaterialModelInputs<dim> &in,
                                                                const double base_viscosity,
                                                                const unsigned int composition_index,
                                                                const unsigned int q) const
      {
        double prefactor = 1;
        switch (viscosity_prefactor_scheme)
          {
            case none:
            {
              break;
            }
            case hk04_olivine_hydration:
            {
              // We use the atomic H_Si_ratio_ppm composition (C_OH) at each point to compute the water fugacity
              // of olivine assuming a composition of 90 mol% Forsterite and 10 mol% Fayalite from Hirth and
              // Kohlstaedt 2004 10.1029/138GM06.
              const unsigned int bound_fluid_idx = this->introspection().compositional_index_for_name("bound_fluid");
              const double A_H2O = 26 / std::pow(1e6, water_fugacity_exponents[composition_index]);
              const double weight_fraction_H2O = in.composition[q][bound_fluid_idx]; // mass fraction of bound water
              const double weight_fraction_ol = 1 - weight_fraction_H2O; // mass fraction of olivine
              const double COH = (weight_fraction_H2O/M_H2O) / (weight_fraction_ol/M_olivine) * 1e6; // COH in H / Si ppm
              const double point_water_fugacity = COH / A_H2O * \
                                                  std::exp((activation_energy_H2O + in.pressure[q]*activation_volume_H2O)/ \
                                                           (constants::gas_constant * in.temperature[q]));
              prefactor = std::pow(point_water_fugacity, -water_fugacity_exponents[composition_index]);
            }
          }

        return prefactor * base_viscosity;
      }


      template <int dim>
      void
      CompositionalViscosityPrefactors<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.declare_entry ("Exponents for water fugacity", "0.0",
                           Patterns::List(Patterns::Double(0)),
                           "List of viscosity prefactors (i.e., multiplicative factors) "
                           "for background material and compositional fields, for a total of N+1 "
                           "where N is the number of all compositional fields or only those "
                           "corresponding to chemical compositions. Units: none.");

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
        // Retrieve the list of composition names
        if (prm.get ("Viscosity prefactor scheme") == "none")
          viscosity_prefactor_scheme = none;
        if (prm.get ("Viscosity prefactor scheme") == "HK04 olivine hydration")
          {
            std::vector<std::string> compositional_field_names = this->introspection().get_composition_names();
            AssertThrow(this->introspection().compositional_name_exists("bound_fluid"),
                        ExcMessage("The HK04 olivine hydration pre-exponential factor only works if "
                                   "there is a compositional field called bound_fluid."));
            viscosity_prefactor_scheme = hk04_olivine_hydration;

            std::vector<std::string> chemical_field_names = this->introspection().chemical_composition_field_names();

            // Establish that a background field is required here
            compositional_field_names.insert(compositional_field_names.begin(), "background");
            chemical_field_names.insert(chemical_field_names.begin(),"background");

            Utilities::MapParsing::Options options(chemical_field_names, "Exponents for water fugacity");

            // Utilities::MapParsing::Options options(compositional_field_names, "background");
            options.list_of_allowed_keys = compositional_field_names;
            water_fugacity_exponents = Utilities::MapParsing::parse_map_to_double_array (prm.get("Exponents for water fugacity"),
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
