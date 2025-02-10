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


#include <aspect/material_model/reaction_model/tian2019_solubility.h>
#include <deal.II/base/parameter_handler.h>


namespace aspect
{
  namespace MaterialModel
  {
    namespace ReactionModel
    {


      template <int dim>
      double
      Tian2019Solubility<dim>::
      melt_fraction (const MaterialModel::MaterialModelInputs<dim> &in,
                     const unsigned int porosity_idx,
                     unsigned int q) const
      {
        double melt_fractions;
        // The bound fluid content is calculated using parametrized phase
        // diagrams for four different rock types: sediment, MORB, gabbro, and
        // peridotite.
        const unsigned int bound_fluid_idx = this->introspection().compositional_index_for_name("bound_fluid");
        const unsigned int sediment_idx = this->introspection().compositional_index_for_name("sediment");
        const unsigned int MORB_idx = this->introspection().compositional_index_for_name("MORB");
        const unsigned int gabbro_idx = this->introspection().compositional_index_for_name("gabbro");
        const unsigned int peridotite_idx = this->introspection().compositional_index_for_name("peridotite");

        // Initialize a vector that stores the compositions (mass fractions) for
        // the four different rock compositions,
        std::vector<double> tracked_rock_mass_fractions(4);
        tracked_rock_mass_fractions[0] = (in.composition[q][peridotite_idx]);
        tracked_rock_mass_fractions[1] = (in.composition[q][gabbro_idx]);
        tracked_rock_mass_fractions[2] = (in.composition[q][MORB_idx]);
        tracked_rock_mass_fractions[3] = (in.composition[q][sediment_idx]);

        // The bound water content (water within the solid phase) for the four different rock types
        std::vector<double> tian_eq_bound_water_content = tian_equilibrium_bound_water_content(in, q);

        // average the water content between the four different rock types
        double average_eq_bound_water_content = MaterialUtilities::average_value (tracked_rock_mass_fractions, tian_eq_bound_water_content, MaterialUtilities::arithmetic);

        // The fluid volume fraction in equilibrium with the solid (stored in the melt_fractions vector)
        // is equal to the sum of the porosity and the change in bound fluid content
        // (current bound fluid - updated average bound fluid).
        melt_fractions = std::max(in.composition[q][bound_fluid_idx] + in.composition[q][porosity_idx] - average_eq_bound_water_content, 0.0);

        return melt_fractions;
      }



      template <int dim>
      std::vector<double>
      Tian2019Solubility<dim>::
      tian_equilibrium_bound_water_content (const MaterialModel::MaterialModelInputs<dim> &in,
                                            unsigned int q) const
      {
        // Create arrays that will store the values of the polynomials at the current pressure
        std::vector<double> LR_values(4);
        std::vector<double> csat_values(4);
        std::vector<double> Td_values(4);

        // Loop over the four rock types (peridotite, gabbro, MORB, sediment) and the polynomial
        // coefficients to fill the vectors defined above. The polynomials for LR are defined in
        // equations 13, B2, B10, and B18. csat polynomials are defined in equations 14, B1, B9, and B17.
        // Td polynomials are defined in equations 15, B3, B11, and B19.
        for (unsigned int i = 0; i<devolatilization_enthalpy_changes.size(); ++i)
          {
            // Pressure, which must be in GPa for the parametrization, or GPa^-1. The polynomials for each lithology
            // breaks down above certain pressures, make sure that we cap the pressure just before this break down.
            // Introduce minimum pressure to avoid a division by 0.
            const double minimum_pressure = 1e-12;
            const double pressure = std::min(std::max(minimum_pressure, in.pressure[q]/1.e9), pressure_cutoffs[i]);
            const double inverse_pressure = 1.0/pressure;
            for (unsigned int j = 0; j<devolatilization_enthalpy_changes[i].size(); ++j)
              {
#if DEAL_II_VERSION_GTE(9, 6, 0)
                LR_values[i] += devolatilization_enthalpy_changes[i][j] * Utilities::pow(inverse_pressure, devolatilization_enthalpy_changes[i].size() - 1 - j);
#else
                LR_values[i] += devolatilization_enthalpy_changes[i][j] * std::pow(inverse_pressure, devolatilization_enthalpy_changes[i].size() - 1 - j);
#endif
              }

            for (unsigned int j = 0; j<water_mass_fractions[i].size(); ++j)
              {
#if DEAL_II_VERSION_GTE(9, 6, 0)
                csat_values[i] += i==3 ? water_mass_fractions[i][j] * Utilities::pow(std::log10(pressure), water_mass_fractions[i].size() - 1 - j) :\
                                  water_mass_fractions[i][j] * Utilities::pow(pressure, water_mass_fractions[i].size() - 1 - j);
#else
                csat_values[i] += i==3 ? water_mass_fractions[i][j] * std::pow(std::log10(pressure), water_mass_fractions[i].size() - 1 - j) :\
                                  water_mass_fractions[i][j] * std::pow(pressure, water_mass_fractions[i].size() - 1 - j);
#endif
              }

            for (unsigned int j = 0; j<devolatilization_onset_temperatures[i].size(); ++j)
              {
#if DEAL_II_VERSION_GTE(9, 6, 0)
                Td_values[i] += devolatilization_onset_temperatures[i][j] * Utilities::pow(pressure, devolatilization_onset_temperatures[i].size() - 1 - j);
#else
                Td_values[i] += devolatilization_onset_temperatures[i][j] * std::pow(pressure, devolatilization_onset_temperatures[i].size() - 1 - j);
#endif
              }
          }

        // Create an array for the equilibrium bound water content that is calculated from these polynomials
        std::vector<double> eq_bound_water_content(4);

        // Define the maximum bound water content allowed for the four different rock compositions
        std::vector<double> max_bound_water_content = {tian_max_peridotite_water, tian_max_gabbro_water, tian_max_MORB_water, tian_max_sediment_water};

        // Loop over all rock compositions and fill the equilibrium bound water content, divide by 100 to convert
        // from percentage to fraction (equation 1)
        for (unsigned int k = 0; k<LR_values.size(); ++k)
          {
            eq_bound_water_content[k] = (std::min(std::exp(csat_values[k]) * \
                                                  std::exp(std::exp(LR_values[k]) * (1/in.temperature[q] - 1/Td_values[k])), \
                                                  max_bound_water_content[k]) / 100.0);
          }
        return eq_bound_water_content;
      }



      template <int dim>
      void
      Tian2019Solubility<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.declare_entry ("Maximum weight percent water in sediment", "3",
                           Patterns::Double (0),
                           "The maximum allowed weight percent that the sediment composition can hold.");
        prm.declare_entry ("Maximum weight percent water in MORB", "2",
                           Patterns::Double (0),
                           "The maximum allowed weight percent that the sediment composition can hold.");
        prm.declare_entry ("Maximum weight percent water in gabbro", "1",
                           Patterns::Double (0),
                           "The maximum allowed weight percent that the sediment composition can hold.");
        prm.declare_entry ("Maximum weight percent water in peridotite", "8",
                           Patterns::Double (0),
                           "The maximum allowed weight percent that the sediment composition can hold.");
      }


      template <int dim>
      void
      Tian2019Solubility<dim>::parse_parameters (ParameterHandler &prm)
      {
        AssertThrow(this->introspection().compositional_name_exists("sediment"),
                    ExcMessage("The Tian approximation only works "
                               "if there is a compositional field called sediment."));
        AssertThrow(this->introspection().compositional_name_exists("MORB"),
                    ExcMessage("The Tian approximation only works "
                               "if there is a compositional field called MORB."));
        AssertThrow(this->introspection().compositional_name_exists("gabbro"),
                    ExcMessage("The Tian approximation only works "
                               "if there is a compositional field called gabbro."));
        AssertThrow(this->introspection().compositional_name_exists("peridotite"),
                    ExcMessage("The Tian approximation only works "
                               "if there is a compositional field called peridotite."));
        tian_max_peridotite_water         = prm.get_double ("Maximum weight percent water in peridotite");
        tian_max_gabbro_water             = prm.get_double ("Maximum weight percent water in gabbro");
        tian_max_MORB_water               = prm.get_double ("Maximum weight percent water in MORB");
        tian_max_sediment_water           = prm.get_double ("Maximum weight percent water in sediment");
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
  namespace ReactionModel \
  { \
    template class Tian2019Solubility<dim>; \
  }

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
