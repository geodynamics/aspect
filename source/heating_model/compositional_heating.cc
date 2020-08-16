/*
  Copyright (C) 2017 - 2020 by the authors of the ASPECT code.

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



#include <aspect/heating_model/compositional_heating.h>
#include <aspect/utilities.h>
#include <aspect/global.h>

namespace aspect
{
  namespace HeatingModel
  {
    template <int dim>
    void
    CompositionalHeating<dim>::
    evaluate (const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
              const MaterialModel::MaterialModelOutputs<dim> &/*material_model_outputs*/,
              HeatingModel::HeatingModelOutputs &heating_model_outputs) const
    {
      for (unsigned int q = 0; q < heating_model_outputs.heating_source_terms.size(); ++q)
        {
          // Compute compositional volume fractions
          const std::vector<double> volume_fractions = MaterialModel::MaterialUtilities::compute_volume_fractions(material_model_inputs.composition[q],
                                                       fields_used_in_heat_production_averaging);

          // Calculate average compositional heat production
          double compositional_heat_production = 0.;
          for (unsigned int c=0; c < volume_fractions.size(); ++c)
            if (c>0 || use_background_field_for_heat_production_averaging)
              compositional_heat_production += volume_fractions[c] * heating_values[c];

          heating_model_outputs.heating_source_terms[q] = compositional_heat_production;
          heating_model_outputs.lhs_latent_heat_terms[q] = 0.0;
        }
    }


    template <int dim>
    void
    CompositionalHeating<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Heating model");
      {
        prm.enter_subsection("Compositional heating");
        {
          prm.declare_entry("Compositional heating values","0.",
                            Patterns::List (Patterns::Double(0.)),
                            "List of heat production per unit volume values for "
                            "background and compositional fields, for a total of "
                            "N+1 values, where the first value corresponds to the "
                            "background material, and N is the number of compositional fields. "
                            "Units: \\si{\\watt\\per\\meter\\cubed}.");
          prm.declare_entry ("Use compositional field for heat production averaging", "1",
                             Patterns::List(Patterns::Integer(0,1)),
                             "A list of integers with as many entries as compositional fields plus one. "
                             "The first entry corresponds to the background material, each following "
                             "entry corresponds to a particular compositional field. If the entry for "
                             "a field is '1' this field is considered during the computation of volume "
                             "fractions, if it is '0' the field is ignored. This is useful "
                             "if some compositional fields are used to track properties like "
                             "finite strain that should not contribute to heat production. The first "
                             "entry determines whether the background field contributes to heat production "
                             "or not (essentially similar to setting its 'Compositional heating values' to "
                             "zero, but included for consistency in the length of the input lists).");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    CompositionalHeating<dim>::parse_parameters (ParameterHandler &prm)
    {
      // increment by one for background:
      const unsigned int n_fields = this->n_compositional_fields() + 1;

      prm.enter_subsection("Heating model");
      {
        prm.enter_subsection("Compositional heating");
        {
          const std::vector<int> used_fields = Utilities::possibly_extend_from_1_to_N (
                                                 Utilities::string_to_int(
                                                   Utilities::split_string_list(
                                                     prm.get("Use compositional field for heat production averaging"))),
                                                 n_fields,
                                                 "Use compositional field for heat production averaging");

          use_background_field_for_heat_production_averaging = static_cast<bool>(used_fields[0]);

          fields_used_in_heat_production_averaging.resize(this->n_compositional_fields());
          for (unsigned int i=0; i<fields_used_in_heat_production_averaging.size(); ++i)
            fields_used_in_heat_production_averaging[i] = static_cast<bool>(used_fields[i+1]);

          heating_values = Utilities::possibly_extend_from_1_to_N (
                             Utilities::string_to_double(
                               Utilities::split_string_list(
                                 prm.get("Compositional heating values"))),
                             n_fields,
                             "Compositional heating values");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

  }
}

// explicit instantiations
namespace aspect
{
  namespace HeatingModel
  {
    ASPECT_REGISTER_HEATING_MODEL(CompositionalHeating,
                                  "compositional heating",
                                  "Implementation of a model in which magnitude of internal heat production "
                                  "is determined from fixed values assigned to each compositional "
                                  "field. These values are interpreted as having units "
                                  "\\si{\\watt\\per\\meter\\cubed}.")
  }
}
