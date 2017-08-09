/*
  Copyright (C) 2015 - 2016 by the authors of the ASPECT code.

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


#include <aspect/heating_model/latent_heat_melt.h>


namespace aspect
{
  namespace HeatingModel
  {
    template <int dim>
    void
    LatentHeatMelt<dim>::
    evaluate (const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
              const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
              HeatingModel::HeatingModelOutputs &heating_model_outputs) const
    {
      Assert(heating_model_outputs.heating_source_terms.size() == material_model_inputs.position.size(),
             ExcMessage ("Heating outputs need to have the same number of entries as the material model inputs."));

      const bool use_operator_split = (this->get_parameters().use_operator_splitting);

      const MaterialModel::ReactionRateOutputs<dim> *reaction_rate_out
        = material_model_outputs.template get_additional_output<MaterialModel::ReactionRateOutputs<dim> >();

      for (unsigned int q=0; q<heating_model_outputs.heating_source_terms.size(); ++q)
        {
          heating_model_outputs.heating_source_terms[q] = 0.0;
          heating_model_outputs.lhs_latent_heat_terms[q] = 0.0;
          heating_model_outputs.rates_of_temperature_change[q] = 0.0;

          if (this->introspection().compositional_name_exists("porosity") &&  this->get_timestep_number() > 0)
            {
              const double porosity_idx = this->introspection().compositional_index_for_name("porosity");
              double melting_rate = 0.0;

              if (!use_operator_split)
                {
                  // with melt migration the reaction term is a mass reaction rate,
                  // without melt migration, the reaction term is a constant value in terms of volume,
                  // and we have to scale it to the correct units
                  if (this->include_melt_transport())
                    melting_rate = material_model_outputs.reaction_terms[q][porosity_idx];
                  else
                    melting_rate = material_model_outputs.reaction_terms[q][porosity_idx]
                                   * material_model_outputs.densities[q] / this->get_timestep();

                  heating_model_outputs.heating_source_terms[q] = melting_entropy_change
                                                                  * melting_rate
                                                                  * material_model_inputs.temperature[q];
                }
              else if (use_operator_split && reaction_rate_out != NULL)
                {
                  // if operator splitting is used in the model, we have to use the reaction rates from the
                  // material model outputs instead of the reaction terms
                  AssertThrow (std::isfinite(reaction_rate_out->reaction_rates[q][porosity_idx]),
                               ExcMessage ("You are trying to use reaction rate outputs from the material "
                                           "model to compute the latent heat of melt in an operator splitting solver scheme, "
                                           "but the material model you use does not actually fill these reaction rate outputs."));

                  // with melt migration the reaction term is a mass reaction rate,
                  // without melt migration, the reaction term is a rate of change of the compositional field,
                  // and we have to scale it to the same units to compute the rate of temperature change
                  melting_rate = reaction_rate_out->reaction_rates[q][porosity_idx];

                  // if operator splitting is used in the model, we want the heating rates due to latent heat of melt
                  // to be part of the reactions (not the advection) in the operator split, and they are changes
                  // in temperature rather than changes in energy
                  heating_model_outputs.rates_of_temperature_change[q] = melting_entropy_change
                                                                         * melting_rate
                                                                         * material_model_inputs.temperature[q]
                                                                         / material_model_outputs.specific_heat[q];
                }
              else if (use_operator_split && reaction_rate_out == NULL)
                {
                  // if operator plit is used, but the reaction rate outputs are not there,
                  // fill the rates of temperature change with NaNs, so that an error is thrown
                  // if they are used anywhere
                  heating_model_outputs.rates_of_temperature_change[q] = std::numeric_limits<double>::quiet_NaN();
                }
            }
        }
    }

    template <int dim>
    void
    LatentHeatMelt<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Heating model");
      {
        prm.enter_subsection("Latent heat melt");
        {
          prm.declare_entry ("Melting entropy change", "-300",
                             Patterns::Double (),
                             "The entropy change for the phase transition "
                             "from solid to melt. "
                             "Units: $J/(kg K)$.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    LatentHeatMelt<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Heating model");
      {
        prm.enter_subsection("Latent heat melt");
        {
          melting_entropy_change = prm.get_double ("Melting entropy change");
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
    ASPECT_REGISTER_HEATING_MODEL(LatentHeatMelt,
                                  "latent heat melt",
                                  "Implementation of a standard model for latent heat "
                                  "of melting. This assumes that there is a compositional field "
                                  "called porosity, and it uses the reaction term of this field "
                                  "(the fraction of material that melted in the current time step) "
                                  "multiplied by a constant entropy change for melting all "
                                  "of the material as source term of the heating model.\n"
                                  "If there is no field called porosity, the heating terms are 0.")
  }
}
