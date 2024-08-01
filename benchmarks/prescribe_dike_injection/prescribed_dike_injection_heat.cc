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

#include "prescribed_dike_injection_heat.h"

namespace aspect
{
  namespace HeatingModel
  {
    template <int dim>
    void
    LatentHeatDikeInjection<dim>::
    evaluate (const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
              const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
              HeatingModel::HeatingModelOutputs &heating_model_outputs) const
    {
      AssertThrow(heating_model_outputs.heating_source_terms.size() == material_model_inputs.position.size(),
                  ExcMessage ("Heating outputs need to have the same number of entries as the material "
                              "model inputs."));

      const MaterialModel::PrescribedPlasticDilation<dim>
      *prescribed_dilation =
        (this->get_parameters().enable_prescribed_dilation)
        ? material_model_outputs.template get_additional_output<MaterialModel::PrescribedPlasticDilation<dim> >()
        : nullptr;

      // Add the latent heat source term corresponding to prescribed injection
      // terms in Stokes equations to the rhs of energy conservation equation.
      double dike_injection_fraction = 0.0;

      for (unsigned int q=0; q<heating_model_outputs.heating_source_terms.size(); ++q)
        {
          heating_model_outputs.heating_source_terms[q] = 0.0;
          heating_model_outputs.lhs_latent_heat_terms[q] = 0.0;
          heating_model_outputs.rates_of_temperature_change[q] = 0.0;
          if (prescribed_dilation != nullptr)
            {
              // User-defined or timestep-dependent injection ratio
              if (this->simulator_is_past_initialization())
                dike_injection_fraction = prescribed_dilation->dilation[q] * this->get_timestep();

              if (dike_material_injection_fraction != 0.0)
                dike_injection_fraction = dike_material_injection_fraction;

              // adding the laten heat source team
              heating_model_outputs.heating_source_terms[q] = dike_injection_fraction * 
                                                              prescribed_dilation->dilation[q] * 
                                                              (latent_heat_of_crystallization + 
                                                              (temperature_of_injected_material - 
                                                              material_model_inputs.temperature[q]) * 
                                                              material_model_outputs.densities[q] * 
                                                              material_model_outputs.specific_heat[q]);
            }
        }
    }

    template <int dim>
    void
    LatentHeatDikeInjection<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Heating model");
      {
        prm.enter_subsection("Latent heat dike injection");
        {
          prm.declare_entry ("Latent heat of crystallization", "1.1e9",
                             Patterns::Double(0),
                             "The latent heat of crystallization that is released when material "
                             "is injected into the model. "
                             "Units: J/m$^3$.");
          prm.declare_entry ("Temperature of the injected material", "1273",
                             Patterns::Double(0),
                             "The temperature of the material injected into the model. "
                             "Units: K.");
          prm.declare_entry("Dike material injection fraction", "0.0",
                            Patterns::Double(0),
                            "Amount of new injected material from the dike. Units: none.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    LatentHeatDikeInjection<dim>::
    create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &material_model_outputs) const
    {
      this->get_material_model().create_additional_named_outputs(material_model_outputs);
    }

    template <int dim>
    void
    LatentHeatDikeInjection<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Heating model");
      {
        prm.enter_subsection("Latent heat dike injection");
        {
          latent_heat_of_crystallization = prm.get_double ("Latent heat of crystallization");
          temperature_of_injected_material = prm.get_double ("Temperature of the injected material");
          dike_material_injection_fraction = prm.get_double ("Dike material injection fraction");
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
    ASPECT_REGISTER_HEATING_MODEL(LatentHeatDikeInjection,
                                  "latent heat dike injection",
                                  "Latent heat releases due to the material injection (e.g., melt) into the model. "
                                  "This heating model takes the source term added to the Stokes "
                                  "equation and adds the corresponding source term to the energy "
                                  "equation. This source term includes both the effect of latent "
                                  "heat release upon crystallization and the fact that injected "
                                  "material might have a different temperature.")
  }
}