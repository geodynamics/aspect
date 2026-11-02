// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2025 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#include <aspect/heating_model/latent_heat.h>


namespace aspect
{
  namespace HeatingModel
  {
    template <int dim>
    void
    LatentHeat<dim>::
    evaluate (const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
              const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
              HeatingModel::HeatingModelOutputs &heating_model_outputs) const
    {
      Assert(heating_model_outputs.heating_source_terms.size() == material_model_inputs.n_evaluation_points(),
             ExcMessage ("Heating outputs need to have the same number of entries as the material model inputs."));

      for (unsigned int q=0; q<heating_model_outputs.heating_source_terms.size(); ++q)
        {
          heating_model_outputs.heating_source_terms[q] = material_model_inputs.temperature[q]
                                                          * material_model_outputs.densities[q]
                                                          * material_model_outputs.entropy_derivative_pressure[q]
                                                          * (material_model_inputs.velocity[q] * material_model_inputs.pressure_gradient[q]);

          heating_model_outputs.lhs_latent_heat_terms[q] = - material_model_outputs.densities[q]
                                                           * material_model_inputs.temperature[q]
                                                           * material_model_outputs.entropy_derivative_temperature[q];
        }
    }



    template <int dim>
    MaterialModel::MaterialProperties::Property
    LatentHeat<dim>::
    get_required_properties () const
    {
      return MaterialModel::MaterialProperties::entropy_derivative_pressure |
             MaterialModel::MaterialProperties::entropy_derivative_temperature |
             MaterialModel::MaterialProperties::density;
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace HeatingModel
  {
    ASPECT_REGISTER_HEATING_MODEL(LatentHeat,
                                  "latent heat",
                                  "Implementation of a standard model for latent heat.")
  }
}
