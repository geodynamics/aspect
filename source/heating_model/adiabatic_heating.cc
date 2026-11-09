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


#include <aspect/heating_model/adiabatic_heating.h>
#include <aspect/gravity_model/interface.h>


namespace aspect
{
  namespace HeatingModel
  {
    template <int dim>
    bool
    AdiabaticHeating<dim>::use_simplified_adiabatic_heating() const
    {
      return simplified_adiabatic_heating;
    }

    template <int dim>
    void
    AdiabaticHeating<dim>::
    evaluate (const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
              const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
              HeatingModel::HeatingModelOutputs &heating_model_outputs) const
    {
      Assert(heating_model_outputs.heating_source_terms.size() == material_model_inputs.n_evaluation_points(),
             ExcMessage ("Heating outputs need to have the same number of entries as the material model inputs."));

      for (unsigned int q=0; q<heating_model_outputs.heating_source_terms.size(); ++q)
        {
          if (!simplified_adiabatic_heating)
            heating_model_outputs.heating_source_terms[q] = (material_model_inputs.velocity[q] * material_model_inputs.pressure_gradient[q])
                                                            * material_model_outputs.thermal_expansion_coefficients[q]
                                                            * material_model_inputs.temperature[q];
          else
            heating_model_outputs.heating_source_terms[q] = (material_model_inputs.velocity[q]
                                                             * this->get_gravity_model().gravity_vector(material_model_inputs.position[q]))
                                                            * material_model_outputs.thermal_expansion_coefficients[q]
                                                            * material_model_inputs.temperature[q]
                                                            * material_model_outputs.densities[q];

          heating_model_outputs.lhs_latent_heat_terms[q] = 0.0;
        }
    }



    template <int dim>
    MaterialModel::MaterialProperties::Property
    AdiabaticHeating<dim>::
    get_required_properties () const
    {
      return MaterialModel::MaterialProperties::thermal_expansion_coefficient |
             MaterialModel::MaterialProperties::density;
    }



    template <int dim>
    void
    AdiabaticHeating<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Heating model");
      {
        prm.enter_subsection("Adiabatic heating");
        {
          prm.declare_entry ("Use simplified adiabatic heating", "false",
                             Patterns::Bool (),
                             "A flag indicating whether the adiabatic heating should be simplified "
                             "from $\\alpha T (\\mathbf u \\cdot \\nabla p)$ to "
                             "$ \\alpha \\rho T (\\mathbf u \\cdot \\mathbf g) $.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    AdiabaticHeating<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Heating model");
      {
        prm.enter_subsection("Adiabatic heating");
        {
          simplified_adiabatic_heating = prm.get_bool ("Use simplified adiabatic heating");
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
    ASPECT_REGISTER_HEATING_MODEL(AdiabaticHeating,
                                  "adiabatic heating",
                                  "Implementation of a standard and a simplified model of "
                                  "adiabatic heating.")
  }
}
