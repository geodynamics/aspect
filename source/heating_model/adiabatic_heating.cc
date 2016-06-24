/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include <aspect/heating_model/adiabatic_heating.h>


namespace aspect
{
  namespace HeatingModel
  {
    template <int dim>
    void
    AdiabaticHeating<dim>::
    evaluate (const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
              const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
              HeatingModel::HeatingModelOutputs &heating_model_outputs) const
    {
      Assert(heating_model_outputs.heating_source_terms.size() == material_model_inputs.position.size(),
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
                                  "Implementation of a standard and a simplified model of"
                                  "adiabatic heating.")
  }
}
