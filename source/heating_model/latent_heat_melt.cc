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
  along with ASPECT; see the file doc/COPYING.  If not see
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

      for (unsigned int q=0; q<heating_model_outputs.heating_source_terms.size(); ++q)
        {
          if (this->introspection().compositional_name_exists("porosity") &&  this->get_timestep_number() > 0)
            {
              // with melt migration the reaction term is a mass reaction rate
              const double porosity_idx = this->introspection().compositional_index_for_name("porosity");
              const double latent_heat = melting_entropy_change * material_model_outputs.reaction_terms[q][porosity_idx];
              heating_model_outputs.heating_source_terms[q] = latent_heat
                                                              * material_model_inputs.temperature[q];

              // without melt migration, the reaction term is a constant value in terms of volume,
              // and we have to scale it to the correct units
              if (!this->include_melt_transport())
                heating_model_outputs.heating_source_terms[q] *= material_model_outputs.densities[q] / this->get_timestep();
            }
          else
            heating_model_outputs.heating_source_terms[q] = 0.0;

          heating_model_outputs.lhs_latent_heat_terms[q] = 0.0;
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
