/*
  Copyright (C) 2011 - 2025 by the authors of the ASPECT code.

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


#include <aspect/heating_model/tidal_heating.h>


namespace aspect
{
  namespace HeatingModel
  {
    template <int dim>
    void
    TidalHeating<dim>::
    evaluate (const MaterialModel::MaterialModelInputs<dim> &/*material_model_inputs*/,
              const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
              HeatingModel::HeatingModelOutputs &heating_model_outputs) const
    {
      /** *
       * H equation is from Tobie et al. (2003) (https://doi.org/10.1029/2003JE002099)
       * H= 2*(viscosity)*(time-averaged tidal strain rate)^2/(1+((viscosity)*(tidal frequency)/(elastic shear modulus))^2))
       * viscosity = material_model_outputs.viscosities
       * time-averaged strain rate = tidal_strain_rate
       * tidal frequency = tidal_frequency
       * elastic shear modulus = elastic_shear_modulus
       */
      for (unsigned int q=0; q<heating_model_outputs.heating_source_terms.size(); ++q)
        {
          heating_model_outputs.heating_source_terms[q] = 2 * material_model_outputs.viscosities[q] * tidal_strain_rate * tidal_strain_rate
                                                          / ( 1 + ( (tidal_frequency * material_model_outputs.viscosities[q]) / elastic_shear_modulus ) * ( (tidal_frequency * material_model_outputs.viscosities[q]) / elastic_shear_modulus ) );
          heating_model_outputs.lhs_latent_heat_terms[q] = 0.0;
        }
    }



    template <int dim>
    MaterialModel::MaterialProperties::Property
    TidalHeating<dim>::
    get_required_properties () const
    {
      return MaterialModel::MaterialProperties::viscosity;
    }



    template <int dim>
    void
    TidalHeating<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Heating model");
      {
        prm.enter_subsection("Tidal heating");
        {
          prm.declare_entry ("Tidal frequency", "2.048e-5",
                             Patterns::Double (0),
                             "The orbital/tidal frequency that produces the heating. "
                             "Default value is the diurnal tidal frequency of Europa, ~3.551 days. "
                             "Units: 1/s.");
          prm.declare_entry ("Elastic shear modulus", "3.3e9",
                             Patterns::Double (0),
                             "Elastic shear modulus of the material. "
                             "For simplicity, this parameter will be used even if elasticity is set in the material model. "
                             "Default value is for Europa's icy shell. "
                             "Units: Pa.");
          prm.declare_entry ("Tidal strain rate", "2e-10",
                             Patterns::Double (0),
                             "Time-averaged strain rate of the diurnal tide. "
                             "Default value is for the diurnal tide of Europa. "
                             "Units: 1/s.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    TidalHeating<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Heating model");
      {
        prm.enter_subsection("Tidal heating");
        {
          tidal_frequency = prm.get_double ("Tidal frequency");
          elastic_shear_modulus = prm.get_double ("Elastic shear modulus");
          tidal_strain_rate = prm.get_double ("Tidal strain rate");
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
    ASPECT_REGISTER_HEATING_MODEL(TidalHeating,
                                  "tidal heating",
                                  "A tidal heating implementation related to diurnal tides. "
                                  "The equation currently ignores regional (radial/lateral) changes. "
                                  "This equation is the Eq.12 from Tobie et al. (2003) (https://doi.org/10.1029/2003JE002099). "
                                  "Unit: W/m^3.")
  }
}
