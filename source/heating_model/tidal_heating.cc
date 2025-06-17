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

#include <aspect/simulator.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/geometry_model/chunk.h>
#include <aspect/geometry_model/two_merged_chunks.h>
#include <aspect/geometry_model/ellipsoidal_chunk.h>
#include <aspect/geometry_model/sphere.h>
#include <aspect/geometry_model/spherical_shell.h>

namespace aspect
{
  namespace HeatingModel
  {
    template <int dim>
    void
    TidalHeating<dim>::initialize ()
    {
      if (strain_rate_distribution == latitudinal_variation)
        {
          // Only spherical geometries will be accepted to have heating rate varying by latitude.
          AssertThrow((Plugins::plugin_type_matches<GeometryModel::Chunk<dim>>(this->get_geometry_model()) ||
                       Plugins::plugin_type_matches<GeometryModel::TwoMergedChunks<dim>>(this->get_geometry_model()) ||
                       Plugins::plugin_type_matches<GeometryModel::EllipsoidalChunk<dim>>(this->get_geometry_model()) ||
                       Plugins::plugin_type_matches<GeometryModel::Sphere<dim>>(this->get_geometry_model()) ||
                       Plugins::plugin_type_matches<GeometryModel::SphericalShell<dim>>(this->get_geometry_model())) &&
                      dim==3,
                      ExcMessage("Latitudinal variation of strain rate can only be used with 3-dimensional geometry models that "
                                 "have either a spherical or ellipsoidal natural coordinate system."));
        }
    }



    template <int dim>
    void
    TidalHeating<dim>::
    evaluate (const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
              const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
              HeatingModel::HeatingModelOutputs &heating_model_outputs) const
    {
      /** *
       * H equation is from Tobie et al. (2003) (https://doi.org/10.1029/2003JE002099)
       * H= 2*(viscosity)*(time-averaged tidal strain rate)^2/(1+((viscosity)*(tidal frequency)/(elastic shear modulus))^2))
       * viscosity = material_model_outputs.viscosities
       * time-averaged strain rate at certain location = local_strain_rate
       * tidal frequency = tidal_frequency
       * elastic shear modulus = elastic_shear_modulus
       * If 'Use latitudinal variation of strain rate' is true, local_strain_rate is calculated
       * with one period cosine function having maximum tidal strain rate at poles and minimum tidal strain rate at equator.
       * This variation of tidal strain rate follows as (maximum_tidal_strain_rate - minimum_tidal_strain_rate)*cos(theta/2)+(maximum_tidal_strain_rate+minimum_tidal_strain_rate)/2
       * where theta is polar angle from spherical coordinates.
       * If false, constant_tidal_strain_rate is used.
       */
      for (unsigned int q=0; q<heating_model_outputs.heating_source_terms.size(); ++q)
        {
          double local_tidal_strain_rate = 0.;
          if (strain_rate_distribution == latitudinal_variation)
            {
              const Point<dim> position = material_model_inputs.position[q];
              const double theta = Utilities::Coordinates::cartesian_to_spherical_coordinates(position)[dim-1];
              local_tidal_strain_rate = ( maximum_tidal_strain_rate - minimum_tidal_strain_rate ) / 2. * std::cos( theta * 2. ) + ( maximum_tidal_strain_rate + minimum_tidal_strain_rate ) / 2.;
            }
          else if (strain_rate_distribution == constant)
            {
              local_tidal_strain_rate = constant_tidal_strain_rate;
            }
          heating_model_outputs.heating_source_terms[q] = 2. * material_model_outputs.viscosities[q] * local_tidal_strain_rate * local_tidal_strain_rate
                                                          / ( 1. + ( (tidal_frequency * material_model_outputs.viscosities[q]) / elastic_shear_modulus ) * ( (tidal_frequency * material_model_outputs.viscosities[q]) / elastic_shear_modulus ) );
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
          prm.declare_entry ("Constant tidal strain rate", "2e-10",
                             Patterns::Double (0),
                             "Time-averaged strain rate by a lunar diurnal tide that is simplified to be constant regardless of location. "
                             "Default value is for the diurnal tide of Europa from Tobie et al. (2003). "
                             "Units: 1/s.");
          prm.declare_entry ("Custom distribution of tidal strain rate", "constant",
                             Patterns::Selection("constant|latitudinal variation"),
                             "Choose how the time-averaged tidal strain rate is distributed. "
                             "If 'constant', the tidal strain rate is fixed to 'Constant tidal strain rate'. "
                             "If 'latitudinal variation', 'Maximum tidal strain rate' and 'Minimum tidal strain rate' are used.");
          prm.declare_entry ("Maximum tidal strain rate", "2.81e-10",
                             Patterns::Double (0),
                             "Maximum time-averaged tidal strain rate by lunar diurnal tide at the poles. "
                             "This parameter will be used when 'Custom distribution of tidal strain rate' is 'latitudinal variation'. "
                             "Default value is for Europa at pole from Nimmo et al. (2007). "
                             "Units: 1/s.");
          prm.declare_entry ("Minimum tidal strain rate", "1.67e-10",
                             Patterns::Double (0),
                             "Minimum time-averaged tidal strain rate by lunar diurnal tide at the equator. "
                             "This parameter will be used when 'Custom distribution of tidal strain rate' is 'latitudinal variation'. "
                             "Default value is for Europa at Equator from Nimmo et al. (2007). "
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
          constant_tidal_strain_rate = prm.get_double ("Constant tidal strain rate");
          if (prm.get ("Custom distribution of tidal strain rate") == "constant")
            strain_rate_distribution = constant;
          else if (prm.get ("Custom distribution of tidal strain rate") == "latitudinal variation")
            strain_rate_distribution = latitudinal_variation;
          else
            AssertThrow (false, ExcMessage ("Not a valid custom distribution of tidal strain rate."));
          maximum_tidal_strain_rate = prm.get_double ("Maximum tidal strain rate");
          minimum_tidal_strain_rate = prm.get_double ("Minimum tidal strain rate");
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
                                  "The default equation ignores regional (radial/lateral) changes. "
                                  "This equation is the Eq.12 from Tobie et al. (2003) (https://doi.org/10.1029/2003JE002099). "
                                  "Selecting 'latitudinal variation' from 'Custom distribution of tidal strain rate' allows simplified latitudinal variation "
                                  "with cosine function, 'Maximum tidal strain rate' and 'Minimum tidal strain rate. "
                                  "Latitudinal variation of tidal strain rate is shown in Fig.3 from Nimmo et al. (2007) (https://doi.org/10.1016/j.icarus.2007.04.021). "
                                  "Unit: W/m^3.")
  }
}
