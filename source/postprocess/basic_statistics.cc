/*
  Copyright (C) 2011 - 2022 by the authors of the ASPECT code.

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


#include <aspect/postprocess/basic_statistics.h>
#include <aspect/material_model/simple.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/boundary_temperature/interface.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/global.h>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    BasicStatistics<dim>::execute (TableHandler &)
    {
      if (!Plugins::plugin_type_matches<const MaterialModel::Simple<dim>>(this->get_material_model()))
        return std::make_pair (std::string(),std::string());

      if ((this->get_timestep_number() == 0) &&
          (this->get_pre_refinement_step() == std::numeric_limits<unsigned int>::max()))
        {
          // Define the representative point (in this case a point with depth 0 -> at the surface)
          const double reference_depth = 0.0;
          const Point<dim> representative_point = this->get_geometry_model().representative_point(reference_depth);
          const double gravity = this->get_gravity_model().gravity_vector(representative_point).norm();
          const double model_depth = this->get_geometry_model().maximal_depth();

          // temperature contrast is only meaningful if boundary temperatures are prescribed, otherwise it is 0
          const double temperature_contrast = (this->has_boundary_temperature()
                                               ?
                                               this->get_boundary_temperature_manager().maximal_temperature(this->get_fixed_temperature_boundary_indicators())
                                               - this->get_boundary_temperature_manager().minimal_temperature(this->get_fixed_temperature_boundary_indicators())
                                               :
                                               0);

          // Compute material properties at the reference point using the adiabatic
          // profile as reference temperature and pressure
          const double temperature = this->get_adiabatic_conditions().temperature(representative_point);
          const double pressure = this->get_adiabatic_conditions().pressure(representative_point);

          MaterialModel::MaterialModelInputs<dim> in (1, this->n_compositional_fields());
          MaterialModel::MaterialModelOutputs<dim> out (1, this->n_compositional_fields());
          in.requested_properties = MaterialModel::MaterialProperties::thermal_conductivity |
                                    MaterialModel::MaterialProperties::equation_of_state_properties |
                                    MaterialModel::MaterialProperties::viscosity;

          in.position[0] = representative_point;
          in.temperature[0] = temperature;
          in.pressure[0] = pressure;
          for (unsigned int c = 0; c < this->n_compositional_fields(); ++c)
            in.composition[0][c] = 0.0;

          this->get_material_model().evaluate(in, out);

          const double thermal_diffusivity = ( (this->get_parameters().formulation_temperature_equation ==
                                                Parameters<dim>::Formulation::TemperatureEquation::reference_density_profile)
                                               ?
                                               out.thermal_conductivities[0] /
                                               (this->get_adiabatic_conditions().density(in.position[0]) * out.specific_heat[0])
                                               :
                                               out.thermal_conductivities[0] / (out.densities[0] * out.specific_heat[0])
                                             );

          // check whether diffusivity is set to 0 (in case of backward advection)
          const double Ra = (thermal_diffusivity != 0) ?
                            out.densities[0] *
                            gravity *
                            out.thermal_expansion_coefficients[0] *
                            temperature_contrast * std::pow(model_depth,3)/
                            (thermal_diffusivity * out.viscosities[0])
                            :
                            std::numeric_limits<double>::infinity();

          this->get_pcout() <<  std::endl;
          this->get_pcout() << "     Model domain depth (m):                        "
                            << model_depth
                            << std::endl;
          this->get_pcout() << "     Temperature contrast across model domain (K):  "
                            << temperature_contrast
                            << std::endl;
          this->get_pcout() << "     Reference depth (m):                           "
                            << reference_depth
                            << std::endl;
          this->get_pcout() << "     Reference temperature (K):                     "
                            << temperature
                            << std::endl;
          this->get_pcout() << "     Reference pressure (Pa):                       "
                            << pressure
                            << std::endl;
          this->get_pcout() << "     Reference gravity (m/s^2):                     "
                            << gravity
                            << std::endl;
          this->get_pcout() << "     Reference density (kg/m^3):                    "
                            << out.densities[0]
                            << std::endl;
          this->get_pcout() << "     Reference thermal expansion coefficient (1/K): "
                            << out.thermal_expansion_coefficients[0]
                            << std::endl;
          this->get_pcout() << "     Reference specific heat capacity (J/(K*kg)):   "
                            << out.specific_heat[0]
                            << std::endl;
          this->get_pcout() << "     Reference thermal conductivity (W/(m*K)):      "
                            << out.thermal_conductivities[0]
                            << std::endl;
          this->get_pcout() << "     Reference viscosity (Pa*s):                    "
                            << out.viscosities[0]
                            << std::endl;
          this->get_pcout() << "     Reference thermal diffusivity (m^2/s):         "
                            << thermal_diffusivity
                            << std::endl;
          this->get_pcout() << "     Rayleigh number:                               "
                            << Ra
                            << std::endl;

          this->get_pcout() << std::endl;
        }
      return std::make_pair (std::string(),std::string());
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(BasicStatistics,
                                  "basic statistics",
                                  "A postprocessor that outputs some simplified statistics "
                                  "like the Rayleigh number and other quantities that only "
                                  "make sense in certain model setups. The output is written "
                                  "after completing initial adaptive refinement steps. "
                                  "The postprocessor assumes a point at the surface at the adiabatic "
                                  "surface temperature and pressure is a reasonable reference condition "
                                  "for computing these properties. Furthermore, the Rayleigh "
                                  "number is computed using the model depth (i.e. not the "
                                  "radius of the Earth), as we need a definition that is geometry "
                                  "independent. Take care when comparing these values to published "
                                  "studies and make sure they use the same definitions.")
  }
}
