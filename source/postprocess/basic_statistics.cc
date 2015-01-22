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


#include <aspect/postprocess/basic_statistics.h>
#include <aspect/material_model/simple.h>
#include <aspect/simulator_access.h>
#include <aspect/global.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    BasicStatistics<dim>::execute (TableHandler &statistics)
    {
      if (this->get_time() == 0e0)
        {
          if (dynamic_cast<const MaterialModel::Simple<dim> *>(&this->get_material_model()) != 0)
            {

              const MaterialModel::Simple<dim> &material_model
                = dynamic_cast<const MaterialModel::Simple<dim> &>(this->get_material_model());

              const double h = this->get_geometry_model().maximal_depth();

              // dT is only meaningful if boundary temperatures are prescribed, otherwise it is 0
              const double dT = (&this->get_boundary_temperature())
                                ?
                                this->get_boundary_temperature().maximal_temperature(this->get_fixed_temperature_boundary_indicators())
                                - this->get_boundary_temperature().minimal_temperature(this->get_fixed_temperature_boundary_indicators())
                                :
                                0;

              // we do not compute the compositions but give the functions below the value 0.0 instead
              std::vector<double> composition_values(this->n_compositional_fields(),0.0);

              Point<dim> representative_point = this->get_geometry_model().representative_point(0);
              const double gravity = this->get_gravity_model().gravity_vector(representative_point).norm();
              const double Ra = material_model.reference_density()*
                                gravity*
                                material_model.reference_thermal_expansion_coefficient()*
                                dT*std::pow(h,3)/
                                (material_model.reference_thermal_diffusivity()*
                                 material_model.reference_viscosity());

              this->get_pcout()<<  std::endl;
              this->get_pcout()<< "     Reference density (kg/m^3):                    "
                               << material_model.reference_density()
                               << std::endl;
              this->get_pcout()<< "     Reference gravity (m/s^2):                     "
                               << gravity
                               << std::endl;
              this->get_pcout()<< "     Reference thermal expansion (1/K):             "
                               << material_model.reference_thermal_expansion_coefficient()
                               << std::endl;
              this->get_pcout()<< "     Temperature contrast across model domain (K):  "
                               << dT
                               << std::endl;
              this->get_pcout()<< "     Model domain depth (m):                        "
                               << h
                               << std::endl;
              this->get_pcout()<< "     Reference thermal diffusivity (m^2/s):         "
                               << material_model.reference_thermal_diffusivity()
                               << std::endl;
              this->get_pcout()<< "     Reference viscosity (Pas):                     "
                               << material_model.reference_viscosity()
                               << std::endl;
              this->get_pcout()<< "     Ra number:                                     "
                               << Ra
                               << std::endl;
              this->get_pcout()<< "     k_value:                                       "
                               << material_model.thermal_conductivity(dT, dT, composition_values, representative_point) //TODO: dT for the pressure is wrong
                               << std::endl;
              this->get_pcout()<< "     reference_cp:                                  "
                               << material_model.reference_cp()
                               << std::endl;
              this->get_pcout()<< "     reference_thermal_diffusivity:                 "
                               << material_model.thermal_conductivity(dT, dT, composition_values, representative_point)/(material_model.reference_density()*material_model.reference_cp()) //TODO: dT for the pressure is wrong
                               << std::endl;
              this->get_pcout()<<  std::endl;
            }
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
                                  "A postprocessor that computes some simplified statistics "
                                  "like the Rayleigh number and other quantities that only "
                                  "make sense in certain model setups.")
  }
}
