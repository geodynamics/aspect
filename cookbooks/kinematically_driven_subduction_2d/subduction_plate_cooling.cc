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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include "subduction_plate_cooling.h"
#include <aspect/geometry_model/box.h>


namespace aspect
{
  namespace InitialTemperature
  {
    template <int dim>
    double
    SubductionPlateCooling<dim>::
    initial_temperature (const Point<dim> &position) const
    {
      // This initial condition only makes sense if the geometry is a
      // box. Throw an error if that is not the case.
      AssertThrow(Plugins::plugin_type_matches<const GeometryModel::Box<dim>>(this->get_geometry_model()),
                  ExcMessage("Subduction plate cooling initial temperature is not implemented "
                             "for any other geometry model than a box. ") );

      // Model domain is 3000 x 670 km.
      // Geometry constants:
      // surface temperature
      const double Ts=273.15;
      // temperature at the base of the plates
      const double Tm = 1573.15;
      // the base of the subducting plate at 110 km depth
      const double y70 = -110000.0;
      // the base of the overriding plate at 82 km depth
      const double y40 = -82000.0;
      // thermal diffusivity
      const double kappa = 1.0e-6;
      // the age of the subducting plate
      const double t70 = 70e6 * year_in_seconds;
      // the age of the overriding plate
      const double t40 = 40e6 * year_in_seconds;

      // the temperature to return
      double temperature = 0.;

      // Calculate the temperature in the portions of the model domain.
      // overriding plate + subduction interface
      if (position[0] < 1500000.0)
        {
          // lithosphere
          if (position[1] >= 588000.0)
            temperature=Ts+(Tm-Ts)*( ((position[1]-670000.0)/y40) + (2.0/numbers::PI)*(std::exp((-kappa*numbers::PI*numbers::PI*t40)/(y40*y40))*std::sin(numbers::PI*(position[1]-670000.0)/y40) +
                                                                                       0.5*std::exp((-kappa*2.0*2.0*numbers::PI*numbers::PI*t40)/(y40*y40))*std::sin(2.0*numbers::PI*(position[1]-670000.0)/y40)));
          // sublithospheric mantle
          else
            temperature=-0.25*(position[1]/1000.0)+1720.15;
        }
      // subducting plate
      else
        {
          // lithosphere
          if (position[1] >= 560000.0)
            temperature=Ts+(Tm-Ts)*( ((position[1]-670000.0)/y70) + (2.0/numbers::PI)*(std::exp((-kappa*numbers::PI*numbers::PI*t70)/(y70*y70))*std::sin(numbers::PI*(position[1]-670000.0)/y70) +
                                                                                       0.5*std::exp((-kappa*2.0*2.0*numbers::PI*numbers::PI*t70)/(y70*y70))*std::sin(2.0*numbers::PI*(position[1]-670000.0)/y70)));
          // sublithospheric mantle
          else
            temperature=-0.25*(position[1]/1000.0)+1713.15;
        }
      return temperature;
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace InitialTemperature
  {
    ASPECT_REGISTER_INITIAL_TEMPERATURE_MODEL(SubductionPlateCooling,
                                              "subduction plate cooling",
                                              "An initial temperature field as specified in Quinquis (2014) "
                                              "using the plate cooling model to prescribe temperatures "
                                              "for an overriding and a subducting plate. ")
  }
}
