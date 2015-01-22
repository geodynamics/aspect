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


#include <aspect/boundary_temperature/box.h>
#include <aspect/geometry_model/box.h>
#include <aspect/postprocess/heat_flux_statistics.h>
#include <aspect/postprocess/pressure_statistics.h>
#include <aspect/simulator_access.h>
#include <aspect/simulator.h>

#include <utility>
#include <limits>


namespace aspect
{
  namespace BoundaryTemperature
  {
    template <int dim>
    class Box2 : public Box<dim>, public aspect::SimulatorAccess<dim>
    {
        public:
        virtual void update();
    };
    template <int dim>
    void Box2<dim>::update()
    {

        const Postprocess::HeatFluxStatistics<dim> * heat_flux_postprocess=
        this->template find_postprocessor<const Postprocess::HeatFluxStatistics<dim> >();
        const Postprocess::PressureStatistics<dim> * pressure_postprocess=
            this->template find_postprocessor<const Postprocess::PressureStatistics<dim> >();
        if(pressure_postprocess==NULL)
	  std::cout << "PressureStatistics is not found!" << std::endl;
        else
	  std::cout << "PressureStatistics is found!" << std::endl;
        if(heat_flux_postprocess==NULL)
	  std::cout << "HeatFluxStatistics is not found!" << std::endl;
        else
	  std::cout << "HeatFluxStatistics is found!" << std::endl;         
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace BoundaryTemperature
  {
    ASPECT_REGISTER_BOUNDARY_TEMPERATURE_MODEL(Box2,
                                               "box2",
                                               "A model in which the temperature is chosen constant on "
                                               "all the sides of a box. For test")
  }
}
  
