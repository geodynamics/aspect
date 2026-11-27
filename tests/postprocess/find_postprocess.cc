// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#include <aspect/boundary_temperature/box.h>
#include <aspect/geometry_model/box.h>
#include <aspect/postprocess/heat_flux_statistics.h>
#include <aspect/postprocess/pressure_statistics.h>


#include <utility>
#include <limits>


namespace aspect
{
  namespace BoundaryTemperature
  {
    template <int dim>
    class Box2 : public Box<dim>
    {
      public:
        virtual void update();
    };
    template <int dim>
    void Box2<dim>::update()
    {
      if (this->get_postprocess_manager().template has_matching_postprocessor<Postprocess::PressureStatistics<dim>>())
        std::cout << "PressureStatistics is found!" << std::endl;
      else
        std::cout << "PressureStatistics is not found!" << std::endl;

      if (this->get_postprocess_manager().template has_matching_postprocessor<Postprocess::HeatFluxStatistics<dim>>())
        std::cout << "HeatFluxStatistics is found!" << std::endl;
      else
        std::cout << "HeatFluxStatistics is not found!" << std::endl;

      try
        {
          this->get_postprocess_manager().template get_matching_postprocessor<Postprocess::PressureStatistics<dim>>();
          std::cout << "PressureStatistics is found!" << std::endl;
        }
      catch (...)
        {
          std::cout << "PressureStatistics is not found!" << std::endl;
        }

      try
        {
          this->get_postprocess_manager().template get_matching_postprocessor<Postprocess::HeatFluxStatistics<dim>>();
          std::cout << "HeatFluxStatistics is found!" << std::endl;
        }
      catch (...)
        {
          std::cout << "HeatFluxStatistics is not found!" << std::endl;
        }
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
