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


#include <aspect/boundary_temperature/box.h>
#include <aspect/geometry_model/box.h>
#include <aspect/mesh_refinement/density.h>
#include <aspect/mesh_refinement/boundary.h>
#include <aspect/simulator.h>

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
      if (this->get_mesh_refinement_manager().template has_matching_mesh_refinement_strategy<MeshRefinement::Density<dim>>())
        std::cout << "Density is found!" << std::endl;
      else
        std::cout << "Density is not found!" << std::endl;

      if (this->get_mesh_refinement_manager().template has_matching_mesh_refinement_strategy<MeshRefinement::Boundary<dim>>())
        std::cout << "Boundary is found!" << std::endl;
      else
        std::cout << "Boundary is not found!" << std::endl;

      try
        {
          this->get_mesh_refinement_manager().template get_matching_mesh_refinement_strategy<MeshRefinement::Density<dim>>();
          std::cout << "Density is found!" << std::endl;
        }
      catch (...)
        {
          std::cout << "Density is not found!" << std::endl;
        }

      try
        {
          this->get_mesh_refinement_manager().template get_matching_mesh_refinement_strategy<MeshRefinement::Boundary<dim>>();
          std::cout << "Boundary is found!" << std::endl;
        }
      catch (...)
        {
          std::cout << "Boundary is not found!" << std::endl;
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
