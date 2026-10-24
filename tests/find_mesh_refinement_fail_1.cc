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
#include <aspect/mesh_refinement/boundary.h>


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
      this->get_mesh_refinement_manager().template get_matching_mesh_refinement_strategy<MeshRefinement::Boundary<dim>>();
      std::cout << "Boundary is found!" << std::endl;
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
