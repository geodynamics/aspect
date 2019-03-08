/*
  Copyright (C) 2011 - 2018 by the authors of the ASPECT code.

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
#include <aspect/boundary_composition/initial_composition.h>
#include <aspect/boundary_composition/box.h>
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
      if (this->get_boundary_composition_manager().template has_matching_boundary_composition_model<BoundaryComposition::InitialComposition<dim> >())
        {
          const BoundaryComposition::InitialComposition<dim> &initial_composition =
            this->get_boundary_composition_manager().template get_matching_boundary_composition_model<BoundaryComposition::InitialComposition<dim> >();

          std::cout << "InitialComposition is Box:"
                    << Plugins::plugin_type_matches<const BoundaryComposition::Box<dim> >(initial_composition)
                    << std::endl;

          Plugins::get_plugin_as_type<const BoundaryComposition::Box<dim> >(initial_composition);
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

