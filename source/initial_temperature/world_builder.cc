/*
  Copyright (C) 2011 - 2017 by the authors of the ASPECT code.

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

#ifdef ASPECT_USE_WORLD_BUILDER
#include <aspect/initial_temperature/world_builder.h>
#include <world_builder/world.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/gravity_model/interface.h>


namespace aspect
{
  namespace InitialTemperature
  {
    template <int dim>
    WorldBuilder<dim>::WorldBuilder ()
    {}

    template <int dim>
    double
    WorldBuilder<dim>::
    initial_temperature (const Point<dim> &position) const
    {
      return this->get_world_builder().temperature(Utilities::convert_point_to_array(position),
                                                   this->get_geometry_model().depth(position),
                                                   this->get_gravity_model().gravity_vector(position).norm());
    }

  }
}

// explicit instantiations
namespace aspect
{
  namespace InitialTemperature
  {
    ASPECT_REGISTER_INITIAL_TEMPERATURE_MODEL(WorldBuilder,
                                              "world builder",
                                              "Specify the initial temperature in through the World Buider."
                                              "More information on the World Builder can be found at "
                                              "https://geodynamicworldbuilder.github.io."
                                              "Make sure to specify the location of the World Builder file "
                                              "in the parameter 'World builder file'.")
  }
}
#endif
