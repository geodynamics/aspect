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


#include <aspect/initial_composition/world_builder.h>
#include <aspect/geometry_model/interface.h>
#include <world_builder/world.h>


namespace aspect
{
  namespace InitialComposition
  {

    template <int dim>
    double
    WorldBuilder<dim>::
    initial_composition (const Point<dim> &position, const unsigned int n_comp) const
    {
      return this->get_world_builder().composition(Utilities::convert_point_to_array(position),
                                                   this->get_geometry_model().depth(position),
                                                   n_comp);
    }


    template <int dim>
    void
    WorldBuilder<dim>::declare_parameters (ParameterHandler &/*prm*/)
    {}


    template <int dim>
    void
    WorldBuilder<dim>::parse_parameters (ParameterHandler &/*prm*/)
    {}

  }
}

// explicit instantiations
namespace aspect
{
  namespace InitialComposition
  {
    ASPECT_REGISTER_INITIAL_COMPOSITION_MODEL(WorldBuilder,
                                              "world builder",
                                              "Specify the composition in terms of an explicit formula. The format of these "
                                              "functions follows the syntax understood by the "
                                              "muparser library, see Section~\\ref{sec:muparser-format}.")
  }
}
