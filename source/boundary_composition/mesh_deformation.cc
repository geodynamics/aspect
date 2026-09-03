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


#include <aspect/boundary_composition/mesh_deformation.h>
#include <aspect/initial_composition/interface.h>


namespace aspect
{
  namespace BoundaryComposition
  {
    // ------------------------------ MeshDeformation -------------------

    template <int dim>
    void
    MeshDeformation<dim>::initialize()
    {
      // Check that mesh deformation is actually active. The mesh deformation plugins
      // are not initialized yet, so we can only check the input parameters.
      AssertThrow(this->get_parameters().mesh_deformation_enabled == true,
                  ExcMessage(
                    "The boundary composition plugin ``mesh deformation'' can only be used when a mesh deformation plugin is active."));
    }



    template <int dim>
    double
    MeshDeformation<dim>::boundary_composition(const types::boundary_id boundary_indicator,
                                               const Point<dim>        &position,
                                               const unsigned int       compositional_field) const
    {
      return this->get_mesh_deformation_handler().boundary_composition(boundary_indicator, position, compositional_field);
    }

  }
}

// explicit instantiations
namespace aspect
{
  namespace BoundaryComposition
  {
    ASPECT_REGISTER_BOUNDARY_COMPOSITION_MODEL(MeshDeformation,
                                               "mesh deformation",
                                               "A model in which the composition at the boundary "
                                               "is retrieved from the active mesh deformation plugins."
                                               "\n\n"
                                               "The active mesh deformation plugins can each return "
                                               "a value for the compositional fields that are prescribed on "
                                               "a boundary; their values are summed per field. If a mesh "
                                               "deformation plugin does not implement the boundary composition "
                                               "function, a default value of zero is returned.")
  }
}
