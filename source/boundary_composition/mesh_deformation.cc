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
                  ExcMessage ("The boundary composition plugin ``mesh deformation'' can only be used when a mesh deformation plugin is active."));
    }



    template <int dim>
    double
    MeshDeformation<dim>::
    boundary_composition (const types::boundary_id boundary_indicator,
                          const Point<dim> &position,
                          const unsigned int compositional_field) const
    {
      return this->get_mesh_deformation_handler().boundary_composition(boundary_indicator, position, compositional_field);
    }



    template <int dim>
    void
    MeshDeformation<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary composition model");
      {
        prm.enter_subsection("Mesh deformation");
        {
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }


    template <int dim>
    void
    MeshDeformation<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary composition model");
      {
        prm.enter_subsection("Mesh deformation");
        {
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
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
                                               "is chosen to be the same as given in the initial "
                                               "conditions."
                                               "\n\n"
                                               "Because this class simply takes what the initial "
                                               "composition had described, this class can not "
                                               "know certain pieces of information such as the "
                                               "minimal and maximal composition on the boundary. "
                                               "For operations that require this, for example in "
                                               "post-processing, this boundary composition model "
                                               "must therefore be told what the minimal and "
                                               "maximal values on the boundary are. This is done "
                                               "using parameters set in section ``Boundary composition model/Initial composition''.")
  }
}
