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
