// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2023 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#include <aspect/postprocess/visualization/boundary_indicator.h>
#include <aspect/geometry_model/interface.h>



namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      BoundaryIndicator<dim>::
      BoundaryIndicator ()
        :
        CellDataVectorCreator<dim>("") // no physical units
      {}



      template <int dim>
      std::pair<std::string, std::unique_ptr<Vector<float>>>
      BoundaryIndicator<dim>::execute() const
      {
        std::pair<std::string, std::unique_ptr<Vector<float>>>
        return_value ("boundary_indicator",
                      std::make_unique<Vector<float>>(this->get_triangulation().n_active_cells()));

        // retrieve the largest used boundary indicator and add one.
        // this value will be set for internal cells
        const types::boundary_id largest_boundary_id_plus_one =
          *this->get_geometry_model().get_used_boundary_indicators().rbegin() + 1;

        // loop over all of the surface cells and extract boundary indicators
        for (const auto &cell : this->get_dof_handler().active_cell_iterators())
          if (cell->is_locally_owned())
            {
              if (cell->at_boundary())
                {
                  types::boundary_id boundary_id = largest_boundary_id_plus_one;
                  for (const unsigned int f : cell->face_indices())
                    {
                      if (cell->face(f)->at_boundary())
                        {
                          boundary_id
                            = cell->face(f)->boundary_id();
                          break;
                        }
                    }
                  (*return_value.second)(cell->active_cell_index()) = static_cast<float> (boundary_id);
                }
              // internal cells are all set to the same boundary indicator value
              // of the largest boundary indicator used for the current geometry plus one.
              else
                {
                  (*return_value.second)(cell->active_cell_index()) = largest_boundary_id_plus_one;
                }
            }

        return return_value;
      }

    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(BoundaryIndicator,
                                                  "boundary indicators",
                                                  "A visualization output object that generates output "
                                                  "about the used boundary indicators. In a loop over the active "
                                                  "cells, if a cell lies at a domain boundary, the boundary indicator "
                                                  "of the face along the boundary is requested. In case the cell "
                                                  "does not lie along any domain boundary, the cell is assigned the "
                                                  "value of the largest used boundary indicator plus one. "
                                                  "When a cell is situated in one of the corners of the domain, "
                                                  "multiple faces will have a boundary indicator. This postprocessor "
                                                  "returns the value of the first face along a boundary that is encountered "
                                                  "in a loop over all the faces."
                                                  "\n\n"
                                                  "Physical units: None.")
    }
  }
}
