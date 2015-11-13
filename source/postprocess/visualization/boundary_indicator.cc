/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include <aspect/postprocess/visualization/boundary_indicator.h>
#include <aspect/simulator_access.h>




namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {

      template <int dim>
      std::pair<std::string, Vector<float> *>
      BoundaryIndicator<dim>::execute() const
      {
        std::pair<std::string, Vector<float> *>
        return_value ("boundary_indicator",
                      new Vector<float>(this->get_triangulation().n_active_cells()));

        // retrieve the largest used boundary indicator and add one.
        // this value will be set for internal cells
        const types::boundary_id largest_boundary_id_plus_one =
          *this->get_geometry_model().get_used_boundary_indicators().rbegin() + 1;

        // loop over all of the surface cells and extract boundary indicators
        typename DoFHandler<dim>::active_cell_iterator
        cell = this->get_dof_handler().begin_active(),
        endc = this->get_dof_handler().end();

        unsigned int cell_index = 0;
        for (; cell!=endc; ++cell,++cell_index)
          if (cell->is_locally_owned())
            {
              if (cell->at_boundary())
                {
                  types::boundary_id boundary_id = largest_boundary_id_plus_one;
                  for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
                    {
                      if (cell->face(f)->at_boundary())
                        {
                          boundary_id
#if DEAL_II_VERSION_GTE(8,3,0)
                            = cell->face(f)->boundary_id();
#else
                            = cell->face(f)->boundary_indicator();
#endif
                          break;
                        }
                    }
                  (*return_value.second)(cell_index) = static_cast<float> (boundary_id);
                }
              // internal cells are all set to the same boundary indicator value
              // of the largest boundary indicator used for the current geometry plus one.
              else
                {
                  (*return_value.second)(cell_index) = largest_boundary_id_plus_one;
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
                                                  "in a loop over all the faces. ")
    }
  }
}
