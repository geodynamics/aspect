/*
  Copyright (C) 2011 - 2021 by the authors of the ASPECT code.

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


#include <aspect/mesh_refinement/boundary.h>
#include <aspect/geometry_model/interface.h>

namespace aspect
{
  namespace MeshRefinement
  {
    template <int dim>
    void
    Boundary<dim>::execute(Vector<float> &indicators) const
    {
      indicators = 0;

      // iterate over all of the cells and choose the ones at the indicated
      // boundaries for refinement (assign the largest error to them)

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned() && cell->at_boundary())
          for (const unsigned int face_no : cell->face_indices())
            if (cell->face(face_no)->at_boundary())
              {
                const types::boundary_id boundary_indicator
                  = cell->face(face_no)->boundary_id();
                if ( boundary_refinement_indicators.find(boundary_indicator) !=
                     boundary_refinement_indicators.end() )
                  {
                    indicators(cell->active_cell_index()) = 1.0;
                    break;  // no need to loop over the rest of the faces
                  }
              }

    }

    template <int dim>
    void
    Boundary<dim>::
    declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Mesh refinement");
      {
        prm.enter_subsection("Boundary");
        {
          prm.declare_entry ("Boundary refinement indicators", "",
                             Patterns::List (Patterns::Anything()),
                             "A comma separated list of names denoting those boundaries "
                             "where there should be mesh refinement."
                             "\n\n"
                             "The names of the boundaries listed here can either be "
                             "numbers (in which case they correspond to the numerical "
                             "boundary indicators assigned by the geometry object), or they "
                             "can correspond to any of the symbolic names the geometry object "
                             "may have provided for each part of the boundary. You may want "
                             "to compare this with the documentation of the geometry model you "
                             "use in your model.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    Boundary<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Mesh refinement");
      {
        prm.enter_subsection("Boundary");
        {
          const std::vector<types::boundary_id> x_boundary_refinement_indicators
            = this->get_geometry_model().translate_symbolic_boundary_names_to_ids(Utilities::split_string_list
                                                                                  (prm.get ("Boundary refinement indicators")));
          boundary_refinement_indicators
            = std::set<types::boundary_id> (x_boundary_refinement_indicators.begin(),
                                            x_boundary_refinement_indicators.end());
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

  }
}

// explicit instantiations
namespace aspect
{
  namespace MeshRefinement
  {
    ASPECT_REGISTER_MESH_REFINEMENT_CRITERION(Boundary,
                                              "boundary",
                                              "A class that implements a mesh refinement criterion which "
                                              "always flags all cells on specified boundaries for refinement. "
                                              "This is useful to provide high accuracy for processes at or "
                                              "close to the edge of the model domain."
                                              "\n\n"
                                              "To use this refinement criterion, you may want to combine "
                                              "it with other refinement criteria, setting the 'Normalize "
                                              "individual refinement criteria' flag and using the `max' "
                                              "setting for 'Refinement criteria merge operation'.")
  }
}
