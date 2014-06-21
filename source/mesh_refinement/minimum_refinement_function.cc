/*
  Copyright (C) 2011 - 2014 by the authors of the ASPECT code.

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



#include <aspect/mesh_refinement/minimum_refinement_function.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>
#include <math.h>

namespace aspect
{
  namespace MeshRefinement
  {
    template <int dim>
    void
    MinimumRefinementFunction<dim>::tag_additional_cells () const
    {
      for (typename Triangulation<dim>::active_cell_iterator
           cell = this->get_triangulation().begin_active();
           cell != this->get_triangulation().end(); ++cell)
        {
          if (cell->is_locally_owned())
            {
              bool refine = false;
              bool clear_coarsen = false;

              for ( unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell;  ++v)
                {
                  const Point<dim> vertex = cell->vertex(v);

                  // TODO: This should be an input parameter for the user to decide
                  // whether to use depth or coordinates
                  const double depth = this->get_geometry_model().depth(vertex);
                  const Point<1> point(depth);

                  if (cell->level() <= rint(min_refinement_level.value(point)))
                    clear_coarsen = true;
                  if (cell->level() <  rint(min_refinement_level.value(point)))
                    {
                      refine = true;
                      break;
                    }
                }

              if (clear_coarsen)
                cell->clear_coarsen_flag ();
              if (refine)
                cell->set_refine_flag ();
            }
        }
    }

    template <int dim>
    void
    MinimumRefinementFunction<dim>::
    declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Mesh refinement");
      {

        prm.enter_subsection("Minimum refinement function");
        {
          /**
           * Let the function that describes the minimal level of refinement
           * as a function of depth declare its parameters.
           * This defines the minimum refinement level each cell should have,
           * and that can not be exceeded by coarsening.
           */
          Functions::ParsedFunction<1>::declare_parameters (prm, 1);
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    MinimumRefinementFunction<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Mesh refinement");
      {
        prm.enter_subsection("Minimum refinement function");
        try
          {
            min_refinement_level.parse_parameters (prm);
          }
        catch (...)
          {
            std::cerr << "ERROR: FunctionParser failed to parse\n"
                      << "\t'Mesh refinement.Minimum refinement function'\n"
                      << "with expression\n"
                      << "\t'" << prm.get("Function expression") << "'";
            throw;
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
    ASPECT_REGISTER_MESH_REFINEMENT_CRITERION(MinimumRefinementFunction,
                                              "minimum refinement function",
                                              "A mesh refinement criterion that ensures a "
                                              "minimum refinement level described by an "
                                              "explicit formula with the depth as argument. "
                                              "After reading in the function, its values are "
                                              "rounded to the nearest integer. ")
  }
}
