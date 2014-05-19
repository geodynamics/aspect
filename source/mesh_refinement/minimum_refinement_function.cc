/*
  Copyright (C) 2014 by the authors of the ASPECT code.

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
/*  $Id: composition.cc 1880 2013-09-11 14:55:19Z dannberg $  */


#include <aspect/mesh_refinement/minimum_refinement_function.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

namespace aspect
{
  namespace MeshRefinement
  {
    template <int dim>
    void
    MinimumRefinementFunction<dim>::tag_additional_cells (unsigned int max_grid_level) const
    {
      // evaluate a single point per cell
      const QMidpoint<dim> quadrature_formula;
      const unsigned int n_q_points = quadrature_formula.size();

      FEValues<dim> fe_values (this->get_mapping(),
        		               this->get_fe(),
                               quadrature_formula,
                               update_values   |
                               update_quadrature_points );

      // ensure minimum refinement level
      const unsigned int max_level = std::min(max_grid_level+1,this->get_triangulation().n_levels());
      for (unsigned int level=0; level < max_level;++level)
        for (typename Triangulation<dim>::active_cell_iterator
             cell = this->get_triangulation().begin_active(level);
             cell != this->get_triangulation().end_active(level); ++cell)
        {
          if (cell->is_locally_owned())
          {
        	fe_values.reinit(cell);
        	const double depth = this->get_geometry_model().depth(fe_values.quadrature_point(0));
        	const Point<1> point(depth);
          	if (level <= std::rint(min_refinement_level.value(point)))
              cell->clear_coarsen_flag ();
          	if (level <  std::rint(min_refinement_level.value(point)))
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
          Functions::ParsedFunction<1>::declare_parameters (prm, 1);
          /*               "The minimum refinement level each cell should have, "
                           "and that can not be exceeded by coarsening. " */
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
        {
        	min_refinement_level.parse_parameters (prm);
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
                                              "explicit formula with the depth as argument.")
  }
}
