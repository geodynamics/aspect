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


#include <aspect/mesh_refinement/topography.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/fe/fe_values.h>

namespace aspect
{
  namespace MeshRefinement
  {
    template <int dim>
    void
    Topography<dim>::execute(Vector<float> &indicators) const
    {
      // For calculating surface topography in the respective postprocessor
      // we use the pressure in the middle of the cell. Thus, to get an
      // accurate result, all the cells at the upper boundary should have
      // the same level of refinement. This postprocessor causes a refinement
      // in the upermost cells, which also makes sure that the upper
      // boundary layer is resolved as good as possible.

      indicators = 0;

      // need quadrature formula to calculate the depth
      // evaluate a single point per cell
      const QMidpoint<dim> quadrature_formula;

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_quadrature_points |
                               update_JxW_values);

      // iterate over all of the cells and choose the ones at the upper
      // boundary for refinement (assign the largest error to them)
      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();
      unsigned int i=0;
      for (; cell!=endc; ++cell, ++i)
        if (cell->is_locally_owned())
          {
            fe_values.reinit (cell);
            const double depth = this->get_geometry_model().depth(fe_values.quadrature_point(0));
            if (cell->at_boundary() && depth < cell->diameter())
              indicators(i) = 1.0;
          }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MeshRefinement
  {
    ASPECT_REGISTER_MESH_REFINEMENT_CRITERION(Topography,
                                              "topography",
                                              "A class that implements a mesh refinement criterion, which "
                                              "always flags all cells in the uppermost layer for refinement. "
                                              "This is useful to provide high accuracy for processes at or "
                                              "close to the surface."
                                              "\n\n"
                                              "To use this refinement criterion, you may want to combine "
                                              "it with other refinement criteria, setting the 'Normalize "
                                              "individual refinement criteria' flag and using the 'max' "
                                              "setting for 'Refinement criteria merge operation'.")
  }
}
