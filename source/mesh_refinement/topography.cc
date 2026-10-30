// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2019 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#include <aspect/mesh_refinement/topography.h>
#include <aspect/geometry_model/interface.h>

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
      // in the uppermost cells, which also makes sure that the upper
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
      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          {
            fe_values.reinit (cell);
            const double depth = this->get_geometry_model().depth(fe_values.quadrature_point(0));
            if (cell->at_boundary() && depth < cell->diameter())
              indicators(cell->active_cell_index()) = 1.0;
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
                                              "individual refinement criteria' flag and using the `max' "
                                              "setting for 'Refinement criteria merge operation'.")
  }
}
