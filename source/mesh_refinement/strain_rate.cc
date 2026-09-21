// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2021 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#include <aspect/mesh_refinement/strain_rate.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/derivative_approximation.h>

namespace aspect
{
  namespace MeshRefinement
  {
    template <int dim>
    void
    StrainRate<dim>::execute(Vector<float> &indicators) const
    {
      indicators = 0;

      const QMidpoint<dim> quadrature;

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature,
                               update_quadrature_points | update_values | update_gradients);

      std::vector<SymmetricTensor<2,dim>> strain_rates (quadrature.size());

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          {
            const unsigned int idx = cell->active_cell_index();
            fe_values.reinit(cell);

            fe_values[this->introspection().extractors.velocities].get_function_symmetric_gradients (this->get_solution(),
                strain_rates);

            indicators(idx) = strain_rates[0].norm();
          }

    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MeshRefinement
  {
    ASPECT_REGISTER_MESH_REFINEMENT_CRITERION(StrainRate,
                                              "strain rate",
                                              "A mesh refinement criterion that computes the "
                                              "refinement indicators equal to the strain rate "
                                              "norm computed at the center of the elements.")
  }
}
