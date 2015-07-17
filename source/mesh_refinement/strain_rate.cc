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

      std::vector<SymmetricTensor<2,dim> > strain_rates (quadrature.size());

      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();
      unsigned int j=0;
      for (; cell!=endc; ++cell, ++j)
        if (cell->is_locally_owned())
          {
            fe_values.reinit(cell);

            fe_values[this->introspection().extractors.velocities].get_function_symmetric_gradients (this->get_solution(),
                strain_rates);

            indicators(j) = strain_rates[0].norm();
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
                                              "A mesh refinement criterion that computes the"
                                              "refinement indicators equal to the strain rate "
                                              "norm computed at the center of the elements.")
  }
}
