/*
  Copyright (C) 2011, 2012 by the authors of the ASPECT code.

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
/*  $Id$  */


#include <aspect/mesh_refinement/composition.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/numerics/error_estimator.h>

namespace aspect
{
  namespace MeshRefinement
  {
    template <int dim>
    void
    Composition<dim>::execute(Vector<float> &indicators) const
    {
      AssertThrow (this->n_compositional_fields() >= 1,
                   ExcMessage ("This refinement criterion can not be used when no "
                               "compositional fields are active!"));
      indicators = 0;

      for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
        {
          Vector<float> this_indicator (indicators.size());

          std::vector<bool> composition_component (dim+2+this->n_compositional_fields(), false);
          composition_component[dim+2+c] = true;
          KellyErrorEstimator<dim>::estimate (this->get_dof_handler(),
//TODO: Replace the 2 by something reasonable, adjusted to the polynomial degree
                                              QGauss<dim-1>(2),
                                              typename FunctionMap<dim>::type(),
                                              this->get_solution(),
                                              this_indicator,
                                              composition_component,
                                              0,
                                              0,
                                              this->get_triangulation().locally_owned_subdomain());
          indicators += this_indicator;
        }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MeshRefinement
  {
    ASPECT_REGISTER_MESH_REFINEMENT_CRITERION(Composition,
                                              "composition",
                                              "A mesh refinement criterion that computes "
                                              "refinement indicators from the compositional fields. "
                                              "If there is more than one compositional field, then "
                                              "it simply takes the sum of the indicators computed "
                                              "from each of the compositional field.")
  }
}
