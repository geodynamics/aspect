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


#include <aspect/mesh_refinement/temperature.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/numerics/error_estimator.h>

namespace aspect
{
  namespace MeshRefinement
  {
    template <int dim>
    void
    Temperature<dim>::execute(Vector<float> &indicators) const
    {
      indicators = 0;

      QGauss<dim-1> quadrature (this->get_fe().base_element(this->introspection().base_elements.temperature).degree+1);

      KellyErrorEstimator<dim>::estimate (this->get_dof_handler(),
                                          quadrature,
                                          typename FunctionMap<dim>::type(),
                                          this->get_solution(),
                                          indicators,
                                          this->introspection().component_masks.temperature,
                                          0,
                                          0,
                                          this->get_triangulation().locally_owned_subdomain());
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MeshRefinement
  {
    ASPECT_REGISTER_MESH_REFINEMENT_CRITERION(Temperature,
                                              "temperature",
                                              "A mesh refinement criterion that computes "
                                              "refinement indicators from the temperature field.")
  }
}
