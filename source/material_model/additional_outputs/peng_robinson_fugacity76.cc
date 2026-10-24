/*
  Copyright (C) 2025 - by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file doc/COPYING. If not see
  <http://www.gnu.org/licenses/>.
*/

#include <aspect/material_model/additional_outputs/peng_robinson_fugacity76.h>

#include <deal.II/base/signaling_nan.h>

namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    PengRobinsonFugacity<dim>::PengRobinsonFugacity(
      const unsigned int n_points)
      :
      NamedAdditionalMaterialOutputs<dim>({"fugacity"}),
                                   fugacities(n_points, numbers::signaling_nan<double>())
    {}



    template <int dim>
    std::vector<double>
    PengRobinsonFugacity<dim>::get_nth_output(
      const unsigned int idx) const
    {
      AssertIndexRange(idx, 1);
      return fugacities;
    }
  }
}



// Explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
#define INSTANTIATE(dim) \
  template class PengRobinsonFugacity<dim>;

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
