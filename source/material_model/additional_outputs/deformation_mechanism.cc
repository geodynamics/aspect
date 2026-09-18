/*
  Copyright (C) 2026 by the authors of the ASPECT code.

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

#include <aspect/material_model/additional_outputs/deformation_mechanism.h>

namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    DeformationMechanismOutputs<dim>::DeformationMechanismOutputs(const unsigned int n_points)
      :
      NamedAdditionalMaterialOutputs<dim>(std::vector<std::string>(1, "dominant_deformation_mechanism")),
      deformation_mechanisms(n_points, DeformationMechanism::uninitialized)
    {}



    template <int dim>
    std::vector<double>
    DeformationMechanismOutputs<dim>::get_nth_output(const unsigned int idx) const
    {
      (void)idx; // suppress warning in release mode
      AssertIndexRange (idx, 1);

      std::vector<double> output(deformation_mechanisms.size());
      for (unsigned int i=0; i<deformation_mechanisms.size(); ++i)
        output[i] = static_cast<unsigned int>(deformation_mechanisms[i]);

      return output;
    }
  }
}



// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
#define INSTANTIATE(dim) \
  template class DeformationMechanismOutputs<dim>;

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
