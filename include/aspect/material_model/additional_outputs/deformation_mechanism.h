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

#ifndef _aspect_material_model_additional_outputs_deformation_mechanism_h
#define _aspect_material_model_additional_outputs_deformation_mechanism_h


#include <aspect/material_model/interface.h>

namespace aspect
{
  namespace MaterialModel
  {
    /**
     * The deformation mechanism selected while evaluating a rheology.
     * The numerical values are used by the corresponding visualization
     * postprocessor.
     */
    enum class DeformationMechanism : unsigned int
    {
      diffusion = 0,
      dislocation = 1,
      plastic_yielding = 2,
      peierls = 3,
      maximum_viscosity = 4,
      grain_boundary_sliding = 5,
      frank_kamenetskii = 6,
      uninitialized = numbers::invalid_unsigned_int
    };

    /**
     * Additional output containing the dominant deformation mechanism at
     * each material model evaluation point.
     */
    template <int dim>
    class DeformationMechanismOutputs : public NamedAdditionalMaterialOutputs<dim>
    {
      public:
        DeformationMechanismOutputs(const unsigned int n_points);

        std::vector<double>
        get_nth_output(const unsigned int idx) const override;

        /**
         * The dominant deformation mechanism at each evaluation point.
         */
        small_vector<DeformationMechanism> deformation_mechanisms;
    };
  }
}

#endif
