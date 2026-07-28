/*
  Copyright (C) 2025 - by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_additional_outputs_viscosity_without_elasticity_h
#define _aspect_material_model_additional_outputs_viscosity_without_elasticity_h

#include <aspect/material_model/visco_plastic.h>


namespace aspect
{
  namespace MaterialModel
  {
    /**
     * Additional output fields for the viscosity of viscous rheology
     * to be added to the MaterialModel::MaterialModelOutputs structure
     * and filled in the MaterialModel::evaluate() function.
     *
     * This structure allows to use viscosity of viscous rheology in tidal heating
     * without the elastic contributions to the effective viscosity.
     * This is because tidal heating is a form of viscoelastic shear dissipation.
     * when elasticity is enabled.
     */
    template <int dim>
    class ViscosityWithoutElasticityAdditionalOutputs : public MaterialModel::NamedAdditionalMaterialOutputs<dim>
    {
      public:
        ViscosityWithoutElasticityAdditionalOutputs(const unsigned int n_points);

        std::vector<double> get_nth_output(const unsigned int idx) const override;

        /**
         * The viscosity of viscous rheology, when elasticity is enabled.
         */
        std::vector<double> viscosity_without_elasticity;
    };
  }
}

#endif
