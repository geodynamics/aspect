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

#ifndef _aspect_material_model_additional_outputs_anisotropic_viscosity_h
#define _aspect_material_model_additional_outputs_anisotropic_viscosity_h


#include <aspect/material_model/interface.h>

namespace aspect
{
  namespace MaterialModel
  {
    /**
      * Additional output fields for anisotropic viscosities to be added to
      * the MaterialModel::MaterialModelOutputs structure and filled in the
      * MaterialModel::Interface::evaluate() function.
      */
    template <int dim>
    class AnisotropicViscosity : public NamedAdditionalMaterialOutputs<dim>
    {
      public:
        AnisotropicViscosity(const unsigned int n_points);

        std::vector<double>
        get_nth_output(const unsigned int idx) const override;

        /**
         * Stress-strain "director" tensors at the given positions. This
         * variable is used to implement anisotropic viscosity.
         *
         * @note The strain rate term in equation (1) of the manual will be
         * multiplied by this tensor *and* the viscosity scalar ($\eta$), as
         * described in the manual section titled "Constitutive laws". This
         * variable is assigned the rank-four identity tensor by default.
         * This leaves the isotropic constitutive law unchanged if the material
         * model does not explicitly assign a value.
         */
        std::vector<SymmetricTensor<4,dim>> stress_strain_directors;

        static Tensor<2,3> euler_angles_to_rotation_matrix(double phi1, double theta, double phi2);
    };
  }
}

#endif
