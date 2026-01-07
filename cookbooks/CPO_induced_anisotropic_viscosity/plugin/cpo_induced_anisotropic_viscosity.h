/*
  Copyright (C) 2019 - 2020 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_CPO_AV_3D_h
#define _aspect_material_model_CPO_AV_3D_h

#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/material_model/simple.h>
#include <aspect/material_model/equation_of_state/interface.h>
#include <aspect/simulator/assemblers/interface.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    template <int dim>
    class CPO_AV_3D : public MaterialModel::Simple<dim>
    {
      public:
        void initialize() override;

        void evaluate (const MaterialModel::MaterialModelInputs<dim> &in,
                       MaterialModel::MaterialModelOutputs<dim> &out) const override;

        bool is_compressible () const override;

        static void declare_parameters (ParameterHandler &prm);

        void parse_parameters (ParameterHandler &prm) override;

        void create_additional_named_outputs(MaterialModel::MaterialModelOutputs<dim> &out) const override;

      private:

        /**
         * Reference viscosity.
         */
        double eta;

        /**
         * Defining a minimum strain rate stabilizes the viscosity calculation,
         * which involves a division by the strain rate. Units: 1/s.
         */
        double min_strain_rate;

        /**
         * These are arrays that store eigenvalues of the olivine textures in a, b, and c axis, and
         * the olivine texture represented by Euler angles. For more details, please refer to the
         * cpo_bingham_average particle property. To use the anisotropic viscosity plugin in this
         * cookbook, the CPO Bingham Average particle property must be included and Use rotation matrix
         * must be set to false. The resulting arrays are:
         * cpo_bingham_avg_a = [phi, eigenvalue 1 for a-axis, eigenvalue 2 for a-axis, eigenvalue 3 for a-axis]
         * cpo_bingham_avg_b = [phi, eigenvalue 1 for b-axis, eigenvalue 2 for b-axis, eigenvalue 3 for b-axis]
         * cpo_bingham_avg_c = [phi, eigenvalue 1 for c-axis, eigenvalue 2 for c-axis, eigenvalue 3 for c-axis]
         * They are used in computing rotation matrix with regards to the CPO reference frame, and
         * the anisotropic Hill coefficients FGHLMN.
         */
        std::vector<double> cpo_bingham_avg_a, cpo_bingham_avg_b, cpo_bingham_avg_c;

        /**
         * These are arrays that store coefficients used to compute the anisotropic Hill coefficients FGHLMN from
         * a certain olivine texture represented with the eigenvalues of its a-, b-, and c-axis. Each array contains
         * 9 coefficients and 1 constant.
         */
        std::vector<double> CnI_F, CnI_G, CnI_H, CnI_L, CnI_M, CnI_N;

        double grain_size;

        EquationOfState::LinearizedIncompressible<dim> equation_of_state;

        void set_assemblers(const SimulatorAccess<dim> &,
                            Assemblers::Manager<dim> &assemblers) const;

        /**
         * This function computes the Moore-Penrose pseudoinverse of a matrix A
         * using Singular Value Decomposition (SVD). It takes a LAPACKFullMatrix A
         * to be inverted and outputs its pseudoinverse A_pinv.
         * SVD Method:
         *   A = U * Sigma * V^T
         *   A_pinv = V * Sigma_pinv * U^T
         * Singular values smaller than a fixed tolerance (1e-12) are treated as zero
         * for numerical stability.
         */
        void pseudoinverse(LAPACKFullMatrix<double> &A,
                           LAPACKFullMatrix<double> &A_pinv) const;

        /**
         * This conversion from euler angles to rotation matrix is different from the function with the same name
         * defined in utilities.cc. Both define the conversion with R3*R2*R1, while our R3 and R1 are the transpose
         * of those defined in utilities.cc. This change is made to fit our negative euler angle values.
         */
        Tensor<2,3> euler_angles_to_rotation_matrix(const double phi1,
                                                    const double theta,
                                                    const double phi2) const;

    };
  }
}

#endif
