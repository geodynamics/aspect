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
         * cookbook, the CPO Bingham Average particle property must be included and Rotation format
         * must be set to euler angles. The resulting arrays are:
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

        /**
         * These variables are used in the flow law to compute anisotropic viscosity and the default values are
         * provided according to Hansen et al. (2016) and Hirth and Kohlstedt (2004)
         */
        double grain_size;
        double stress_exponent;
        double activation_energy;
        double fluidity_constant;
        double grain_size_exponent;

        /**
         * The parameter determines which method to use to get the anisotropic tensor for viscosity.
         * When false: use a pseudo inverse to invert the anisotropic tensor for fluidity and inverts iteratively for the scalar viscosity
         * When true: use an analytical inversion based on the orthotropic symmetry invariants of the strain-rate (see rathmann et.al.2021)
         * The advantages of using the analytical inversion is that it converges to cases with high stress exponent,
         * and increases computational efficiency by avoiding an iterative inversion.
         */
        bool use_analytical_inversion;

        /**
         * The iteration for computing scalar viscosity is terminated when
         * 1) the relative change falls below the relative tolerance;
         * 2) the number of iterations exceeds the maximum number of iterations.
         */
        double relative_tolerance;
        unsigned int max_iteration;

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
         * This function recasts a dim 6 rank-2 SymmetricTensor into a Dealii::FullMatrix and discard the out of plane components.
         * With that we can convert it from kelvin notation back to rank-4 tensor of dim specified by the model.
         */
        SymmetricTensor<4,dim> kelvin_to_r4_tensor(const Tensor<2,6> &V) const;

        /**
         * This function computes the anisotropic tensor for viscosity in the cpo reference in kelvin notation.
         * When rotating into the model frame and converting back from kelvin notation it is a rank-4 tensor,
         * which relates all components of the strain-rate with all components of the stress.
         * In the isotropic case this tensor multiplied by a strain-rate returns the deviatoric strain-rate.
         */
        SymmetricTensor<2,6> viscosity_tensor_cpo_frame(const double F, const double G, const double H,
                                                        const double L, const double M, const double N ) const;

        /**
         * This function computes the orthotropic strain-rate invariant of a function
         * based on strain rate in the cpo reference frame and the hill coefficients.
         * In the isotropic case (F,G,H = 0.5 and L,M,N=1.5) this function returns
         * $sqrt{\dot{\epslion}{ij}^d\dot{\epsilon}{ij}^d}$, where $\dot{\epsilon}^d$
         * is the deviatoric strain-rate.
         */
        double orthotropic_strain_rate_invariant( const Tensor<2,3> &strain_rate_cpo_frame,
                                                  const double F, const double G, const double H,
                                                  const double L, const double M, const double N) const;

    };
  }
}

#endif
