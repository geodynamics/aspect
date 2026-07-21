
/*
 Copyright (C) 2015 - 2024 by the authors of the ASPECT code.

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

#include "cpo_induced_anisotropic_viscosity.h"
#include <aspect/material_model/additional_outputs/anisotropic_viscosity.h>
#include <aspect/simulator/assemblers/stokes.h>
#include <aspect/simulator/assemblers/stokes_anisotropic_viscosity.h>
#include <aspect/simulator_signals.h>

#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/physics/notation.h>

namespace aspect
{

  namespace MaterialModel
  {

    template <int dim>
    void
    CPO_AV_3D<dim>::set_assemblers(const SimulatorAccess<dim> &,
                                   Assemblers::Manager<dim> &assemblers) const
    {
      // Search for the regular Stokes assembler and preconditioner assembler
      // and replace them with the versions for anisotropic viscosity
      for (unsigned int i=0; i<assemblers.stokes_preconditioner.size(); ++i)
        {
          if (Plugins::plugin_type_matches<Assemblers::StokesPreconditioner<dim>>(*(assemblers.stokes_preconditioner[i])))
            assemblers.stokes_preconditioner[i] = std::make_unique<Assemblers::StokesPreconditionerAnisotropicViscosity<dim>> ();
        }

      for (unsigned int i=0; i<assemblers.stokes_system.size(); ++i)
        {
          if (Plugins::plugin_type_matches<Assemblers::StokesIncompressibleTerms<dim>>(*(assemblers.stokes_system[i])))
            assemblers.stokes_system[i] = std::make_unique<Assemblers::StokesIncompressibleTermsAnisotropicViscosity<dim>> ();
        }
    }



    template<int dim>
    void
    CPO_AV_3D<dim>::pseudoinverse(LAPACKFullMatrix<double> &A,
                                  LAPACKFullMatrix<double> &A_pinv) const
    {
      Assert(A.m() == A.n(),
             ExcMessage("Pseudoinverse is only implemented for square matrices."));

      // Get the number of matrix rows(=columns) m
      const unsigned int m = A.m();

      // Compute SVD: A = U * Sigma * V^T
      A.compute_svd();
      const double tol = 1e-12;
      Vector<double> Sigma_pinv(m);
      for (unsigned int i=0; i<m; ++i)
        {
          Sigma_pinv[i] = (std::abs(A.singular_value(i)) > tol ? 1.0/A.singular_value(i) : 0.0);
        }

      // A^+ = V * Sigma^+ * U^T
      const LAPACKFullMatrix<double> U = A.get_svd_u();
      const LAPACKFullMatrix<double> VT = A.get_svd_vt();
      LAPACKFullMatrix<double> UT(m,m);
      U.transpose(UT);
      VT.Tmmult(A_pinv, UT, Sigma_pinv);
    }



    template<int dim>
    Tensor<2,3>
    CPO_AV_3D<dim>::euler_angles_to_rotation_matrix(const double phi1,
                                                    const double theta,
                                                    const double phi2) const
    {
      Tensor<2,3> rot_matrix;

      rot_matrix[0][0] = std::cos(phi2)*std::cos(phi1) - std::cos(theta)*std::sin(phi1)*std::sin(phi2);
      rot_matrix[0][1] = -std::cos(phi2)*std::sin(phi1) - std::cos(theta)*std::cos(phi1)*std::sin(phi2);
      rot_matrix[0][2] = std::sin(phi2)*std::sin(theta);
      rot_matrix[1][0] = std::sin(phi2)*std::cos(phi1) + std::cos(theta)*std::sin(phi1)*std::cos(phi2);
      rot_matrix[1][1] = -std::sin(phi2)*std::sin(phi1) + std::cos(theta)*std::cos(phi1)*std::cos(phi2);
      rot_matrix[1][2] = -std::cos(phi2)*std::sin(theta);
      rot_matrix[2][0] = std::sin(theta)*std::sin(phi1);
      rot_matrix[2][1] = std::sin(theta)*std::cos(phi1);
      rot_matrix[2][2] = std::cos(theta);

      AssertThrow(rot_matrix[2][2] <= 1.0, ExcMessage("Invalid rotation matrix: cos(theta) > 1"));

      return rot_matrix;
    }



    template <int dim>
    void
    CPO_AV_3D<dim>::
    initialize()
    {
      this->get_signals().set_assemblers.connect (std::bind(&CPO_AV_3D<dim>::set_assemblers,
                                                            std::cref(*this),
                                                            std::placeholders::_1,
                                                            std::placeholders::_2));
      // AssertThrow(dim==3,
      //            ExcMessage("Olivine has 3 independent slip systems, allowing for deformation in 3 independent directions, hence these models only work in 3D"));

      cpo_bingham_avg_a.push_back (this->introspection().compositional_index_for_name("phi1"));
      cpo_bingham_avg_a.push_back (this->introspection().compositional_index_for_name("eigvalue_a1"));
      cpo_bingham_avg_a.push_back (this->introspection().compositional_index_for_name("eigvalue_a2"));
      cpo_bingham_avg_a.push_back (this->introspection().compositional_index_for_name("eigvalue_a3"));

      cpo_bingham_avg_b.push_back (this->introspection().compositional_index_for_name("theta"));
      cpo_bingham_avg_b.push_back (this->introspection().compositional_index_for_name("eigvalue_b1"));
      cpo_bingham_avg_b.push_back (this->introspection().compositional_index_for_name("eigvalue_b2"));
      cpo_bingham_avg_b.push_back (this->introspection().compositional_index_for_name("eigvalue_b3"));

      cpo_bingham_avg_c.push_back (this->introspection().compositional_index_for_name("phi2"));
      cpo_bingham_avg_c.push_back (this->introspection().compositional_index_for_name("eigvalue_c1"));
      cpo_bingham_avg_c.push_back (this->introspection().compositional_index_for_name("eigvalue_c2"));
      cpo_bingham_avg_c.push_back (this->introspection().compositional_index_for_name("eigvalue_c3"));
    }

    // template <>
    // void
    // CPO_AV_3D<2>::evaluate (const MaterialModel::MaterialModelInputs<2> &,
    //                         MaterialModel::MaterialModelOutputs<2> &) const
    // {
    //   Assert (false, ExcNotImplemented());
    // }

    template<int dim>
    SymmetricTensor<4,dim>
    CPO_AV_3D<dim>::kelvin_to_r4_tensor(const Tensor<2,6> V) const
    {
      // Converts rank 2 kelvin notation viscosity tensor to full rank 4 tensor necessary in the assembler
      // recast into FullMatrix type with shape 6x6 for 3d or 3x3 for 2d applications
      FullMatrix<double> V_mat((dim==3) ? 6 : 3 , (dim==3) ? 6 : 3);
      if (dim == 3)
        {
          for (unsigned int vi=0; vi<6; ++vi)
            {
              for (unsigned int vj=0; vj<6; ++vj)
                {
                  V_mat[vi][vj] = V[vi][vj];
                }
            }
        }
      else // (dim == 2) // discard out of plane components
        {
          for (unsigned int vi=0; vi<3; ++vi)
            {
              for (unsigned int vj=0; vj<3; ++vj)
                {
                  V_mat[vi][vj] = V[((vi == 2) ? 5 : vi)][((vj == 2) ? 5 : vj)];
                }
            }
        }

      SymmetricTensor<4,dim> V_r4;
      dealii::Physics::Notation::Kelvin::to_tensor(V_mat, V_r4);
      return V_r4;
    }


    template<int dim>
    SymmetricTensor<2,6>
    CPO_AV_3D<dim>::viscosity_tensor_cpo_frame( const double F, const double G, const double H,
                                                const double L, const double M, const double N) const
    {
      std::vector<int> ji = {1,2,0}; // tuple of indices shifted by one
      std::vector<int> ki = {2,0,1}; // tuple of indices shifted by two

      SymmetricTensor<2,6> visc_tensor;
      std::vector<double> Hi = {F, G, H, M, N, L};

      const double gam = 4*(Hi[2]*Hi[1] + Hi[0]*Hi[2] + Hi[0]*Hi[1]);

      for (unsigned int i=0; i < 3; ++i)
        {
          // upper left part
          visc_tensor[i][i] = 2/(3*gam)*(4*Hi[i] + Hi[ji[i]] + Hi[ki[i]]);
          visc_tensor[ji[i]][ki[i]] = 2/(3*gam)*(Hi[i] - 2*Hi[ji[i]] - 2*Hi[ki[i]]);
          // sym: visc_tensor[ki[i]][ji[i]] = visc_tensor[ji[i]][ki[i]];

          // lower right part
          visc_tensor[i+3][i+3] = 3/(2*Hi[i+3]);
        }
      return visc_tensor;
    }

    template<int dim>
    double
    CPO_AV_3D<dim>::orthotropic_strain_rate_invariant(const Tensor<2,3> strain_rate_cpo_frame,
                                                      const double F, const double G, const double H,
                                                      const double L, const double M, const double N) const
    {
      // based on rathmann 2021 (https://doi.org/10.1017/jog.2022.33)
      // formulated in terms of Hill coefficients

      std::vector<int> ji = {1,2,0}; // tuple of indices shifted by one
      std::vector<int> ki = {2,0,1}; // tuple of indices shifted by two

      std::vector<double> Hi = {F, G, H, L, M, N};
      const double gam = 4*(Hi[2]*Hi[1] + Hi[0]*Hi[2] + Hi[0]*Hi[1]);

      double anisotropic_invariant = 0.0;

      for (unsigned int i=0; i < 3; ++i)
        {
          anisotropic_invariant += ( 2/(3*gam)*Utilities::fixed_power<2>(strain_rate_cpo_frame[i][i])*(
                                       4*Hi[i]+Hi[ji[i]]+Hi[ki[i]])
                                     + 2/(3*gam)*2*strain_rate_cpo_frame[ji[i]][ji[i]]*strain_rate_cpo_frame[ki[i]][ki[i]]*(
                                       Hi[i] - 2*Hi[ji[i]] - 2*Hi[ki[i]])
                                     + 3/Hi[i+3]*Utilities::fixed_power<2>(strain_rate_cpo_frame[ji[i]][ki[i]])
                                   );
        }

      return std::pow(anisotropic_invariant, 0.5);
    }


    template <int dim>
    void
    CPO_AV_3D<dim>::evaluate (const MaterialModel::MaterialModelInputs<dim> &in,
                              MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      const std::shared_ptr<MaterialModel::AnisotropicViscosity<dim>> anisotropic_viscosity = out.template get_additional_output_object<MaterialModel::AnisotropicViscosity<dim>>();

      EquationOfStateOutputs<dim> eos_outputs (1);
      const unsigned int viscosity_field_index = this->introspection().compositional_index_for_name("scalar_viscosity");

      const double n = stress_exponent;

      for (unsigned int q=0; q<in.n_evaluation_points(); ++q)
        {
          equation_of_state.evaluate(in, q, eos_outputs);

          // Get parameters for effective viscosity
          out.densities[q] = eos_outputs.densities[0];
          out.viscosities[q] = eta;
          out.thermal_expansion_coefficients[q] = eos_outputs.thermal_expansion_coefficients[0];
          out.specific_heat[q] = eos_outputs.specific_heat_capacities[0];
          out.thermal_conductivities[q] = 1;
          out.compressibilities[q] = eos_outputs.compressibilities[0];
          out.entropy_derivative_pressure[q] = eos_outputs.entropy_derivative_pressure[0];
          out.entropy_derivative_temperature[q] = eos_outputs.entropy_derivative_temperature[0];

          // Calculate effective viscosity
          const SymmetricTensor<2,dim> strain_rate = in.strain_rate[q];

          // even in 2d we need 3d strain-rates and velocity gradient tensors. So we make them 3d by
          // adding an extra dimension which is zero.
          SymmetricTensor<2,3> strain_rate_3d;
          strain_rate_3d[0][0] = strain_rate[0][0];
          strain_rate_3d[0][1] = strain_rate[0][1];
          strain_rate_3d[1][1] = strain_rate[1][1];
          if (dim == 3)
            {
              strain_rate_3d[0][2] = strain_rate[0][2];
              strain_rate_3d[1][2] = strain_rate[1][2];
              strain_rate_3d[2][2] = strain_rate[2][2];
            }

          const SymmetricTensor<2,3> deviatoric_strain_rate
            = (this->get_material_model().is_compressible()
               ?
               strain_rate_3d - 1./3. * trace(strain_rate_3d) * unit_symmetric_tensor<3>()
               :
               strain_rate_3d);

          // Create constant value to use for AV
          // const double A_o = fluidity_constant*std::exp(-530000/(8.314*in.temperature[q]));
          // // const double n = 3.5; //n=3 for test against VPSC, n=3.5 for D-Rex in ASPECT
          // // The values of A_o and 0.73 were picked so that Gamma = 3.5322e-15[1/(s*Pa^n)] if T=1600K and d=1000 microns
          // const double Gamma = (A_o/(std::pow(grain_size,0.73)));

          const double A_o = fluidity_constant*std::exp(-activation_energy/(8.314*std::max(in.temperature[q],1.0e-10)));
          // 1.1e5*std::exp(-530000/(8.314*in.temperature[q]));
          // The values of A_o and 0.73 were picked so that Gamma = 3.5322e-15[1/(s*Pa^n)] if T=1600K and d=1000 microns
          const double Gamma = (A_o/(std::pow(grain_size, grain_size_exponent)));

          // The computation of the viscosity tensor is only necessary after the simulator has been initialized
          // and when the condition allows dislocation creep
          if  ((this->simulator_is_past_initialization()) && (this->get_timestep_number() > 0) && (std::isfinite(determinant(deviatoric_strain_rate))) && (anisotropic_viscosity != nullptr)) // && (in.temperature[q]>1000)
            {
              // Get eigenvalues from compositional fields
              const std::vector<double> &composition = in.composition[q];
              const double eigvalue_a1 = composition[cpo_bingham_avg_a[1]];
              const double eigvalue_b1 = composition[cpo_bingham_avg_b[1]];
              const double eigvalue_c1 = composition[cpo_bingham_avg_c[1]];
              const double eigvalue_a2 = composition[cpo_bingham_avg_a[2]];
              const double eigvalue_b2 = composition[cpo_bingham_avg_b[2]];
              const double eigvalue_c2 = composition[cpo_bingham_avg_c[2]];
              const double eigvalue_a3 = composition[cpo_bingham_avg_a[3]];
              const double eigvalue_b3 = composition[cpo_bingham_avg_b[3]];
              const double eigvalue_c3 = composition[cpo_bingham_avg_c[3]];

              // Calculate the rotation matrix from the euler angles
              const double phi1 = composition[cpo_bingham_avg_a[0]];
              const double theta = composition[cpo_bingham_avg_b[0]];
              const double phi2 = composition[cpo_bingham_avg_c[0]];

              const Tensor<2,3> R = euler_angles_to_rotation_matrix(phi1, theta, phi2);

              // initialize scalar viscosity
              double scalar_viscosity = composition[viscosity_field_index];

              // Compute Hill Parameters FGHLMN from the eigenvalues of a,b,c axis
              // CPO2Hill v3 model:
              const double F = Utilities::fixed_power<2>(eigvalue_a1)*CnI_F[0] + eigvalue_a1*CnI_F[1] + eigvalue_a2*CnI_F[2] + (1/eigvalue_a3)*CnI_F[3] + Utilities::fixed_power<2>(eigvalue_b1)*CnI_F[4] + eigvalue_b1*CnI_F[5] + eigvalue_b2*CnI_F[6] + (1/eigvalue_b3)*CnI_F[7] + Utilities::fixed_power<2>(eigvalue_c1)*CnI_F[8] + eigvalue_c1*CnI_F[9] + eigvalue_c2*CnI_F[10] + (1/eigvalue_c3)*CnI_F[11] + CnI_F[12];
              const double G = Utilities::fixed_power<2>(eigvalue_a1)*CnI_G[0] + eigvalue_a1*CnI_G[1] + eigvalue_a2*CnI_G[2] + (1/eigvalue_a3)*CnI_G[3] + Utilities::fixed_power<2>(eigvalue_b1)*CnI_G[4] + eigvalue_b1*CnI_G[5] + eigvalue_b2*CnI_G[6] + (1/eigvalue_b3)*CnI_G[7] + Utilities::fixed_power<2>(eigvalue_c1)*CnI_G[8] + eigvalue_c1*CnI_G[9] + eigvalue_c2*CnI_G[10] + (1/eigvalue_c3)*CnI_G[11] + CnI_G[12];
              const double H = Utilities::fixed_power<2>(eigvalue_a1)*CnI_H[0] + eigvalue_a1*CnI_H[1] + eigvalue_a2*CnI_H[2] + (1/eigvalue_a3)*CnI_H[3] + Utilities::fixed_power<2>(eigvalue_b1)*CnI_H[4] + eigvalue_b1*CnI_H[5] + eigvalue_b2*CnI_H[6] + (1/eigvalue_b3)*CnI_H[7] + Utilities::fixed_power<2>(eigvalue_c1)*CnI_H[8] + eigvalue_c1*CnI_H[9] + eigvalue_c2*CnI_H[10] + (1/eigvalue_c3)*CnI_H[11] + CnI_H[12];
              const double L = std::abs(Utilities::fixed_power<2>(eigvalue_a1)*CnI_L[0] + eigvalue_a1*CnI_L[1] + eigvalue_a2*CnI_L[2] + (1/eigvalue_a3)*CnI_L[3] + Utilities::fixed_power<2>(eigvalue_b1)*CnI_L[4] + eigvalue_b1*CnI_L[5] + eigvalue_b2*CnI_L[6] + (1/eigvalue_b3)*CnI_L[7] + Utilities::fixed_power<2>(eigvalue_c1)*CnI_L[8] + eigvalue_c1*CnI_L[9] + eigvalue_c2*CnI_L[10] + (1/eigvalue_c3)*CnI_L[11] + CnI_L[12]);
              const double M = std::abs(Utilities::fixed_power<2>(eigvalue_a1)*CnI_M[0] + eigvalue_a1*CnI_M[1] + eigvalue_a2*CnI_M[2] + (1/eigvalue_a3)*CnI_M[3] + Utilities::fixed_power<2>(eigvalue_b1)*CnI_M[4] + eigvalue_b1*CnI_M[5] + eigvalue_b2*CnI_M[6] + (1/eigvalue_b3)*CnI_M[7] + Utilities::fixed_power<2>(eigvalue_c1)*CnI_M[8] + eigvalue_c1*CnI_M[9] + eigvalue_c2*CnI_M[10] + (1/eigvalue_c3)*CnI_M[11] + CnI_M[12]);
              const double N = std::abs(Utilities::fixed_power<2>(eigvalue_a1)*CnI_N[0] + eigvalue_a1*CnI_N[1] + eigvalue_a2*CnI_N[2] + (1/eigvalue_a3)*CnI_N[3] + Utilities::fixed_power<2>(eigvalue_b1)*CnI_N[4] + eigvalue_b1*CnI_N[5] + eigvalue_b2*CnI_N[6] + (1/eigvalue_b3)*CnI_N[7] + Utilities::fixed_power<2>(eigvalue_c1)*CnI_N[8] + eigvalue_c1*CnI_N[9] + eigvalue_c2*CnI_N[10] + (1/eigvalue_c3)*CnI_N[11] + CnI_N[12]);

              // rotation matrix in Kelvin/Mandel notation
              Tensor<2,6> R_CPO_K;
              R_CPO_K[0][0] = Utilities::fixed_power<2>(R[0][0]);
              R_CPO_K[0][1] = Utilities::fixed_power<2>(R[0][1]);
              R_CPO_K[0][2] = Utilities::fixed_power<2>(R[0][2]);
              R_CPO_K[0][3] = numbers::SQRT2*R[0][1]*R[0][2];
              R_CPO_K[0][4] = numbers::SQRT2*R[0][0]*R[0][2];
              R_CPO_K[0][5] = numbers::SQRT2*R[0][0]*R[0][1];

              R_CPO_K[1][0] = Utilities::fixed_power<2>(R[1][0]);
              R_CPO_K[1][1] = Utilities::fixed_power<2>(R[1][1]);
              R_CPO_K[1][2] = Utilities::fixed_power<2>(R[1][2]);
              R_CPO_K[1][3] = numbers::SQRT2*R[1][1]*R[1][2];
              R_CPO_K[1][4] = numbers::SQRT2*R[1][0]*R[1][2];
              R_CPO_K[1][5] = numbers::SQRT2*R[1][0]*R[1][1];

              R_CPO_K[2][0] = Utilities::fixed_power<2>(R[2][0]);
              R_CPO_K[2][1] = Utilities::fixed_power<2>(R[2][1]);
              R_CPO_K[2][2] = Utilities::fixed_power<2>(R[2][2]);
              R_CPO_K[2][3] = numbers::SQRT2*R[2][1]*R[2][2];
              R_CPO_K[2][4] = numbers::SQRT2*R[2][0]*R[2][2];
              R_CPO_K[2][5] = numbers::SQRT2*R[2][0]*R[2][1];

              R_CPO_K[3][0] = numbers::SQRT2*R[1][0]*R[2][0];
              R_CPO_K[3][1] = numbers::SQRT2*R[1][1]*R[2][1];
              R_CPO_K[3][2] = numbers::SQRT2*R[1][2]*R[2][2];
              R_CPO_K[3][3] = R[1][1]*R[2][2]+R[1][2]*R[2][1];
              R_CPO_K[3][4] = R[1][0]*R[2][2]+R[1][2]*R[2][0];
              R_CPO_K[3][5] = R[1][0]*R[2][1]+R[1][1]*R[2][0];

              R_CPO_K[4][0] = numbers::SQRT2*R[0][0]*R[2][0];
              R_CPO_K[4][1] = numbers::SQRT2*R[0][1]*R[2][1];
              R_CPO_K[4][2] = numbers::SQRT2*R[0][2]*R[2][2];
              R_CPO_K[4][3] = R[0][1]*R[2][2]+R[0][2]*R[2][1];
              R_CPO_K[4][4] = R[0][0]*R[2][2]+R[0][2]*R[2][0];
              R_CPO_K[4][5] = R[0][0]*R[2][1]+R[0][1]*R[2][0];

              R_CPO_K[5][0] = numbers::SQRT2*R[0][0]*R[1][0];
              R_CPO_K[5][1] = numbers::SQRT2*R[0][1]*R[1][1];
              R_CPO_K[5][2] = numbers::SQRT2*R[0][2]*R[1][2];
              R_CPO_K[5][3] = R[0][1]*R[1][2]+R[0][2]*R[1][1];
              R_CPO_K[5][4] = R[0][0]*R[1][2]+R[0][2]*R[1][0];
              R_CPO_K[5][5] = R[0][0]*R[1][1]+R[0][1]*R[1][0];

              // using the analytical inversion based on principal invariants of orthotropic symmetry group
              if (use_analytical_inversion == true)
                {
                  // rotate 3d strain-rate into cpo frame and compute scalar anisotropic viscosity
                  Tensor<2,3> strain_rate_cpo = R*strain_rate_3d*transpose(R);
                  scalar_viscosity =  std::pow(Gamma,(-1/n))*std::pow(CPO_AV_3D::orthotropic_strain_rate_invariant(strain_rate_cpo, F,G,H, M,N,L), ((1-n)/n) );

                  // calculate the anisotropic tensor for viscosity in CPO frame and rotate into model frame
                  const Tensor<2,6> viscosity_tensor = transpose(R_CPO_K) * CPO_AV_3D::viscosity_tensor_cpo_frame(F,G,H, M,N,L) * R_CPO_K;

                  // save viscosity tensor in stress-strain director to be used in
                  anisotropic_viscosity->stress_strain_directors[q] = CPO_AV_3D::kelvin_to_r4_tensor(viscosity_tensor);
                }
              else // using the iterative inversion based on a pseudo inverse of the anisotropic tensor for fluidity
                {
                  SymmetricTensor<2,6> A;
                  A[0][0] = (2./3.) * (G+H);
                  A[0][1] = (2./3.) * (-H);
                  A[0][2] = (2./3.) * (-G);
                  A[1][1] = (2./3.) * (H+F);
                  A[1][2] = (2./3.) * (-F);
                  A[2][2] = (2./3.) * (F+G);
                  A[3][3] = (2./3.) * (L);
                  A[4][4] = (2./3.) * (M);
                  A[5][5] = (2./3.) * (N);

                  // A is the anisotropic tensor for the fluidity. We need its inverse, but it's not invertible due to singularity.
                  // Thus we compute the Moore-Penrose pseudo inverse using SVD
                  LAPACKFullMatrix<double> A_mat_lapack(6,6), pinvA_mat_lapack(6,6);
                  for (unsigned int ai=0; ai<6; ++ai)
                    {
                      for (unsigned int aj=0; aj<6; ++aj)
                        {
                          A_mat_lapack(ai,aj) = A[ai][aj];
                        }
                    }
                  pseudoinverse(A_mat_lapack, pinvA_mat_lapack);

                  SymmetricTensor<2,6> invA;
                  for (unsigned int ai=0; ai<6; ++ai)
                    {
                      for (unsigned int aj=0; aj<6; ++aj)
                        {
                          invA[ai][aj] = pinvA_mat_lapack(ai,aj);
                        }
                    }

                  // Calculate the viscosity tensor in the CPO frame
                  const Tensor<2,6> viscosity_tensor = transpose(R_CPO_K) * invA * R_CPO_K;

                  // save 2D (3D) viscosity tensor in stress strain-director
                  anisotropic_viscosity->stress_strain_directors[q] = CPO_AV_3D::kelvin_to_r4_tensor(viscosity_tensor);

                  // Convert rank 2 viscosity tensor to rank 4 necessary for iterative inversion
                  FullMatrix<double> viscosity_mat(6,6);
                  for (unsigned int vi=0; vi<6; ++vi)
                    {
                      for (unsigned int vj=0; vj<6; ++vj)
                        {
                          viscosity_mat[vi][vj] = viscosity_tensor[vi][vj];
                        }
                    }

                  // 3d viscosity tensor not always passed on to assembler
                  SymmetricTensor<4,3> viscosity_tensor_3D_r4;
                  dealii::Physics::Notation::Kelvin::to_tensor(viscosity_mat, viscosity_tensor_3D_r4);

                  // calculating anisotropic scalar viscosity from iterative inversion

                  // In the first time step using the actual strain rate can lead to convergence issue if the strain rate varies significantly within the model domain.
                  // Thus for the first timestep we calculate an initial viscosity based on the strain rate.
                  if (this->get_timestep_number() == 1)
                    {
                      double edot_ii=std::max(std::max(deviatoric_strain_rate.norm(), 0.),
                                              min_strain_rate);
                      scalar_viscosity= std::pow(Gamma,-1/n)*std::pow(edot_ii,((1. - n)/n));
                    }

                  unsigned int n_iterations = 0;
                  double residual = scalar_viscosity;
                  double threshold = relative_tolerance*scalar_viscosity;
                  // Here we convert stress to MPa to be consistent with the constitutive equation defined in Signorelli et al. (2021),
                  // in which the stress is in MPa.
                  SymmetricTensor<2,3> stress = scalar_viscosity * viscosity_tensor_3D_r4 * deviatoric_strain_rate; // 2 * / 1e6;

                  while (std::abs(residual) > threshold && n_iterations < max_iteration)
                    {
                      stress = (1./2.) * (stress + scalar_viscosity * viscosity_tensor_3D_r4 * deviatoric_strain_rate); // 2* / 1e6);

                      const Tensor<2,3> S_CPO= R * stress * transpose(R);

                      double Jhill = 2.0/3.0*(F*Utilities::fixed_power<2>(S_CPO[2][2]-S_CPO[1][1]) + G*Utilities::fixed_power<2>(S_CPO[0][0]-S_CPO[2][2]) + H*Utilities::fixed_power<2>(S_CPO[1][1]-S_CPO[0][0]) + 2*L*Utilities::fixed_power<2>(S_CPO[1][2]) + 2*M*Utilities::fixed_power<2>(S_CPO[0][2]) + 2*N*Utilities::fixed_power<2>(S_CPO[0][1]));
                      if (Jhill < 0)
                        {
                          Jhill = 2.0/3.0*(std::abs(F)*Utilities::fixed_power<2>(S_CPO[2][2]-S_CPO[1][1]) + std::abs(G)*Utilities::fixed_power<2>(S_CPO[0][0]-S_CPO[2][2]) + std::abs(H)*Utilities::fixed_power<2>(S_CPO[1][1]-S_CPO[0][0]) + 2*L*Utilities::fixed_power<2>(S_CPO[1][2]) + 2*M*Utilities::fixed_power<2>(S_CPO[0][2]) + 2*N*Utilities::fixed_power<2>(S_CPO[0][1]));
                        }

                      AssertThrow(std::isfinite(Jhill),
                                  ExcMessage("Jhill should be finite"));
                      AssertThrow(Jhill >= 0,
                                  ExcMessage("Jhill should not be negative"));

                      const double scalar_viscosity_new = (1 / (Gamma * std::pow(Jhill,(n-1)/2)));
                      residual = std::abs(scalar_viscosity_new - scalar_viscosity);
                      scalar_viscosity = scalar_viscosity_new;
                      threshold = relative_tolerance*scalar_viscosity;
                      n_iterations++;
                    }
                }

              out.viscosities[q] = scalar_viscosity;

              AssertThrow(std::isfinite(out.viscosities[q]),
                          ExcMessage("Viscosity should be finite"));
              AssertThrow(out.viscosities[q] > 0,
                          ExcMessage("Viscosity should be positive"));

            }
          else // timestep == 0 or no anisotropic viscosity
            {
              if ((this->simulator_is_past_initialization()) && (std::isfinite(determinant(deviatoric_strain_rate))))
                {
                  // for the zero-th timestep calculating the scalar viscosity based on the strain-rate -> i.e. isotropic response
                  double edot_ii=std::max(std::max(deviatoric_strain_rate.norm(), 0.),
                                          min_strain_rate);
                  out.viscosities[q] = std::pow(Gamma, (-1/n))*std::pow(edot_ii,((1. - n)/n)); //

                }

              if (anisotropic_viscosity != nullptr)
                {
                  if ((this->simulator_is_past_initialization()) && (std::isfinite(determinant(deviatoric_strain_rate))))
                    {
                      // for the zero-th timestep calculating the scalar viscosity based on the strain-rate -> i.e. isotropic response
                      // double edot_ii=std::max(std::max(deviatoric_strain_rate.norm(), 0.),
                      //                               min_strain_rate);
                      // out.viscosities[q] = 1/Gamma * std::pow(edot_ii,((1. - n)/n)); //
                      double edot_ii=std::max(std::sqrt(std::max(-second_invariant(deviator(strain_rate)), 0.)),
                                              min_strain_rate);
                      out.viscosities[q] = 1/Gamma * std::pow(edot_ii,((1. - n)/n));
                    }
                  // Assign an isotropic viscosity tensor
                  SymmetricTensor<2,6> viscosity_tensor;
                  viscosity_tensor[0][0] = 2.0/3.0;  // 4.0/9.0;
                  viscosity_tensor[0][1] = -1.0/3.0; // -2.0/9.0;
                  viscosity_tensor[0][2] = -1.0/3; // -2.0/9.0;
                  viscosity_tensor[1][1] = 2.0/3.0;  // 4.0/9.0;
                  viscosity_tensor[1][2] = -1.0/3.0; // -2.0/9.0;
                  viscosity_tensor[2][2] = 2.0/3.0;
                  viscosity_tensor[3][3] = 1.0;
                  viscosity_tensor[4][4] = 1.0;
                  viscosity_tensor[5][5] = 1.0;

                  // save viscosity tensor in stress-strain director to be used in
                  anisotropic_viscosity->stress_strain_directors[q] = CPO_AV_3D::kelvin_to_r4_tensor(viscosity_tensor);

                }
            }

          // Prescribe the scalar viscosity to compositional field for access in the next time step
          if (const std::shared_ptr<PrescribedFieldOutputs<dim>> prescribed_field_out
              = out.template get_additional_output_object<PrescribedFieldOutputs<dim>>())
            {
              prescribed_field_out->prescribed_field_outputs[q][viscosity_field_index] = out.viscosities[q];
            }
        }
    }



    template <int dim>
    bool
    CPO_AV_3D<dim>::is_compressible () const
    {
      return false;
    }



    template <int dim>
    void
    CPO_AV_3D<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("CPO-induced Anisotropic Viscosity");
        {
          EquationOfState::LinearizedIncompressible<dim>::declare_parameters (prm);

          prm.declare_entry ("Coefficients and intercept for F", "0.5920219168461529, -0.831936049, -0.000937583, -0.00029648, 0.380413345, -0.533048795, 0.46835862365145753, -0.000965503, -1.249340274, 1.0748554477472438, -0.167662132, 0.003358407, 0.5215020195972386",
                             Patterns::List(Patterns::Double()),
                             "12 Coefficients and 1 intercept to compute the Hill Parameter F "
                             "according to the linear regression relation provided in the cookbook documentation. "
                             "The first 3 coefficients are multiplied respectively by: "
                             "the square of the largest eigenvalue, the second-largest eigenvalue, "
                             "and the inverse of the smallest eigenvalue of the a-axis orientation tensor. "
                             "The next 3 coefficients are used in the same way for the b-axis, "
                             "and the final 3 for the c-axis. Together with the intercept, "
                             "these values form the full regression expression for F.");
          prm.declare_entry ("Coefficients and intercept for G", "-1.6951323, 1.3364976547800977, -0.18410694, 4.918308228973878e-05, 0.7501414478371807, 0.691412915, 0.37696216069289673, -0.001537058, -0.66969478, -0.551507796, -0.428462988, 0.003403174, 0.2602865312446454",
                             Patterns::List(Patterns::Double()),
                             "12 Coefficients and 1 intercept to compute the Hill Parameter G in the same way as above.");
          prm.declare_entry ("Coefficients and intercept for H", "-1.139612048, 1.353113344145978, 0.7510486623018213, -0.001656848, -0.255721452, -1.006433455, -0.11595106, 0.003177149, 0.6837240306536184, -0.031162568, -0.080356281, 0.005621241, 0.26788414888354445",
                             Patterns::List(Patterns::Double()),
                             "12 Coefficients and 1 intercept to compute the Hill Parameter H in the same way as above.");
          prm.declare_entry ("Coefficients and intercept for L", "-3.510950516, 2.6864808033885543, 0.035838123, -0.000504338, 3.9483598066088383, -3.816102334, -0.778569714, 0.003688104, 4.122460734346824, -2.482527095, 1.3200590504614733, -0.002399896, 2.0027068994912076",
                             Patterns::List(Patterns::Double()),
                             "12 Coefficients and 1 intercept to compute the Hill Parameter L in the same way as above.");
          prm.declare_entry ("Coefficients and intercept for M", "4.536980494567378, -3.227568914, 0.27609495132676254, 0.007436169, -7.446913908, 5.763821882498737, -1.4026181, 0.032132134, 2.9678024468288697, -3.434721081, -2.265560577, 0.1215179917888699, 2.4409386788998457",
                             Patterns::List(Patterns::Double()),
                             "12 Coefficients and 1 intercept to compute the Hill Parameter M in the same way as above.");
          prm.declare_entry ("Coefficients and intercept for N", "7.872922986831338, -7.933513948, -2.588175191, 0.029843804645040883, 7.605864624755694, -5.469451775, -0.347688637, 0.06395131, -1.78763278, 2.2550636824173584, 3.023166891521831, -0.102862765, 3.6958741003498234",
                             Patterns::List(Patterns::Double()),
                             "12 Coefficients and 1 intercept to compute the Hill Parameter N in the same way as above.");

          prm.declare_entry ("Reference viscosity", "1e9",
                             Patterns::Double(),
                             "Magnitude of reference viscosity.");
          prm.declare_entry ("Minimum strain rate", "1.4e-20", Patterns::Double(),
                             "Stabilizes strain dependent viscosity. Units: \\si{\\per\\second}");
          prm.declare_entry ("Grain size", "1.0e-3",
                             Patterns::Double(),
                             "Olivine anisotropic viscosity is dependent of grain size. Value is given in meters");
          prm.declare_entry ("Fluidity constant", "1.1e5",
                             Patterns::Double(),
                             "Prefactor for Arhenius temperature activation");
          prm.declare_entry ("Activation energy", "530000.0",
                             Patterns::Double(),
                             "Activation energy for Arhenius temperature dependence of rheology");
          prm.declare_entry ("Grain size exponent", "0.73",
                             Patterns::Double(),
                             "Exponent for grainsize dependence");
          prm.declare_entry ("Stress exponent", "3.5",
                             Patterns::Double(),
                             "Stress exponent for non-linear rheology");
          prm.declare_entry ("Use analytical inversion", "false",
                             Patterns::Bool (),
                             "Whether to use the analytical or the iterative inversion for the anisotropic scalar viscosity.");
          prm.declare_entry ("Relative tolerance for iteration", "0.0001",
                             Patterns::Double(),
                             "The iteration for computing scalar viscosity is terminated when the relative change falls below the relative tolerance.");
          prm.declare_entry ("Maximum number of iterations", "100",
                             Patterns::Integer(),
                             "To prevent excessive computation, the number of iterations with a maximum number of iterations.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

    }



    template <int dim>
    void
    CPO_AV_3D<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("CPO-induced Anisotropic Viscosity");
        {
          equation_of_state.parse_parameters (prm);
          eta = prm.get_double("Reference viscosity");
          min_strain_rate = prm.get_double("Minimum strain rate");
          grain_size = prm.get_double("Grain size");
          stress_exponent = prm.get_double("Stress exponent");
          fluidity_constant = prm.get_double("Fluidity constant");
          grain_size_exponent = prm.get_double("Grain size exponent");
          activation_energy = prm.get_double("Activation energy");
          use_analytical_inversion  = prm.get_bool ("Use analytical inversion");
          relative_tolerance = prm.get_double("Relative tolerance for iteration");
          max_iteration = prm.get_integer("Maximum number of iterations");
          CnI_F = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(prm.get("Coefficients and intercept for F")));
          CnI_G = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(prm.get("Coefficients and intercept for G")));
          CnI_H = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(prm.get("Coefficients and intercept for H")));
          CnI_L = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(prm.get("Coefficients and intercept for L")));
          CnI_M = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(prm.get("Coefficients and intercept for M")));
          CnI_N = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(prm.get("Coefficients and intercept for N")));
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // Declare dependence
      this->model_dependence.density = NonlinearDependence::compositional_fields;
    }



    template <int dim>
    void
    CPO_AV_3D<dim>::create_additional_named_outputs(MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      if (out.template get_additional_output_object<AnisotropicViscosity<dim>>() == nullptr)
        {
          const unsigned int n_points = out.n_evaluation_points();
          out.additional_outputs.push_back(
            std::make_unique<MaterialModel::AnisotropicViscosity<dim>> (n_points));
        }

      if (out.template get_additional_output_object<PrescribedFieldOutputs<dim>>() == NULL)
        {
          const unsigned int n_points = out.n_evaluation_points();
          out.additional_outputs.push_back(
            std::make_unique<MaterialModel::PrescribedFieldOutputs<dim>> (n_points,this->n_compositional_fields()));
        }
    }
  }
}



// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(CPO_AV_3D,
                                   "CPO-induced anisotropic viscosity",
                                   "Olivine CPO related viscous anisotropy based on the Simple material model. "
                                   "For more details see the documentation of this class in its header file and the cookbook documentation.")
  }
}
