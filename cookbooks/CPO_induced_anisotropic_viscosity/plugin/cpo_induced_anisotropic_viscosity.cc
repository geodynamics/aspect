
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
      AssertThrow(dim==3,
                  ExcMessage("Olivine has 3 independent slip systems, allowing for deformation in 3 independent directions, hence these models only work in 3D"));

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



    template <>
    void
    CPO_AV_3D<2>::evaluate (const MaterialModel::MaterialModelInputs<2> &,
                            MaterialModel::MaterialModelOutputs<2> &) const
    {
      Assert (false, ExcNotImplemented());
    }



    template <>
    void
    CPO_AV_3D<3>::evaluate (const MaterialModel::MaterialModelInputs<3> &in,
                            MaterialModel::MaterialModelOutputs<3> &out) const
    {
      const int dim=3;
      const std::shared_ptr<MaterialModel::AnisotropicViscosity<dim>> anisotropic_viscosity =
        out.template get_additional_output_object<MaterialModel::AnisotropicViscosity<dim>>();
      EquationOfStateOutputs<dim> eos_outputs (1);
      const unsigned int viscosity_field_index = this->introspection().compositional_index_for_name("scalar_viscosity");

      for (unsigned int q=0; q<in.n_evaluation_points(); ++q)
        {
          equation_of_state.evaluate(in, q, eos_outputs);

          // Get parameters for compute the effective viscosity
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
          const SymmetricTensor<2,dim> deviatoric_strain_rate
            = (this->get_material_model().is_compressible()
               ?
               strain_rate - 1./3. * trace(strain_rate) * unit_symmetric_tensor<dim>()
               :
               strain_rate);

          // The computation of the viscosity tensor is only necessary after the simulator has been initialized
          // and when the condition allows dislocation creep
          if  ((this->simulator_is_past_initialization()) && (this->get_timestep_number() > 0) && (in.temperature[q]>1000) && (std::isfinite(determinant(deviatoric_strain_rate))) && (anisotropic_viscosity != nullptr))
            {
              // Create constant value to use for AV
              const double A_o = 1.1e5*std::exp(-530000/(8.314*in.temperature[q]));
              const double n = 3.5;
              // The values of A_o and 0.73 were picked so that Gamma = 3.5322e-15[1/(s*Pa^n)] if T=1600K and d=1000 microns
              const double Gamma = (A_o/(std::pow(grain_size,0.73)));

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

              const Tensor<2,3> R = transpose(euler_angles_to_rotation_matrix(phi1, theta, phi2));

              // Compute Hill Parameters FGHLMN from the eigenvalues of a,b,c axis
              // CPO2Hill v5 model:
              const double F = Utilities::fixed_power<2>(eigvalue_a1)*CnI_F[0] + eigvalue_a2*CnI_F[1] + (1/eigvalue_a3)*CnI_F[2] + Utilities::fixed_power<2>(eigvalue_b1)*CnI_F[3] + eigvalue_b2*CnI_F[4] + (1/eigvalue_b3)*CnI_F[5] + Utilities::fixed_power<2>(eigvalue_c1)*CnI_F[6] + eigvalue_c2*CnI_F[7] + (1/eigvalue_c3)*CnI_F[8] + CnI_F[9];
              const double G = Utilities::fixed_power<2>(eigvalue_a1)*CnI_G[0] + eigvalue_a2*CnI_G[1] + (1/eigvalue_a3)*CnI_G[2] + Utilities::fixed_power<2>(eigvalue_b1)*CnI_G[3] + eigvalue_b2*CnI_G[4] + (1/eigvalue_b3)*CnI_G[5] + Utilities::fixed_power<2>(eigvalue_c1)*CnI_G[6] + eigvalue_c2*CnI_G[7] + (1/eigvalue_c3)*CnI_G[8] + CnI_G[9];
              const double H = Utilities::fixed_power<2>(eigvalue_a1)*CnI_H[0] + eigvalue_a2*CnI_H[1] + (1/eigvalue_a3)*CnI_H[2] + Utilities::fixed_power<2>(eigvalue_b1)*CnI_H[3] + eigvalue_b2*CnI_H[4] + (1/eigvalue_b3)*CnI_H[5] + Utilities::fixed_power<2>(eigvalue_c1)*CnI_H[6] + eigvalue_c2*CnI_H[7] + (1/eigvalue_c3)*CnI_H[8] + CnI_H[9];
              const double L = std::abs(Utilities::fixed_power<2>(eigvalue_a1)*CnI_L[0] + eigvalue_a2*CnI_L[1] + (1/eigvalue_a3)*CnI_L[2] + Utilities::fixed_power<2>(eigvalue_b1)*CnI_L[3] + eigvalue_b2*CnI_L[4] + (1/eigvalue_b3)*CnI_L[5] + Utilities::fixed_power<2>(eigvalue_c1)*CnI_L[6] + eigvalue_c2*CnI_L[7] + (1/eigvalue_c3)*CnI_L[8] + CnI_L[9]);
              const double M = std::abs(Utilities::fixed_power<2>(eigvalue_a1)*CnI_M[0] + eigvalue_a2*CnI_M[1] + (1/eigvalue_a3)*CnI_M[2] + Utilities::fixed_power<2>(eigvalue_b1)*CnI_M[3] + eigvalue_b2*CnI_M[4] + (1/eigvalue_b3)*CnI_M[5] + Utilities::fixed_power<2>(eigvalue_c1)*CnI_M[6] + eigvalue_c2*CnI_M[7] + (1/eigvalue_c3)*CnI_M[8] + CnI_M[9]);
              const double N = std::abs(Utilities::fixed_power<2>(eigvalue_a1)*CnI_N[0] + eigvalue_a2*CnI_N[1] + (1/eigvalue_a3)*CnI_N[2] + Utilities::fixed_power<2>(eigvalue_b1)*CnI_N[3] + eigvalue_b2*CnI_N[4] + (1/eigvalue_b3)*CnI_N[5] + Utilities::fixed_power<2>(eigvalue_c1)*CnI_N[6] + eigvalue_c2*CnI_N[7] + (1/eigvalue_c3)*CnI_N[8] + CnI_N[9]);

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

              SymmetricTensor<2,6> A;
              A[0][0] = 2.0/3.0*(F+H);
              A[0][1] = 2.0/3.0*(-F);
              A[0][2] = 2.0/3.0*(-H);
              A[1][1] = 2.0/3.0*(G+F);
              A[1][2] = 2.0/3.0*(-G);
              A[2][2] = 2.0/3.0*(H+G);
              A[3][3] = 2.0/3.0*L;
              A[4][4] = 2.0/3.0*M;
              A[5][5] = 2.0/3.0*N;

              // A is the anisotropic tensor for the fluidity. We need its inverse, but it's not invertible due to singularity.
              // Thus we compute the Moore-Penrose pseudo inverse using SVD
              LAPACKFullMatrix<double> A_mat_lapack(6, 6), pinvA_mat_lapack(6,6);
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

              // Calculate the fluidity tensor in the CPO frame
              const Tensor<2,6> V = R_CPO_K * invA * transpose(R_CPO_K);

              // Convert rank 2 viscosity tensor to rank 4
              FullMatrix<double> V_mat(6,6);
              for (unsigned int vi=0; vi<6; ++vi)
                {
                  for (unsigned int vj=0; vj<6; ++vj)
                    {
                      V_mat[vi][vj] = V[vi][vj];
                    }
                }
              SymmetricTensor<4,dim> V_r4;
              dealii::Physics::Notation::Kelvin::to_tensor(V_mat, V_r4);
              anisotropic_viscosity->stress_strain_directors[q] = V_r4;

              double scalar_viscosity = composition[viscosity_field_index];

              // In the first time step using the actual strain rate can lead to convergence issue if the strain rate varies significantly within the model domain.
              // Thus for the first timestep we calculate an initial viscosity based on the strain rate.
              if (this->get_timestep_number() == 1)
                {
                  const double edot_ii=std::max(std::sqrt(std::max(-second_invariant(deviator(strain_rate)), 0.)),
                                                min_strain_rate);
                  scalar_viscosity= 1/Gamma * std::pow(edot_ii,((1. - n)/n));
                }

              unsigned int n_iterations = 0;
              const unsigned int max_iteration = 100;
              double residual = scalar_viscosity;
              double threshold = 0.0001*scalar_viscosity;
              // Here we convert stress to MPa to be consistent with the constitutive equation defined in Signorelli et al. (2021),
              // in which the stress is in MPa.
              SymmetricTensor<2,dim> stress = 2 * scalar_viscosity * V_r4 * deviatoric_strain_rate / 1e6;

              const Tensor<2,dim> R_T = transpose(R);
              while (std::abs(residual) > threshold && n_iterations < max_iteration)
                {
                  stress = (1./2.) * (stress + 2*scalar_viscosity * V_r4 * deviatoric_strain_rate / 1e6);

                  const Tensor<2,dim> S_CPO= R_T * stress * R;

                  double Jhill = F*Utilities::fixed_power<2>(S_CPO[0][0]-S_CPO[1][1]) + G*Utilities::fixed_power<2>(S_CPO[1][1]-S_CPO[2][2]) + H*Utilities::fixed_power<2>(S_CPO[2][2]-S_CPO[0][0]) + 2*L*Utilities::fixed_power<2>(S_CPO[1][2]) + 2*M*Utilities::fixed_power<2>(S_CPO[0][2]) + 2*N*Utilities::fixed_power<2>(S_CPO[0][1]);
                  if (Jhill < 0)
                    {
                      Jhill = std::abs(F)*Utilities::fixed_power<2>(S_CPO[0][0]-S_CPO[1][1]) + std::abs(G)*Utilities::fixed_power<2>(S_CPO[1][1]-S_CPO[2][2]) + std::abs(H)*Utilities::fixed_power<2>(S_CPO[2][2]-S_CPO[0][0]) + 2*L*Utilities::fixed_power<2>(S_CPO[1][2]) + 2*M*Utilities::fixed_power<2>(S_CPO[0][2]) + 2*N*Utilities::fixed_power<2>(S_CPO[0][1]);
                    }

                  AssertThrow(std::isfinite(Jhill),
                              ExcMessage("Jhill should be finite"));
                  AssertThrow(Jhill >= 0,
                              ExcMessage("Jhill should not be negative"));

                  const double scalar_viscosity_new = (1 / (Gamma * std::pow(Jhill,(n-1)/2)));
                  residual = std::abs(scalar_viscosity_new - scalar_viscosity);
                  scalar_viscosity = scalar_viscosity_new;
                  threshold = 0.001*scalar_viscosity;
                  n_iterations++;

                }
              // Store the scalar viscosity in out.viscosities
              out.viscosities[q] = scalar_viscosity;

              AssertThrow(std::isfinite(out.viscosities[q]),
                          ExcMessage("Viscosity should be finite"));
              AssertThrow(out.viscosities[q] > 0,
                          ExcMessage("Viscosity should be positive"));

            }
          else // timestep == 0 or no anisotropic viscosity
            {
              if (anisotropic_viscosity != nullptr)
                {
                  // Assign an isotropic viscosity tensor
                  SymmetricTensor<2,6> V;
                  V[0][0] = 2.0/3.0;
                  V[0][1] = -1.0/3.0;
                  V[0][2] = -1.0/3.0;
                  V[1][1] = 2.0/3.0;
                  V[1][2] = -1.0/3.0;
                  V[2][2] = 2.0/3.0;
                  V[3][3] = 1;
                  V[4][4] = 1;
                  V[5][5] = 1;

                  // Convert rank 2 viscosity tensor to rank 4
                  FullMatrix<double> V_mat(6,6);
                  for (unsigned int vi=0; vi<6; ++vi)
                    {
                      for (unsigned int vj=0; vj<6; ++vj)
                        {
                          V_mat[vi][vj] = V[vi][vj];
                        }
                    }
                  SymmetricTensor<4,dim> V_r4;
                  dealii::Physics::Notation::Kelvin::to_tensor(V_mat, V_r4);
                  anisotropic_viscosity->stress_strain_directors[q] = V_r4;
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

          prm.declare_entry ("Coefficients and intercept for F", "1.0390459583037057,  -0.767458622,  0.003066208,  0.19651133418307049,  0.413093763,  0.015463162,  -0.935925291,  -2.392877563,  0.051834768, 1.0799807050187482",
                             Patterns::List(Patterns::Double()),
                             "9 Coefficients and 1 intercept to compute the Hill Parameter F "
                             "according to the linear regression relation provided in the cookbook documentation. "
                             "The first 3 coefficients are multiplied respectively by: "
                             "the square of the largest eigenvalue, the second-largest eigenvalue, "
                             "and the inverse of the smallest eigenvalue of the a-axis orientation tensor. "
                             "The next 3 coefficients are used in the same way for the b-axis, "
                             "and the final 3 for the c-axis. Together with the intercept, "
                             "these values form the full regression expression for F.");
          prm.declare_entry ("Coefficients and intercept for G", "-2.836270315,  -1.632453092,  0.000687606,  0.2671850239576621,  -0.993392913,  0.002699241,  1.9689530759060374,  2.314442451425019,  -0.018655905, 0.6887411607403755",
                             Patterns::List(Patterns::Double()),
                             "9 Coefficients and 1 intercept to compute the Hill Parameter G in the same way as above.");
          prm.declare_entry ("Coefficients and intercept for H", "1.6687493021559732,  0.5797579293682223,  0.003241593,  0.701661336,  0.2513824481429968,  0.000229291,  -2.003227619,  -2.57032429,  0.071454541, 0.7490268673620638",
                             Patterns::List(Patterns::Double()),
                             "9 Coefficients and 1 intercept to compute the Hill Parameter H in the same way as above.");
          prm.declare_entry ("Coefficients and intercept for L", "-0.325145943,  0.7284642859944138,  0.000404879,  -0.665446098,  0.5152847961409479,  0.002722782,  -1.026786493,  -1.262574542,  0.009168498, 1.595422603",
                             Patterns::List(Patterns::Double()),
                             "9 Coefficients and 1 intercept to compute the Hill Parameter L in the same way as above.");
          prm.declare_entry ("Coefficients and intercept for M", "1.6427437063774875,  0.8777500120437522,  0.004651732,  2.489417876177839,  0.8162729707609052,  -0.010736521,  -2.49420455,  -0.511446494,  -0.009362491, 0.893677343",
                             Patterns::List(Patterns::Double()),
                             "9 Coefficients and 1 intercept to compute the Hill Parameter M in the same way as above.");
          prm.declare_entry ("Coefficients and intercept for N", "0.8122098589701904,  0.15663795996228266,  0.001500252,  -1.648578168,  0.19362392490527092,  -0.009650519,  1.6796559729985163,  -0.103640482,  0.01971017, 1.2132200780065174",
                             Patterns::List(Patterns::Double()),
                             "9 Coefficients and 1 intercept to compute the Hill Parameter N in the same way as above.");

          prm.declare_entry ("Reference viscosity", "1e9",
                             Patterns::Double(),
                             "Magnitude of reference viscosity.");
          prm.declare_entry ("Minimum strain rate", "1.4e-20", Patterns::Double(),
                             "Stabilizes strain dependent viscosity. Units: \\si{\\per\\second}");
          prm.declare_entry ("Grain size", "1e-3",
                             Patterns::Double(),
                             "Olivine anisotropic viscosity is dependent of grain size. Value is given in meters");
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
