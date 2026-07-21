/*
  Copyright (C) 2022 - 2024 by the authors of the ASPECT code.

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

#include "common.h"
#include <aspect/material_model/utilities.h>
#include <deal.II/base/parameter_handler.h>
#include <aspect/utilities.h>

// /**
//  * Compare two symmetric tensors
//  */
// inline void compare_tensors_approx(
//     const dealii::Tensor<2,3> &computed,
//     const dealii::Tensor<2,3> &expected)
// {
//   INFO("tensors are not the same: \n" <<
//        "expected = " << expected[0][0] << " " << expected[0][1] << " " << expected[0][2] << "\n" <<
//        "           " << expected[1][1] << " " << expected[1][2] << "\n" <<
//        "           " << expected[2][2] << "\n" <<
//        "computed = " << computed[0][0] << " " << computed[0][1] << " " << computed[0][2] << "\n" <<
//        "           " << computed[1][1] << " " << computed[1][2] << "\n" <<
//        "           " << computed[2][2] << "\n" );
//   const double tol = 1e-14;
//   CHECK(
//         ((computed[0][0] == Approx(expected[0][0]) || std::fabs(computed[0][0]) < tol)
//           && (computed[0][1] == Approx(expected[0][1]) || std::fabs(computed[0][1]) < tol)
//           && (computed[0][2] == Approx(expected[0][2]) || std::fabs(computed[0][2]) < tol)
//           && (computed[1][0] == Approx(expected[1][0]) || std::fabs(computed[1][0]) < tol)
//           && (computed[1][1] == Approx(expected[1][1]) || std::fabs(computed[1][1]) < tol)
//           && (computed[1][2] == Approx(expected[1][2]) || std::fabs(computed[1][2]) < tol)
//           && (computed[2][0] == Approx(expected[2][0]) || std::fabs(computed[2][0]) < tol)
//           && (computed[2][1] == Approx(expected[2][1]) || std::fabs(computed[2][1]) < tol)
//           && (computed[2][2] == Approx(expected[2][2]) || std::fabs(computed[2][2]) < tol))
//         );
// }

// /**
//  * Compare two symmetric rank 4 tensors in kelvin notation
//  */
// inline void compare_kelvin_sym_tensors_approx(
//     const dealii::SymmetricTensor<2,6> &computed,
//     const dealii::SymmetricTensor<2,6> &expected)
// {
//   INFO("tensors are not the same: \n" <<
//        "expected = " << expected[0][0] << " " << expected[0][1] << " " << expected[0][2] << expected[0][3] << " " << expected[0][4] << " " << expected[0][5] << "\n" <<
//        "           " <<  expected[1][1] << " " << expected[1][2] << expected[1][3] << " " << expected[1][4] << " " << expected[1][5] <<"\n" <<
//        "           " << expected[2][2] << expected[2][3] << " " << expected[2][4] << " " << expected[2][5] << "\n" <<
//        "           " << expected[3][3] << expected[3][4] << " " << expected[3][5] << "\n" <<
//        "           " << expected[4][4] << expected[4][5] << "\n" <<
//        "           " << expected[5][5] << "\n" <<

//        "computed = " << computed[0][0] << " " << computed[0][1] << " " << computed[0][2] << computed[0][3] << " " << computed[0][4] << " " << computed[0][5] << "\n" <<
//        "           " << computed[1][1] << " " << computed[1][2] << computed[1][3] << " " << computed[1][4] << " " << computed[1][5] <<"\n" <<
//        "           " << computed[2][2] << computed[2][3] << " " << computed[2][4] << " " << computed[2][5] << "\n" <<
//        "           " << computed[3][3] << computed[3][4] << " " << computed[3][5] << "\n" <<
//        "           " << computed[4][4] << computed[4][5] << "\n" <<
//        "           " << computed[5][5] << "\n" );
//   const double tol = 1e-14;
//   CHECK(
//         ((    computed[0][0] == Approx(expected[0][0]) || std::fabs(computed[0][0]) < tol)
//           && (computed[0][1] == Approx(expected[0][1]) || std::fabs(computed[0][1]) < tol)
//           && (computed[0][2] == Approx(expected[0][2]) || std::fabs(computed[0][2]) < tol)
//           && (computed[0][3] == Approx(expected[0][3]) || std::fabs(computed[0][3]) < tol)
//           && (computed[0][4] == Approx(expected[0][4]) || std::fabs(computed[0][4]) < tol)
//           && (computed[0][5] == Approx(expected[0][5]) || std::fabs(computed[0][5]) < tol)
          
//           && (computed[1][1] == Approx(expected[1][1]) || std::fabs(computed[1][1]) < tol)
//           && (computed[1][2] == Approx(expected[1][2]) || std::fabs(computed[1][2]) < tol)
//           && (computed[1][3] == Approx(expected[1][3]) || std::fabs(computed[1][3]) < tol)
//           && (computed[1][4] == Approx(expected[1][4]) || std::fabs(computed[1][4]) < tol)
//           && (computed[1][5] == Approx(expected[1][5]) || std::fabs(computed[1][5]) < tol)
          
//           && (computed[2][2] == Approx(expected[2][2]) || std::fabs(computed[2][2]) < tol)
//           && (computed[2][3] == Approx(expected[2][3]) || std::fabs(computed[2][3]) < tol)
//           && (computed[2][4] == Approx(expected[2][4]) || std::fabs(computed[2][4]) < tol)
//           && (computed[2][5] == Approx(expected[2][5]) || std::fabs(computed[2][5]) < tol)
          
//           && (computed[3][3] == Approx(expected[3][3]) || std::fabs(computed[3][3]) < tol)
//           && (computed[3][4] == Approx(expected[3][4]) || std::fabs(computed[3][4]) < tol)
//           && (computed[3][5] == Approx(expected[3][5]) || std::fabs(computed[3][5]) < tol)
          
//           && (computed[4][4] == Approx(expected[4][4]) || std::fabs(computed[4][4]) < tol)
//           && (computed[4][5] == Approx(expected[4][5]) || std::fabs(computed[4][5]) < tol)
          
//           && (computed[5][5] == Approx(expected[5][5]) || std::fabs(computed[5][5]) < tol)
//         )
          
//       );
// }


TEST_CASE("OrthotropicRheology Functions")
{
  using namespace aspect;
  using namespace dealii; // for symmetric tensors
  
  // std::vector<int> ji = {1,2,0}; // tuple of indices shifted by one
  // std::vector<int> ki = {2,0,1}; // tuple of indices shifted by two

  // SymmetricTensor<2,3> test_stress; 
  // SymmetricTensor<2,3> test_strain_rate; 
  // SymmetricTensor<2,6> isotropic_tensor; 
  const double F_iso = 0.5; const double G_iso = 0.5; const double H_iso = 0.5; 
  const double L_iso = 1.5; const double M_iso = 1.5; const double N_iso = 1.5; 

  // for (unsigned int i=0; i < 3; ++i)
  //         {
  //           test_stress[i][i]         = i + 1.0; 
  //           test_stress[ji[i]][ki[i]] = i + 4.0;
  //           // test_stress[ki[i]][ji[i]] = i + 4.0;

  //           test_strain_rate[i][i]          = i + 1.5; 
  //           test_strain_rate[ji[i]][ki[i]]  = i + 4.5; 
  //           // test_strain_rate[ki[i]][ji[i]]  = i + 4.5;

  //           isotropic_tensor[i][i]          = 2.0/3.0;
  //           isotropic_tensor[ji[i]][ki[i]]  = -1.0/3.0;
  //           isotropic_tensor[i+3][i+3]      = 1.0;
  //         }; 

  aspect::MaterialModel::MaterialUtilities::OrthotropicRheology<3> ortho_rheo_3d;

  const SymmetricTensor<2,6> iso_fluidity_tensor = ortho_rheo_3d.fluidity_tensor_cpo_frame(F_iso, G_iso, H_iso, L_iso, M_iso, N_iso);
  // const SymmetricTensor<2,6> iso_viscosity_tensor = ortho_rheo_3d.viscosity_tensor_cpo_frame(F_iso, G_iso, H_iso, L_iso, M_iso, N_iso);
  
  // //compare_kelvin_sym_tensors_approx(iso_fluidity_tensor, isotropic_tensor);
  // //compare_kelvin_sym_tensors_approx(iso_viscosity_tensor, isotropic_tensor); 

  // // second invariant as defined in composite
  // const double isotropic_strain_rate_invariant = std::sqrt(std::max(-aspect::Utilities::Tensors::consistent_second_invariant_of_deviatoric_tensor(aspect::Utilities::Tensors::consistent_deviator(test_strain_rate)), 0.)); 

  // CHECK(isotropic_strain_rate_invariant > 0);

  // const double orthotropic_strain_rate_inv = ortho_rheo_3d.strain_rate_invariant(test_strain_rate, F_iso, G_iso, H_iso, L_iso, M_iso, N_iso);
  // const double orthotropic_stess_inv = ortho_rheo_3d.stress_invariant(test_stress, F_iso, G_iso, H_iso, L_iso, M_iso, N_iso);
  
  // CHECK( orthotropic_strain_rate_inv == Approx(isotropic_strain_rate_invariant));
}

















