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

/**
 * Compare two symmetric tensors
 */
inline void compare_tensors_approx(
    const dealii::Tensor<2,3> &computed,
    const dealii::Tensor<2,3> &expected)
{
  INFO("tensors are not the same: \n" <<
       "expected = " << expected[0][0] << "\n" <<
       "           " << expected[0][1] << " " << expected[1][1] << "\n" <<
       "           " << expected[0][2] << " " << expected[1][2] << " " << expected[2][2] << "\n" <<
       "computed = " << computed[0][0] << "\n" <<
       "           " << computed[0][1] << " " << computed[1][1] << "\n" <<
       "           " << computed[0][2] << " " << computed[1][2] << " " << computed[2][2] << "\n" );
  const double tol = 1e-14;
  CHECK(
        ((computed[0][0] == Approx(expected[0][0]) || std::fabs(computed[0][0]) < tol)
          && (computed[0][1] == Approx(expected[0][1]) || std::fabs(computed[0][1]) < tol)
          && (computed[0][2] == Approx(expected[0][2]) || std::fabs(computed[0][2]) < tol)
          && (computed[1][0] == Approx(expected[1][0]) || std::fabs(computed[1][0]) < tol)
          && (computed[1][1] == Approx(expected[1][1]) || std::fabs(computed[1][1]) < tol)
          && (computed[1][2] == Approx(expected[1][2]) || std::fabs(computed[1][2]) < tol)
          && (computed[2][0] == Approx(expected[2][0]) || std::fabs(computed[2][0]) < tol)
          && (computed[2][1] == Approx(expected[2][1]) || std::fabs(computed[2][1]) < tol)
          && (computed[2][2] == Approx(expected[2][2]) || std::fabs(computed[2][2]) < tol))
        );
}

/**
 * Compare two symmetric rank 4 tensors in kelvin notation
 */
inline void compare_kelvin_sym_tensors_approx(
    const dealii::SymmetricTensor<2,6> &computed,
    const dealii::SymmetricTensor<2,6> &expected)
{
  INFO("tensors are not the same: \n" <<
       "expected = " << expected[0][0] << "\n" <<
       "           " << expected[1][0] << " " << expected[1][1] <<"\n" <<
       "           " << expected[2][0] << " " << expected[2][1] << " " << expected[2][2] << "\n" <<
       "           " << expected[3][0] << " " << expected[3][1] << " " << expected[3][2] << " " << expected[3][3] << "\n" <<
       "           " << expected[4][0] << " " << expected[4][1] << " " << expected[4][2] << " " << expected[4][3] << " " << expected[4][4] << "\n" <<
       "           " << expected[5][0] << " " << expected[5][1] << " " << expected[5][2] << " " << expected[5][3] << " " << expected[5][4] << " " << expected[5][5] << "\n" <<

       "computed = " << computed[0][0] << "\n" <<
       "           " << computed[1][0] << " " << computed[1][1] <<"\n" <<
       "           " << computed[2][0] << " " << computed[2][1] << " " << computed[2][2] << "\n" <<
       "           " << computed[3][0] << " " << computed[3][1] << " " << computed[3][2] << " " << computed[3][3] << "\n" <<
       "           " << computed[4][0] << " " << computed[4][1] << " " << computed[4][2] << " " << computed[4][3] << " " << computed[4][4] << "\n" <<
       "           " << computed[5][0] << " " << computed[5][1] << " " << computed[5][2] << " " << computed[5][3] << " " << computed[5][4] << " " << computed[5][5] << "\n" );

  const double tol = 1e-14;
  CHECK(
        ((    computed[0][0] == Approx(expected[0][0]) || std::fabs(computed[0][0]) < tol)
          && (computed[0][1] == Approx(expected[0][1]) || std::fabs(computed[0][1]) < tol)
          && (computed[0][2] == Approx(expected[0][2]) || std::fabs(computed[0][2]) < tol)
          && (computed[0][3] == Approx(expected[0][3]) || std::fabs(computed[0][3]) < tol)
          && (computed[0][4] == Approx(expected[0][4]) || std::fabs(computed[0][4]) < tol)
          && (computed[0][5] == Approx(expected[0][5]) || std::fabs(computed[0][5]) < tol)
          
          && (computed[1][1] == Approx(expected[1][1]) || std::fabs(computed[1][1]) < tol)
          && (computed[1][2] == Approx(expected[1][2]) || std::fabs(computed[1][2]) < tol)
          && (computed[1][3] == Approx(expected[1][3]) || std::fabs(computed[1][3]) < tol)
          && (computed[1][4] == Approx(expected[1][4]) || std::fabs(computed[1][4]) < tol)
          && (computed[1][5] == Approx(expected[1][5]) || std::fabs(computed[1][5]) < tol)
          
          && (computed[2][2] == Approx(expected[2][2]) || std::fabs(computed[2][2]) < tol)
          && (computed[2][3] == Approx(expected[2][3]) || std::fabs(computed[2][3]) < tol)
          && (computed[2][4] == Approx(expected[2][4]) || std::fabs(computed[2][4]) < tol)
          && (computed[2][5] == Approx(expected[2][5]) || std::fabs(computed[2][5]) < tol)
          
          && (computed[3][3] == Approx(expected[3][3]) || std::fabs(computed[3][3]) < tol)
          && (computed[3][4] == Approx(expected[3][4]) || std::fabs(computed[3][4]) < tol)
          && (computed[3][5] == Approx(expected[3][5]) || std::fabs(computed[3][5]) < tol)
          
          && (computed[4][4] == Approx(expected[4][4]) || std::fabs(computed[4][4]) < tol)
          && (computed[4][5] == Approx(expected[4][5]) || std::fabs(computed[4][5]) < tol)
          
          && (computed[5][5] == Approx(expected[5][5]) || std::fabs(computed[5][5]) < tol)
        )
          
      );
}


TEST_CASE("OrthotropicRheology Functions in 3D")
{
  using namespace aspect;
  using namespace dealii; // for symmetric tensors
  
  aspect::MaterialModel::MaterialUtilities::OrthotropicRheology<3> ortho_rheo_3d;


  std::vector<int> ji = {1,2,0}; // tuple of indices shifted by one
  std::vector<int> ki = {2,0,1}; // tuple of indices shifted by two


  const double test_prefactor = 1.157e6;
  const double test_stress_exponent  = 3.76;

  SymmetricTensor<2,3> test_stress; 
  SymmetricTensor<2,3> test_strain_rate; 
  SymmetricTensor<2,6> isotropic_tensor; 
  const double F_iso = 0.5; const double G_iso = 0.5; const double H_iso = 0.5; 
  const double L_iso = 1.5; const double M_iso = 1.5; const double N_iso = 1.5; 

  const double F_rand = 0.3; const double G_rand = 0.654; const double H_rand = 1.5 - F_rand - G_rand; 
  const double L_rand = 2; const double M_rand = 1.2; const double N_rand = 1.6; 

  for (unsigned int i=0; i < 3; ++i)
          {
            test_stress[i][i]         = i + 1.0; 
            test_stress[ji[i]][ki[i]] = i + 4.0;
            // test_stress[ki[i]][ji[i]] = i + 4.0;

            test_strain_rate[i][i]          = i + 1.5; 
            test_strain_rate[ji[i]][ki[i]]  = i + 4.5; 
            // test_strain_rate[ki[i]][ji[i]]  = i + 4.5;

            isotropic_tensor[i][i]          = 2.0/3.0;
            isotropic_tensor[ji[i]][ki[i]]  = -1.0/3.0;
            isotropic_tensor[i+3][i+3]      = 1.0;
          }; 

  const SymmetricTensor<2,3> test_stress_dev = aspect::Utilities::Tensors::consistent_deviator(test_stress);
  const SymmetricTensor<2,3> test_strain_rate_dev = aspect::Utilities::Tensors::consistent_deviator(test_strain_rate);

  // second invariant as defined in composite rheology
  const double isotropic_strain_rate_invariant = std::sqrt(std::max( -aspect::Utilities::Tensors::consistent_second_invariant_of_deviatoric_tensor(test_strain_rate_dev), 0.)); 
  
  const double isotropic_stress_invariant = std::sqrt(std::max(-aspect::Utilities::Tensors::consistent_second_invariant_of_deviatoric_tensor(test_stress_dev), 0.)); 


  INFO("isotropic invariant <= 0");
  CHECK(isotropic_strain_rate_invariant > 0);

  /*
  * Test cases in cpo reference frame
  */

  // Testing consistency with isotropic case
  SymmetricTensor<2,6> fluidity_tensor = ortho_rheo_3d.fluidity_tensor_cpo_frame(F_iso, G_iso, H_iso, L_iso, M_iso, N_iso);
  SymmetricTensor<2,6> viscosity_tensor = ortho_rheo_3d.viscosity_tensor_cpo_frame(F_iso, G_iso, H_iso, L_iso, M_iso, N_iso);
  
  INFO("isotropic tensor is not the same as isotropic tensor from orthotropic rheology");
  compare_kelvin_sym_tensors_approx(fluidity_tensor, isotropic_tensor);
  compare_kelvin_sym_tensors_approx(viscosity_tensor, isotropic_tensor); 

  // check conversion to rank 4 tensor and consistency with deviator
  SymmetricTensor<4,3> full_fluidity_tensor = ortho_rheo_3d.kelvin_to_r4_tensor(fluidity_tensor); 
  compare_tensors_approx(full_fluidity_tensor*test_stress, test_stress_dev);

  SymmetricTensor<4,3> full_viscosity_tensor = ortho_rheo_3d.kelvin_to_r4_tensor(viscosity_tensor);  
  compare_tensors_approx(full_viscosity_tensor*test_strain_rate, test_strain_rate_dev);


  double orthotropic_strain_rate_inv = ortho_rheo_3d.strain_rate_invariant(test_strain_rate, F_iso, G_iso, H_iso, L_iso, M_iso, N_iso);
  double orthotropic_stress_inv = ortho_rheo_3d.stress_invariant(test_stress, F_iso, G_iso, H_iso, L_iso, M_iso, N_iso);
  
  INFO("isotropic invariant and orthotropic invariant for isotropic case are not the same");
  CHECK( orthotropic_strain_rate_inv == Approx(isotropic_strain_rate_invariant));
  CHECK( orthotropic_stress_inv == Approx(isotropic_stress_invariant));

  /*
  * Consistency of anisotropic rheology in CPO reference frame
  */

  const double orthotropic_strain_rate_inv_test = ortho_rheo_3d.strain_rate_invariant(test_strain_rate, F_rand, G_rand, H_rand, L_rand, M_rand, N_rand);
  const double orthotropic_stress_inv_test = ortho_rheo_3d.stress_invariant(test_stress, F_rand, G_rand, H_rand, L_rand, M_rand, N_rand);

  fluidity_tensor = ortho_rheo_3d.fluidity_tensor_cpo_frame(F_rand, G_rand, H_rand, L_rand, M_rand, N_rand);
  viscosity_tensor = ortho_rheo_3d.viscosity_tensor_cpo_frame(F_rand, G_rand, H_rand, L_rand, M_rand, N_rand);

  // rank 4 symmetric tensor still in cpo frame
  full_fluidity_tensor = ortho_rheo_3d.kelvin_to_r4_tensor(fluidity_tensor);
  full_viscosity_tensor = ortho_rheo_3d.kelvin_to_r4_tensor(viscosity_tensor);  

  // forward rheology then backward rheology
  SymmetricTensor<2,3> orthotropic_strain_rate = test_prefactor*std::pow(orthotropic_stress_inv_test, test_stress_exponent-1)*full_fluidity_tensor*test_stress;

  orthotropic_strain_rate_inv = ortho_rheo_3d.strain_rate_invariant(orthotropic_strain_rate, F_rand, G_rand, H_rand, L_rand, M_rand, N_rand);

  SymmetricTensor<2,3> orthotropic_stress = (std::pow(test_prefactor, -1/test_stress_exponent)*std::pow(orthotropic_strain_rate_inv, (1-test_stress_exponent)/test_stress_exponent)*full_viscosity_tensor*orthotropic_strain_rate);

  INFO("forward rheology then backward rheology not consistent");
  compare_tensors_approx(orthotropic_stress, test_stress_dev);

  // backward rheology then forward rheology 

  orthotropic_stress = (std::pow(test_prefactor, -1/test_stress_exponent)*std::pow(orthotropic_strain_rate_inv_test, (1-test_stress_exponent)/test_stress_exponent)*full_viscosity_tensor*test_strain_rate);
  
  orthotropic_stress_inv = ortho_rheo_3d.stress_invariant(orthotropic_stress, F_rand, G_rand, H_rand, L_rand, M_rand, N_rand);

  orthotropic_strain_rate = test_prefactor*std::pow(orthotropic_stress_inv, test_stress_exponent-1)*full_fluidity_tensor*orthotropic_stress; 

  INFO("backward rheology then forward rheology not consistent");
  compare_tensors_approx(orthotropic_strain_rate, test_strain_rate_dev);

  /*
  * checking consistency in the Lab frame
  */

  // define random CPO reference frame with passive rotation matrix (what one gets from DREX)
  std::vector<double> EA = {50.0, 15.0, 40.0};
  Tensor<2,3> rot_to_cpo   = aspect::Utilities::zxz_euler_angles_to_rotation_matrix(EA[0], EA[1], EA[2]);
  // rot_from_cpo = rot_to_cpo^-1 = transpose(rot_to_cpo)

  // move test stress and strain into cpo frame 
  SymmetricTensor<2,3> test_strain_rate_cpo = symmetrize(rot_to_cpo*test_strain_rate*transpose(rot_to_cpo));
  SymmetricTensor<2,3> test_stress_cpo = symmetrize(rot_to_cpo*test_stress*transpose(rot_to_cpo));
  
  /*
  * check rotation of viscosity tensor
  */ 

  // tensor part of correct rotation in lab frame
  SymmetricTensor<2,3> orthotropic_stress_tensor_part = symmetrize(transpose(rot_to_cpo)*(full_viscosity_tensor*test_strain_rate_cpo)*rot_to_cpo);

  // rotation using full stiffness_matrix
  SymmetricTensor<4,3> full_viscosity_tensor_lab = aspect::Utilities::Tensors::rotate_full_stiffness_tensor(transpose(rot_to_cpo), full_viscosity_tensor);

  INFO("compare tensor part of orthotropic stress to check rotation of viscosity tensor")
  compare_tensors_approx(full_viscosity_tensor_lab*test_strain_rate, orthotropic_stress_tensor_part);
  
  // rotation using rotation in kelvin notation (computationally more efficient)
  // not yet implemented
  // SymmetricTensor<2,6> viscosity_tensor_lab = aspect::Utilities::rotate_kelvin_tensor(transpose(rot_to_cpo), viscosity_tensor)
  
  // full_viscosity_tensor_lab = ortho_rheo_3d.kelvin_to_r4_tensor(viscosity_tensor_lab)

  /*
  * backwards then forward computation (example usage)
  */

  // computing the invariant, notice that strain_rate is in cpo frame
  orthotropic_strain_rate_inv = ortho_rheo_3d.strain_rate_invariant(test_strain_rate_cpo, F_rand, G_rand, H_rand, L_rand, M_rand, N_rand);

  // computing anisotropic stress
  orthotropic_stress = ( std::pow(test_prefactor, -1/test_stress_exponent)*std::pow(orthotropic_strain_rate_inv, (1-test_stress_exponent)/test_stress_exponent)*full_viscosity_tensor_lab*test_strain_rate);
  
  // computing stress invariant
  SymmetricTensor<2,3> orthotropic_stress_cpo = symmetrize(rot_to_cpo*orthotropic_stress*transpose(rot_to_cpo));

  orthotropic_stress_inv = ortho_rheo_3d.stress_invariant(orthotropic_stress_cpo, F_rand, G_rand, H_rand, L_rand, M_rand, N_rand);

  SymmetricTensor<4,3> full_fluidity_tensor_lab = aspect::Utilities::Tensors::rotate_full_stiffness_tensor(transpose(rot_to_cpo), full_fluidity_tensor);
  
  // computing orthotropic strain_rate 
  orthotropic_strain_rate = test_prefactor*std::pow( orthotropic_stress_inv, test_stress_exponent-1 )*full_fluidity_tensor_lab*orthotropic_stress ;
  
  INFO("check backward forward rheology in lab frame");
  compare_tensors_approx(orthotropic_strain_rate, test_strain_rate_dev) ;

  /*
  * forward then backward computation
  */
  
  // computing the invariant, notice that strain_rate is in cpo frame
  orthotropic_stress_inv = ortho_rheo_3d.stress_invariant(test_stress_cpo, F_rand, G_rand, H_rand, L_rand, M_rand, N_rand);

  // computing anisotropic strain:rate
  orthotropic_strain_rate = test_prefactor*std::pow( orthotropic_stress_inv, test_stress_exponent-1 )*full_fluidity_tensor_lab*test_stress ;
  
  // computing strain_rate invariant
  SymmetricTensor<2,3> orthotropic_strain_rate_cpo = symmetrize(rot_to_cpo*orthotropic_strain_rate*transpose(rot_to_cpo));

  orthotropic_strain_rate_inv = ortho_rheo_3d.strain_rate_invariant(orthotropic_strain_rate_cpo, F_rand, G_rand, H_rand, L_rand, M_rand, N_rand);
  
  // computing orthotropic stress
  orthotropic_stress = ( std::pow(test_prefactor, -1/test_stress_exponent)*std::pow(orthotropic_strain_rate_inv, (1-test_stress_exponent)/test_stress_exponent)*full_viscosity_tensor_lab*orthotropic_strain_rate);
  
  INFO("check forward backward rheology in lab frame");
  compare_tensors_approx(orthotropic_stress, test_stress_dev) ;

}

















