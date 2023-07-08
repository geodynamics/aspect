/*
  Copyright (C) 2022 - 2023 by the authors of the ASPECT code.

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
#include <aspect/particle/property/crystal_preferred_orientation.h>
#include <aspect/particle/property/cpo_elastic_tensor.h>
#include <deal.II/base/parameter_handler.h>


/**
 * Compare the given two std::array<double,3> entries with an epsilon (using Catch::Approx)
 */
inline void compare_3d_arrays_approx(
  const std::array<double,3> &computed,
  const std::array<double,3> &expected)
{
  CHECK(computed.size() == expected.size());
  for (unsigned int i=0; i< computed.size(); ++i)
    {
      INFO("vector index i=" << i << ": ");
      CHECK(computed[i] == Approx(expected[i]));
    }
}


/**
 * Compare the given two std::array<double,3> entries with an epsilon (using Catch::Approx)
 */
inline void compare_3d_arrays_approx(
  const std::vector<double> &computed,
  const std::vector<double> &expected)
{
  CHECK(computed.size() == expected.size());
  for (unsigned int i=0; i< computed.size(); ++i)
    {
      INFO("vector index i=" << i << ": ");
      CHECK(computed[i] == Approx(expected[i]));
    }
}


/**
 * Compare two rotation matrices
 */
inline void compare_rotation_matrices_approx(
  const std::array<std::array<double,3>,3> &computed,
  const std::array<std::array<double,3>,3> &expected)
{
  // sign of eigenvector is not important
  INFO("rotation matrices are not the same: \n" <<
       "expected = " << expected[0][0] << " " << expected[0][1] << " " << expected[0][2] << "\n" <<
       "           " << expected[1][0] << " " << expected[1][1] << " " << expected[1][2] << "\n" <<
       "           " << expected[2][0] << " " << expected[2][1] << " " << expected[2][2] << "\n" <<
       "computed = " << computed[0][0] << " " << computed[0][1] << " " << computed[0][2] << "\n" <<
       "           " << computed[1][0] << " " << computed[1][1] << " " << computed[1][2] << "\n" <<
       "           " << computed[2][0] << " " << computed[2][1] << " " << computed[2][2] << "\n" );
  CHECK((
          (computed[0][0] == Approx(expected[0][0]) && computed[0][1] == Approx(expected[0][1]) && computed[0][2] == Approx(expected[0][2]) &&
           computed[1][0] == Approx(expected[1][0]) && computed[1][1] == Approx(expected[1][1]) && computed[1][2] == Approx(expected[1][2]) &&
           computed[2][0] == Approx(expected[2][0]) && computed[2][1] == Approx(expected[2][1]) && computed[2][2] == Approx(expected[2][2]))
          ||
          (computed[0][0] == Approx(-expected[0][0]) && computed[0][1] == Approx(-expected[0][1]) && computed[0][2] == Approx(-expected[0][2]) &&
           computed[1][0] == Approx(-expected[1][0]) && computed[1][1] == Approx(-expected[1][1]) && computed[1][2] == Approx(-expected[1][2]) &&
           computed[2][0] == Approx(-expected[2][0]) && computed[2][1] == Approx(-expected[2][1]) && computed[2][2] == Approx(-expected[2][2]))));
}

/**
 * Compare two rotation matrices
 */
inline void compare_rotation_matrices_approx(
  const dealii::Tensor<2,3> &computed,
  const dealii::Tensor<2,3> &expected)
{
  // sign of eigenvector is not important
  INFO("rotation matrices are not the same: \n" <<
       "expected = " << expected[0][0] << " " << expected[0][1] << " " << expected[0][2] << "\n" <<
       "           " << expected[1][0] << " " << expected[1][1] << " " << expected[1][2] << "\n" <<
       "           " << expected[2][0] << " " << expected[2][1] << " " << expected[2][2] << "\n" <<
       "computed = " << computed[0][0] << " " << computed[0][1] << " " << computed[0][2] << "\n" <<
       "           " << computed[1][0] << " " << computed[1][1] << " " << computed[1][2] << "\n" <<
       "           " << computed[2][0] << " " << computed[2][1] << " " << computed[2][2] << "\n" );
  const double tol = 1e-14;
  CHECK((
          ((computed[0][0] == Approx(expected[0][0]) || std::fabs(computed[0][0]) < tol)
           && (computed[0][1] == Approx(expected[0][1]) || std::fabs(computed[0][1]) < tol)
           && (computed[0][2] == Approx(expected[0][2]) || std::fabs(computed[0][2]) < tol)
           && (computed[1][0] == Approx(expected[1][0]) || std::fabs(computed[1][0]) < tol)
           && (computed[1][1] == Approx(expected[1][1]) || std::fabs(computed[1][1]) < tol)
           && (computed[1][2] == Approx(expected[1][2]) || std::fabs(computed[1][2]) < tol)
           && (computed[2][0] == Approx(expected[2][0]) || std::fabs(computed[2][0]) < tol)
           && (computed[2][1] == Approx(expected[2][1]) || std::fabs(computed[2][1]) < tol)
           && (computed[2][2] == Approx(expected[2][2]) || std::fabs(computed[2][2]) < tol))
          ||
          ((computed[0][0] == Approx(-expected[0][0]) || std::fabs(computed[0][0]) < tol)
           && (computed[0][1] == Approx(-expected[0][1]) || std::fabs(computed[0][1]) < tol)
           && (computed[0][2] == Approx(-expected[0][2]) || std::fabs(computed[0][2]) < tol)
           && (computed[1][0] == Approx(-expected[1][0]) || std::fabs(computed[1][0]) < tol)
           && (computed[1][1] == Approx(-expected[1][1]) || std::fabs(computed[1][1]) < tol)
           && (computed[1][2] == Approx(-expected[1][2]) || std::fabs(computed[1][2]) < tol)
           && (computed[2][0] == Approx(-expected[2][0]) || std::fabs(computed[2][0]) < tol)
           && (computed[2][1] == Approx(-expected[2][1]) || std::fabs(computed[2][1]) < tol)
           && (computed[2][2] == Approx(-expected[2][2]) || std::fabs(computed[2][2]) < tol))));
}

TEST_CASE("CPO core: Store and Load")
{
  // test store and load functions
  // the first and last element should not be changed, because we start at
  // data position = 1.
  aspect::Particle::Property::CrystalPreferredOrientation<3> cpo;
  aspect::ParameterHandler prm;
  cpo.declare_parameters(prm);
  prm.enter_subsection("Postprocess");
  {
    prm.enter_subsection("Particles");
    {
      prm.enter_subsection("Crystal Preferred Orientation");
      {
        prm.set("Random number seed","1");
        prm.set("Number of grains per particle","3");
        prm.set("CPO derivatives algorithm","Spin tensor");
        prm.set("Property advection method","Backward Euler");
        prm.enter_subsection("Initial grains");
        {
          prm.set("Model name","Uniform grains and random uniform rotations");
          // Let the minerals just passively rotate with the rotation of
          // the particle caused by the flow.
          prm.set("Minerals","Passive,Passive");
          prm.set("Volume fractions minerals","0.7,0.3");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection ();
    }
    prm.leave_subsection ();
  }
  prm.leave_subsection ();

  cpo.parse_parameters(prm);
  cpo.initialize();

  unsigned int cpo_data_position = 1;
  std::vector<double> data_array(70,-1.);
  dealii::ArrayView<double> data(&data_array[0],70);

  cpo.set_volume_fraction_mineral(cpo_data_position,data,0,0.7);
  cpo.set_volume_fraction_mineral(cpo_data_position,data,1,0.3);

  using namespace dealii;
  Tensor<2,3> rotation_matrix;
  rotation_matrix[TableIndices<2>(0,0)] = 0;
  rotation_matrix[TableIndices<2>(0,1)] = 1./1000.;
  rotation_matrix[TableIndices<2>(0,2)] = 2./1000.;
  rotation_matrix[TableIndices<2>(1,0)] = 3./1000.;
  rotation_matrix[TableIndices<2>(1,1)] = 4./1000.;
  rotation_matrix[TableIndices<2>(1,2)] = 5./1000.;
  rotation_matrix[TableIndices<2>(2,0)] = 6./1000.;
  rotation_matrix[TableIndices<2>(2,1)] = 7./1000.;
  rotation_matrix[TableIndices<2>(2,2)] = 8./1000.;
  cpo.set_rotation_matrix_grains(cpo_data_position,data,0,0,rotation_matrix);

  rotation_matrix[TableIndices<2>(0,0)] = 9./1000.;
  rotation_matrix[TableIndices<2>(0,1)] = 10./1000.;
  rotation_matrix[TableIndices<2>(0,2)] = 11./1000.;
  rotation_matrix[TableIndices<2>(1,0)] = 12./1000.;
  rotation_matrix[TableIndices<2>(1,1)] = 13./1000.;
  rotation_matrix[TableIndices<2>(1,2)] = 14./1000.;
  rotation_matrix[TableIndices<2>(2,0)] = 15./1000.;
  rotation_matrix[TableIndices<2>(2,1)] = 16./1000.;
  rotation_matrix[TableIndices<2>(2,2)] = 17./1000.;
  cpo.set_rotation_matrix_grains(cpo_data_position,data,0,1,rotation_matrix);

  rotation_matrix[TableIndices<2>(0,0)] = 18./1000.;
  rotation_matrix[TableIndices<2>(0,1)] = 19./1000.;
  rotation_matrix[TableIndices<2>(0,2)] = 20./1000.;
  rotation_matrix[TableIndices<2>(1,0)] = 21./1000.;
  rotation_matrix[TableIndices<2>(1,1)] = 22./1000.;
  rotation_matrix[TableIndices<2>(1,2)] = 23./1000.;
  rotation_matrix[TableIndices<2>(2,0)] = 24./1000.;
  rotation_matrix[TableIndices<2>(2,1)] = 25./1000.;
  rotation_matrix[TableIndices<2>(2,2)] = 26./1000.;
  cpo.set_rotation_matrix_grains(cpo_data_position,data,0,2,rotation_matrix);

  cpo.set_volume_fractions_grains(cpo_data_position,data,0,0,0.1);
  cpo.set_volume_fractions_grains(cpo_data_position,data,0,1,0.2);
  cpo.set_volume_fractions_grains(cpo_data_position,data,0,2,0.3);


  rotation_matrix[TableIndices<2>(0,0)] = 27./1000.;
  rotation_matrix[TableIndices<2>(0,1)] = 28./1000.;
  rotation_matrix[TableIndices<2>(0,2)] = 29./1000.;
  rotation_matrix[TableIndices<2>(1,0)] = 30./1000.;
  rotation_matrix[TableIndices<2>(1,1)] = 31./1000.;
  rotation_matrix[TableIndices<2>(1,2)] = 32./1000.;
  rotation_matrix[TableIndices<2>(2,0)] = 33./1000.;
  rotation_matrix[TableIndices<2>(2,1)] = 34./1000.;
  rotation_matrix[TableIndices<2>(2,2)] = 35./1000.;
  cpo.set_rotation_matrix_grains(cpo_data_position,data,1,0,rotation_matrix);

  rotation_matrix[TableIndices<2>(0,0)] = 36./1000.;
  rotation_matrix[TableIndices<2>(0,1)] = 37./1000.;
  rotation_matrix[TableIndices<2>(0,2)] = 38./1000.;
  rotation_matrix[TableIndices<2>(1,0)] = 39./1000.;
  rotation_matrix[TableIndices<2>(1,1)] = 40./1000.;
  rotation_matrix[TableIndices<2>(1,2)] = 41./1000.;
  rotation_matrix[TableIndices<2>(2,0)] = 42./1000.;
  rotation_matrix[TableIndices<2>(2,1)] = 43./1000.;
  rotation_matrix[TableIndices<2>(2,2)] = 44./1000.;
  cpo.set_rotation_matrix_grains(cpo_data_position,data,1,1,rotation_matrix);

  rotation_matrix[TableIndices<2>(0,0)] = 45./1000.;
  rotation_matrix[TableIndices<2>(0,1)] = 46./1000.;
  rotation_matrix[TableIndices<2>(0,2)] = 47./1000.;
  rotation_matrix[TableIndices<2>(1,0)] = 48./1000.;
  rotation_matrix[TableIndices<2>(1,1)] = 49./1000.;
  rotation_matrix[TableIndices<2>(1,2)] = 50./1000.;
  rotation_matrix[TableIndices<2>(2,0)] = 51./1000.;
  rotation_matrix[TableIndices<2>(2,1)] = 52./1000.;
  rotation_matrix[TableIndices<2>(2,2)] = 53./1000.;
  cpo.set_rotation_matrix_grains(cpo_data_position,data,1,2,rotation_matrix);

  cpo.set_volume_fractions_grains(cpo_data_position,data,1,0,0.4);
  cpo.set_volume_fractions_grains(cpo_data_position,data,1,1,0.5);
  cpo.set_volume_fractions_grains(cpo_data_position,data,1,2,0.6);

  data[0] = 20847932.2;
  data[65] = 6541684.3;

  cpo.set_deformation_type(cpo_data_position,data,0,(double)aspect::Particle::Property::DeformationType::passive);
  cpo.set_deformation_type(cpo_data_position,data,1,(double)aspect::Particle::Property::DeformationType::passive);


  CHECK(data[0] == Approx(20847932.2)); // before data position
  // mineral 1
  CHECK(data[1] ==  Approx(0.0)); // deformation type
  CHECK(data[2] ==  Approx(0.7)); //Mineral volume fraction
  CHECK(data[3] ==  Approx(0.1)); // grain 0 volume fraction
  double counter_rotation = 0.;
  for (size_t iii = 4; iii < 13; iii++)
    {
      CHECK(data[iii] == Approx(counter_rotation));
      counter_rotation += 1./1000.;
    }
  CHECK(data[13] ==  Approx(0.2)); // grain 1 volume fraction
  for (size_t iii = 14; iii < 23; iii++)
    {
      CHECK(data[iii] == Approx(counter_rotation));
      counter_rotation += 1./1000.;
    }
  CHECK(data[23] ==  Approx(0.3)); // grain 2 volume fraction
  for (size_t iii = 24; iii < 33; iii++)
    {
      CHECK(data[iii] == Approx(counter_rotation));
      counter_rotation += 1./1000.;
    }
  // mineral 2
  CHECK(data[33] ==  Approx(0.0)); // deformation type
  CHECK(data[34] ==  Approx(0.3)); //Mineral volume fraction
  CHECK(data[35] ==  Approx(0.4)); // grain 0 volume fraction
  for (size_t iii = 36; iii < 45; iii++)
    {
      CHECK(data[iii] == Approx(counter_rotation));
      counter_rotation += 1./1000.;
    }
  CHECK(data[45] ==  Approx(0.5)); // grain 1 volume fraction
  for (size_t iii = 46; iii < 55; iii++)
    {
      CHECK(data[iii] == Approx(counter_rotation));
      counter_rotation += 1./1000.;
    }
  CHECK(data[55] ==  Approx(0.6)); // grain 2 volume fraction
  for (size_t iii = 56; iii < 65; iii++)
    {
      CHECK(data[iii] == Approx(counter_rotation));
      counter_rotation += 1./1000.;
    }
  CHECK(data[65] == Approx(6541684.3)); // after data position
  CHECK(data[66] == Approx(-1.0)); // after data position

}

TEST_CASE("CPO core: Spin tensor")
{
  using namespace dealii;
  using namespace aspect;

  {
    // test initialization in 3d.

    aspect::Particle::Property::CrystalPreferredOrientation<3> cpo_3d;
    aspect::ParameterHandler prm;
    cpo_3d.declare_parameters(prm);
    prm.enter_subsection("Postprocess");
    {
      prm.enter_subsection("Particles");
      {
        prm.enter_subsection("Crystal Preferred Orientation");
        {
          prm.set("Random number seed","1");
          prm.set("Number of grains per particle","5");
          prm.set("CPO derivatives algorithm","Spin tensor");
          prm.set("Property advection method","Forward Euler");
          prm.enter_subsection("Initial grains");
          {
            prm.set("Model name","Uniform grains and random uniform rotations");
            // Let the minerals just passively rotate with the rotation of
            // the particle caused by the flow.
            prm.set("Minerals","Passive,Passive");
            prm.set("Volume fractions minerals","0.5,0.5");
          }
          prm.leave_subsection();
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }
    prm.leave_subsection ();

    cpo_3d.parse_parameters(prm);
    cpo_3d.initialize();

    Point<3> dummy_point;
    std::vector<double> data;
    cpo_3d.initialize_one_particle_property(dummy_point, data);
    // The CPO particles are initialized. With the same seed, the outcome should
    // always be the same, so test that for seed = 1. Furthermore, in the data
    // I can only really test that the first entry is the water content (0) and
    // that every first entry of each particle is 1/n_grains = 1/10 = 0.1.
    CHECK(data[0] == Approx(-1.0));//isnan(data[0])); // default fabric type which is only computed on a update
    CHECK(data[1] == Approx(0.5)); // default volume fraction olivine
    CHECK(data[2] == Approx(0.2));
    CHECK(data[3] == Approx(0.8794381492));
    CHECK(data[4] == Approx(-0.4757482965));
    CHECK(data[5] == Approx(-0.015877661));
    CHECK(data[6] == Approx(-0.197154988));
    CHECK(data[7] == Approx(-0.394402755));
    CHECK(data[8] == Approx(0.8975390674));
    CHECK(data[9] == Approx(-0.4332648756));
    CHECK(data[10] == Approx(-0.7861997362));
    CHECK(data[11] == Approx(-0.4406489786));
    CHECK(data[12] == Approx(0.2));
    CHECK(data[13] == Approx(-0.911475687));
    CHECK(data[14] == Approx(0.4113536692));
    CHECK(data[15] == Approx(0.0004804047));
    CHECK(data[16] == Approx(-0.3056521145));
    CHECK(data[17] == Approx(-0.6780433033));
    CHECK(data[18] == Approx(0.6684564786));
    CHECK(data[19] == Approx(0.2752977604));
    CHECK(data[20] == Approx(0.6091349914));
    CHECK(data[21] == Approx(0.7437511045));
    CHECK(data[22] == Approx(0.2));
    CHECK(data[32] == Approx(0.2));
    CHECK(data[42] == Approx(0.2));
    CHECK(data[52] == Approx(-1.0));// default fabric type which is only computed on a update
    CHECK(data[53] == Approx(0.5)); // default volume fraction olivine
    CHECK(data[54] == Approx(0.2));
    CHECK(data[64] == Approx(0.2));
    CHECK(data[74] == Approx(0.2));
    CHECK(data[84] == Approx(0.2));
    CHECK(data[94] == Approx(0.2));
    CHECK(data[103] == Approx(-0.0688278144));

    std::vector<double> volume_fractions(5,0.2);
    std::vector<dealii::Tensor<2,3>> rotation_matrices_minerals(5);
    rotation_matrices_minerals[0][0][0] = 0.5;
    rotation_matrices_minerals[0][0][1] = 0.5;
    rotation_matrices_minerals[0][0][2] = 0.5;
    rotation_matrices_minerals[0][1][0] = 0.5;
    rotation_matrices_minerals[0][1][1] = 0.5;
    rotation_matrices_minerals[0][1][2] = 0.5;
    rotation_matrices_minerals[0][2][0] = 0.5;
    rotation_matrices_minerals[0][2][1] = 0.5;
    rotation_matrices_minerals[0][2][2] = 0.5;

    rotation_matrices_minerals[1][0][0] = 0.1;
    rotation_matrices_minerals[1][0][1] = 0.2;
    rotation_matrices_minerals[1][0][2] = 0.3;
    rotation_matrices_minerals[1][1][0] = 0.4;
    rotation_matrices_minerals[1][1][1] = 0.5;
    rotation_matrices_minerals[1][1][2] = 0.6;
    rotation_matrices_minerals[1][2][0] = 0.7;
    rotation_matrices_minerals[1][2][1] = 0.8;
    rotation_matrices_minerals[1][2][2] = 0.9;

    rotation_matrices_minerals[2][0][0] = 0.1;
    rotation_matrices_minerals[2][0][1] = 0.2;
    rotation_matrices_minerals[2][0][2] = 0.3;
    rotation_matrices_minerals[2][1][0] = 0.4;
    rotation_matrices_minerals[2][1][1] = 0.5;
    rotation_matrices_minerals[2][1][2] = 0.6;
    rotation_matrices_minerals[2][2][0] = 0.7;
    rotation_matrices_minerals[2][2][1] = 0.8;
    rotation_matrices_minerals[2][2][2] = 0.9;

    rotation_matrices_minerals[3][0][0] = 0.1;
    rotation_matrices_minerals[3][0][1] = 0.2;
    rotation_matrices_minerals[3][0][2] = 0.3;
    rotation_matrices_minerals[3][1][0] = 0.4;
    rotation_matrices_minerals[3][1][1] = 0.5;
    rotation_matrices_minerals[3][1][2] = 0.6;
    rotation_matrices_minerals[3][2][0] = 0.7;
    rotation_matrices_minerals[3][2][1] = 0.8;
    rotation_matrices_minerals[3][2][2] = 0.9;

    rotation_matrices_minerals[4][0][0] = 0.1;
    rotation_matrices_minerals[4][0][1] = 0.2;
    rotation_matrices_minerals[4][0][2] = 0.3;
    rotation_matrices_minerals[4][1][0] = 0.4;
    rotation_matrices_minerals[4][1][1] = 0.5;
    rotation_matrices_minerals[4][1][2] = 0.6;
    rotation_matrices_minerals[4][2][0] = 0.7;
    rotation_matrices_minerals[4][2][1] = 0.8;
    rotation_matrices_minerals[4][2][2] = 0.9;

    SymmetricTensor<2,3> strain_rate_nondimensional;
    strain_rate_nondimensional[0][1] = 0.5959;

    SymmetricTensor<2,3> compressible_strain_rate;

    Tensor<2,3> velocity_gradient_tensor_nondimensional;
    velocity_gradient_tensor_nondimensional[0][1] = 2.0* 0.5959;
    velocity_gradient_tensor_nondimensional[1][0] = 2.0* 0.5959;

    const Point<3> position;
    const double temperature = 0;
    const double pressure = 0;
    const Tensor<1,3> velocity;
    const std::vector<double> compositions;
    const SymmetricTensor<2,3> strain_rate;
    const double water_content = 0;

    std::pair<std::vector<double>, std::vector<Tensor<2,3>>> derivatives;
    derivatives = cpo_3d.compute_derivatives(0,
                                             data,
                                             0,
                                             strain_rate_nondimensional,
                                             velocity_gradient_tensor_nondimensional,
                                             position,
                                             temperature,
                                             pressure,
                                             velocity,
                                             compositions,
                                             strain_rate,
                                             compressible_strain_rate,
                                             water_content);

    // The correct analytical solution to check against
    // Note that this still has to be multiplied with with volume fraction
    // of each grain to get the same solution as D-Rex would get.
    double solution[5] = {0.0, 0.0, 0.0, 0.0 ,0.0};
    for (unsigned int iii = 0; iii < derivatives.first.size(); ++iii)
      CHECK(derivatives.first[iii] == Approx(solution[iii]));

    Tensor<2,3> deriv_direction_solution_2_to_5;
    deriv_direction_solution_2_to_5[0][0] = 0;
    deriv_direction_solution_2_to_5[0][1] = 0.0;
    deriv_direction_solution_2_to_5[0][2] = 0.0;
    deriv_direction_solution_2_to_5[1][0] = 0.0;
    deriv_direction_solution_2_to_5[1][1] = 0.0;
    deriv_direction_solution_2_to_5[1][2] = 0.0;
    deriv_direction_solution_2_to_5[2][0] = 0.0;
    deriv_direction_solution_2_to_5[2][1] = 0.0;
    deriv_direction_solution_2_to_5[2][2] = 0.0;
    Tensor<2,3> deriv_direction_solution[5] = {Tensor<2,3>(),deriv_direction_solution_2_to_5,deriv_direction_solution_2_to_5,deriv_direction_solution_2_to_5,deriv_direction_solution_2_to_5};
    for (size_t index = 0; index < 5; index++)
      {
        for (size_t iii = 0; iii < 3; iii++)
          {
            for (size_t jjj = 0; jjj < 3; jjj++)
              {
                CHECK(deriv_direction_solution[index][iii][jjj] == Approx(derivatives.second[index][iii][jjj]));
              }
          }
      }
    // now check if the value is the same as the D-Rex output when contracted with the direction matrix
    Tensor<2,3> deriv_direction_full_solution_2_to_5;
    deriv_direction_full_solution_2_to_5[0][0] = 0.0 ;
    deriv_direction_full_solution_2_to_5[0][1] = 0.0;
    deriv_direction_full_solution_2_to_5[0][2] = 0.0;
    deriv_direction_full_solution_2_to_5[1][0] = 0.0;
    deriv_direction_full_solution_2_to_5[1][1] = 0.0;
    deriv_direction_full_solution_2_to_5[1][2] = 0.0;
    deriv_direction_full_solution_2_to_5[2][0] = 0.0;
    deriv_direction_full_solution_2_to_5[2][1] = 0.0;
    deriv_direction_full_solution_2_to_5[2][2] = 0.0;
    Tensor<2,3> deriv_direction_full_solution[5] = {Tensor<2,3>(),deriv_direction_full_solution_2_to_5,deriv_direction_full_solution_2_to_5,deriv_direction_full_solution_2_to_5,deriv_direction_full_solution_2_to_5};
    for (size_t index = 0; index < 5; index++)
      {
        auto full_solution = rotation_matrices_minerals[index] * deriv_direction_solution[index];
        for (size_t iii = 0; iii < 3; iii++)
          {
            for (size_t jjj = 0; jjj < 3; jjj++)
              {
                CHECK(full_solution[iii][jjj] == Approx(deriv_direction_full_solution[index][iii][jjj]));
              }
          }
      }
  }
}


TEST_CASE("Fabric determination function")
{
  using namespace aspect;
  using namespace Particle::Property;
  CrystalPreferredOrientation<3> cpo;
  double MPa = 1e6;

  CHECK(cpo.determine_deformation_type_karato_2008(379.*MPa, 0.) == DeformationType::olivine_a_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(381.*MPa, 0.) == DeformationType::olivine_d_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(0.*MPa, 100.) == DeformationType::olivine_a_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(100.*MPa, 50.) == DeformationType::olivine_a_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(360.*MPa, 50.) == DeformationType::olivine_a_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(379.*MPa, 50.) == DeformationType::olivine_d_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(480.*MPa, 49.) == DeformationType::olivine_d_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(480.*MPa, 75.) == DeformationType::olivine_b_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(100.*MPa, 100.) == DeformationType::olivine_a_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(360.*MPa, 100.) == DeformationType::olivine_a_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(379.*MPa, 100.) == DeformationType::olivine_b_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(480.*MPa, 100.) == DeformationType::olivine_b_fabric);

  CHECK(cpo.determine_deformation_type_karato_2008(20.*MPa, 200.) == DeformationType::olivine_a_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(100.*MPa, 200.) == DeformationType::olivine_a_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(200.*MPa, 200.) == DeformationType::olivine_e_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(360.*MPa, 200.) == DeformationType::olivine_e_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(380.*MPa, 200.) == DeformationType::olivine_b_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(480.*MPa, 200.) == DeformationType::olivine_b_fabric);

  CHECK(cpo.determine_deformation_type_karato_2008(20.*MPa, 300.) == DeformationType::olivine_e_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(100.*MPa, 300.) == DeformationType::olivine_e_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(200.*MPa, 300.) == DeformationType::olivine_e_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(360.*MPa, 300.) == DeformationType::olivine_e_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(380.*MPa, 300.) == DeformationType::olivine_b_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(480.*MPa, 300.) == DeformationType::olivine_b_fabric);

  CHECK(cpo.determine_deformation_type_karato_2008(20.*MPa, 380.) == DeformationType::olivine_e_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(100.*MPa, 380.) == DeformationType::olivine_e_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(200.*MPa, 380.) == DeformationType::olivine_e_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(340.*MPa, 380.) == DeformationType::olivine_e_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(360.*MPa, 380.) == DeformationType::olivine_e_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(380.*MPa, 380.) == DeformationType::olivine_b_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(480.*MPa, 380.) == DeformationType::olivine_b_fabric);

  CHECK(cpo.determine_deformation_type_karato_2008(20.*MPa, 400.) == DeformationType::olivine_e_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(100.*MPa, 400.) == DeformationType::olivine_e_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(200.*MPa, 400.) == DeformationType::olivine_e_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(340.*MPa, 400.) == DeformationType::olivine_e_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(360.*MPa, 400.) == DeformationType::olivine_c_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(380.*MPa, 400.) == DeformationType::olivine_b_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(480.*MPa, 400.) == DeformationType::olivine_b_fabric);

  CHECK(cpo.determine_deformation_type_karato_2008(20.*MPa, 600.) == DeformationType::olivine_e_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(100.*MPa, 600.) == DeformationType::olivine_e_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(200.*MPa, 600.) == DeformationType::olivine_e_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(340.*MPa, 600.) == DeformationType::olivine_c_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(360.*MPa, 600.) == DeformationType::olivine_b_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(380.*MPa, 600.) == DeformationType::olivine_b_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(480.*MPa, 600.) == DeformationType::olivine_b_fabric);

  CHECK(cpo.determine_deformation_type_karato_2008(20.*MPa, 800.) == DeformationType::olivine_e_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(100.*MPa, 800.) == DeformationType::olivine_e_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(200.*MPa, 800.) == DeformationType::olivine_c_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(340.*MPa, 800.) == DeformationType::olivine_c_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(360.*MPa, 800.) == DeformationType::olivine_b_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(380.*MPa, 800.) == DeformationType::olivine_b_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(480.*MPa, 800.) == DeformationType::olivine_b_fabric);

  CHECK(cpo.determine_deformation_type_karato_2008(20.*MPa, 1000.) == DeformationType::olivine_e_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(100.*MPa, 1000.) == DeformationType::olivine_c_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(200.*MPa, 1000.) == DeformationType::olivine_c_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(340.*MPa, 1000.) == DeformationType::olivine_b_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(360.*MPa, 1000.) == DeformationType::olivine_b_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(380.*MPa, 1000.) == DeformationType::olivine_b_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(480.*MPa, 1000.) == DeformationType::olivine_b_fabric);

  CHECK(cpo.determine_deformation_type_karato_2008(20.*MPa, 1200.) == DeformationType::olivine_c_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(100.*MPa, 1200.) == DeformationType::olivine_c_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(200.*MPa, 1200.) == DeformationType::olivine_c_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(340.*MPa, 1200.) == DeformationType::olivine_b_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(360.*MPa, 1200.) == DeformationType::olivine_b_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(380.*MPa, 1200.) == DeformationType::olivine_b_fabric);
  CHECK(cpo.determine_deformation_type_karato_2008(480.*MPa, 1200.) == DeformationType::olivine_b_fabric);
}


TEST_CASE("CPO")
{
  using namespace dealii;
  using namespace aspect;

  {
    // test initialization 3d.
    const int dim3=3;

    Particle::Property::CrystalPreferredOrientation<dim3> lpo_3d;
    ParameterHandler prm;
    lpo_3d.declare_parameters(prm);

    prm.enter_subsection("Postprocess");
    {
      prm.enter_subsection("Particles");
      {
        prm.enter_subsection("Crystal Preferred Orientation");
        {
          prm.set("Random number seed","1");
          prm.set("Number of grains per particle","5");


        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }
    prm.leave_subsection ();

    lpo_3d.parse_parameters(prm);
    lpo_3d.initialize();


    Point<dim3> dummy_point;
    std::vector<double> data;
    lpo_3d.initialize_one_particle_property(dummy_point, data);

    // The LPO particles are initialized. With the same seed, the outcome should
    // always be the same, so test that for seed = 1. Forthermore, in the data
    // I can only really test that the first entry is the water content (0) and
    // that every first entry of each particle is 1/n_grains = 1/10 = 0.1.
    CHECK(data[0] == Approx(-1.0));//isnan(data[0])); // default fabric type which is only computed on a update
    CHECK(data[1] == Approx(0.5)); // default volume fraction olivine
    CHECK(data[2] == Approx(0.2));
    CHECK(data[3] == Approx(0.8794381492));
    CHECK(data[4] == Approx(-0.4757482965));
    CHECK(data[5] == Approx(-0.015877661));
    CHECK(data[6] == Approx(-0.197154988));
    CHECK(data[7] == Approx(-0.394402755));
    CHECK(data[8] == Approx(0.8975390674));
    CHECK(data[9] == Approx(-0.4332648756));
    CHECK(data[10] == Approx(-0.7861997362));
    CHECK(data[11] == Approx(-0.4406489786));
    CHECK(data[12] == Approx(0.2));
    CHECK(data[13] == Approx(-0.911475687));
    CHECK(data[14] == Approx(0.4113536692));
    CHECK(data[15] == Approx(0.0004804047));
    CHECK(data[16] == Approx(-0.3056521145));
    CHECK(data[17] == Approx(-0.6780433033));
    CHECK(data[18] == Approx(0.6684564786));
    CHECK(data[19] == Approx(0.2752977604));
    CHECK(data[20] == Approx(0.6091349914));
    CHECK(data[21] == Approx(0.7437511045));
    CHECK(data[22] == Approx(0.2));
    CHECK(data[32] == Approx(0.2));
    CHECK(data[42] == Approx(0.2));
    CHECK(data[52] == Approx(-1.0));// default fabric type which is only computed on a update
    CHECK(data[53] == Approx(0.5)); // default volume fraction olivine
    CHECK(data[54] == Approx(0.2));
    CHECK(data[64] == Approx(0.2));
    CHECK(data[74] == Approx(0.2));
    CHECK(data[84] == Approx(0.2));
    CHECK(data[94] == Approx(0.2));
    CHECK(data[103] == Approx(-0.0688278144));

    std::vector<double> volume_fractions(5,0.2);
    std::vector<dealii::Tensor<2,3>> a_cosine_matrices(5);
    a_cosine_matrices[0][0][0] = 0.5;
    a_cosine_matrices[0][0][1] = 0.5;
    a_cosine_matrices[0][0][2] = 0.5;
    a_cosine_matrices[0][1][0] = 0.5;
    a_cosine_matrices[0][1][1] = 0.5;
    a_cosine_matrices[0][1][2] = 0.5;
    a_cosine_matrices[0][2][0] = 0.5;
    a_cosine_matrices[0][2][1] = 0.5;
    a_cosine_matrices[0][2][2] = 0.5;

    a_cosine_matrices[1][0][0] = 0.1;
    a_cosine_matrices[1][0][1] = 0.2;
    a_cosine_matrices[1][0][2] = 0.3;
    a_cosine_matrices[1][1][0] = 0.4;
    a_cosine_matrices[1][1][1] = 0.5;
    a_cosine_matrices[1][1][2] = 0.6;
    a_cosine_matrices[1][2][0] = 0.7;
    a_cosine_matrices[1][2][1] = 0.8;
    a_cosine_matrices[1][2][2] = 0.9;

    a_cosine_matrices[2][0][0] = 0.1;
    a_cosine_matrices[2][0][1] = 0.2;
    a_cosine_matrices[2][0][2] = 0.3;
    a_cosine_matrices[2][1][0] = 0.4;
    a_cosine_matrices[2][1][1] = 0.5;
    a_cosine_matrices[2][1][2] = 0.6;
    a_cosine_matrices[2][2][0] = 0.7;
    a_cosine_matrices[2][2][1] = 0.8;
    a_cosine_matrices[2][2][2] = 0.9;

    a_cosine_matrices[3][0][0] = 0.1;
    a_cosine_matrices[3][0][1] = 0.2;
    a_cosine_matrices[3][0][2] = 0.3;
    a_cosine_matrices[3][1][0] = 0.4;
    a_cosine_matrices[3][1][1] = 0.5;
    a_cosine_matrices[3][1][2] = 0.6;
    a_cosine_matrices[3][2][0] = 0.7;
    a_cosine_matrices[3][2][1] = 0.8;
    a_cosine_matrices[3][2][2] = 0.9;

    a_cosine_matrices[4][0][0] = 0.1;
    a_cosine_matrices[4][0][1] = 0.2;
    a_cosine_matrices[4][0][2] = 0.3;
    a_cosine_matrices[4][1][0] = 0.4;
    a_cosine_matrices[4][1][1] = 0.5;
    a_cosine_matrices[4][1][2] = 0.6;
    a_cosine_matrices[4][2][0] = 0.7;
    a_cosine_matrices[4][2][1] = 0.8;
    a_cosine_matrices[4][2][2] = 0.9;

    lpo_3d.set_volume_fractions_grains(0,data,0,0,0.2);
    lpo_3d.set_volume_fractions_grains(0,data,0,1,0.2);
    lpo_3d.set_volume_fractions_grains(0,data,0,2,0.2);
    lpo_3d.set_volume_fractions_grains(0,data,0,3,0.2);
    lpo_3d.set_volume_fractions_grains(0,data,0,4,0.2);

    lpo_3d.set_rotation_matrix_grains(0,data,0,0,a_cosine_matrices[0]);
    lpo_3d.set_rotation_matrix_grains(0,data,0,1,a_cosine_matrices[1]);
    lpo_3d.set_rotation_matrix_grains(0,data,0,2,a_cosine_matrices[2]);
    lpo_3d.set_rotation_matrix_grains(0,data,0,3,a_cosine_matrices[3]);
    lpo_3d.set_rotation_matrix_grains(0,data,0,4,a_cosine_matrices[4]);

    SymmetricTensor<2,3> strain_rate_nondimensional;
    strain_rate_nondimensional[0][1] = 0.5959;

    Tensor<2,3> velocity_gradient_tensor_nondimensional;
    velocity_gradient_tensor_nondimensional[0][1] = 2.0* 0.5959;
    velocity_gradient_tensor_nondimensional[1][0] = 2.0* 0.5959;

    std::array<double,4> ref_resolved_shear_stress;
    ref_resolved_shear_stress[0] = 1;
    ref_resolved_shear_stress[1] = 2;
    ref_resolved_shear_stress[2] = 3;
    ref_resolved_shear_stress[3] = 1e60; // can't really use nummerical limits max or infinite, because need to be able to square it without becoming infinite. This is the value fortran D-Rex uses.


    std::pair<std::vector<double>, std::vector<Tensor<2,3>>> derivatives;
    derivatives = lpo_3d.compute_derivatives_drex_2004(0,
                                                       data,
                                                       0,
                                                       strain_rate_nondimensional,
                                                       velocity_gradient_tensor_nondimensional,
                                                       ref_resolved_shear_stress,
                                                       true);


    // The correct analytical solution to check against
    double solution[5] = {3.150563756, -0.787640939, -0.787640939, -0.787640939 ,-0.787640939};
    for (unsigned int i = 0; i < derivatives.first.size(); ++i)
      CHECK(derivatives.first[i] == Approx(solution[i]));


    Tensor<2,3> deriv_direction_solution_2_to_5;
    deriv_direction_solution_2_to_5[0][0] = 0;
    deriv_direction_solution_2_to_5[0][1] = 0.00501912;
    deriv_direction_solution_2_to_5[0][2] = 0.0100382;
    deriv_direction_solution_2_to_5[1][0] = -0.00501912;
    deriv_direction_solution_2_to_5[1][1] = 0;
    deriv_direction_solution_2_to_5[1][2] = 0.00501912;
    deriv_direction_solution_2_to_5[2][0] = -0.0100382;
    deriv_direction_solution_2_to_5[2][1] = -0.00501912;
    deriv_direction_solution_2_to_5[2][2] = 0;

    Tensor<2,3> deriv_direction_solution[5] = {Tensor<2,3>(),deriv_direction_solution_2_to_5,deriv_direction_solution_2_to_5,deriv_direction_solution_2_to_5,deriv_direction_solution_2_to_5};
    for (size_t index = 0; index < 5; index++)
      {
        for (size_t iii = 0; iii < 3; iii++)
          {
            for (size_t jjj = 0; jjj < 3; jjj++)
              {
                CHECK(deriv_direction_solution[index][iii][jjj] == Approx(derivatives.second[index][iii][jjj]));
              }
          }
      }
    // now check if the value is the same as the D-Rex output when contracted with the direction matrix
    Tensor<2,3> deriv_direction_full_solution_2_to_5;
    deriv_direction_full_solution_2_to_5[0][0] = -0.004015284 ;
    deriv_direction_full_solution_2_to_5[0][1] = -0.0010038242;
    deriv_direction_full_solution_2_to_5[0][2] = 0.0020076485;
    deriv_direction_full_solution_2_to_5[1][0] = -0.00853248;
    deriv_direction_full_solution_2_to_5[1][1] = -0.0010038242;
    deriv_direction_full_solution_2_to_5[1][2] = 0.0065248576;
    deriv_direction_full_solution_2_to_5[2][0] = -0.0130497151;
    deriv_direction_full_solution_2_to_5[2][1] = -0.0010038242;
    deriv_direction_full_solution_2_to_5[2][2] = 0.0110420667;
    Tensor<2,3> deriv_direction_full_solution[5] = {Tensor<2,3>(),deriv_direction_full_solution_2_to_5,deriv_direction_full_solution_2_to_5,deriv_direction_full_solution_2_to_5,deriv_direction_full_solution_2_to_5};
    for (size_t index = 0; index < 5; index++)
      {
        auto full_solution = a_cosine_matrices[index] * deriv_direction_solution[index];
        for (size_t iii = 0; iii < 3; iii++)
          {
            for (size_t jjj = 0; jjj < 3; jjj++)
              {
                CHECK(full_solution[iii][jjj] == Approx(deriv_direction_full_solution[index][iii][jjj]));
              }
          }
      }
  }

  {
    // thirdly test 3d lpo with different strain-rate and velocity gradient tensors.
    // The solution of this test has not been analytically confirmed.
    const int dim3=3;

    Particle::Property::CrystalPreferredOrientation<dim3> lpo_3d;
    ParameterHandler prm;
    lpo_3d.declare_parameters(prm);

    prm.enter_subsection("Postprocess");
    {
      prm.enter_subsection("Particles");
      {
        prm.enter_subsection("Crystal Preferred Orientation");
        {
          prm.set("Random number seed","1");
          prm.set("Number of grains per particle","5");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }
    prm.leave_subsection ();

    lpo_3d.parse_parameters(prm);
    lpo_3d.initialize();


    Point<dim3> dummy_point;
    std::vector<double> data;
    lpo_3d.initialize_one_particle_property(dummy_point, data);

    // The LPO particles are initialized. With the same seed, the outcome should
    // always be the same, so test that for seed = 1. Forthermore, in the data
    // I can only really test that the first entry is the water content (0) and
    // that every first entry of each particle is 1/n_grains = 1/10 = 0.1.
    CHECK(data[0] == Approx(-1.0));//isnan(data[0])); // default fabric type which is only computed on a update
    CHECK(data[1] == Approx(0.5)); // default volume fraction olivine
    CHECK(data[2] == Approx(0.2));
    CHECK(data[3] == Approx(0.8794381492));
    CHECK(data[4] == Approx(-0.4757482965));
    CHECK(data[5] == Approx(-0.015877661));
    CHECK(data[6] == Approx(-0.197154988));
    CHECK(data[7] == Approx(-0.394402755));
    CHECK(data[8] == Approx(0.8975390674));
    CHECK(data[9] == Approx(-0.4332648756));
    CHECK(data[10] == Approx(-0.7861997362));
    CHECK(data[11] == Approx(-0.4406489786));
    CHECK(data[12] == Approx(0.2));
    CHECK(data[13] == Approx(-0.911475687));
    CHECK(data[14] == Approx(0.4113536692));
    CHECK(data[15] == Approx(0.0004804047));
    CHECK(data[16] == Approx(-0.3056521145));
    CHECK(data[17] == Approx(-0.6780433033));
    CHECK(data[18] == Approx(0.6684564786));
    CHECK(data[19] == Approx(0.2752977604));
    CHECK(data[20] == Approx(0.6091349914));
    CHECK(data[21] == Approx(0.7437511045));
    CHECK(data[22] == Approx(0.2));
    CHECK(data[32] == Approx(0.2));
    CHECK(data[42] == Approx(0.2));
    CHECK(data[52] == Approx(-1.0));// default fabric type which is only computed on a update
    CHECK(data[53] == Approx(0.5)); // default volume fraction olivine
    CHECK(data[54] == Approx(0.2));
    CHECK(data[64] == Approx(0.2));
    CHECK(data[74] == Approx(0.2));
    CHECK(data[84] == Approx(0.2));
    CHECK(data[94] == Approx(0.2));
    CHECK(data[103] == Approx(-0.0688278144));

    std::vector<double> volume_fractions(5,0.2);
    std::vector<dealii::Tensor<2,3>> a_cosine_matrices(5);
    a_cosine_matrices[0][0][0] = 0.5;
    a_cosine_matrices[0][0][1] = 0.5;
    a_cosine_matrices[0][0][2] = 0.5;
    a_cosine_matrices[0][1][0] = 0.5;
    a_cosine_matrices[0][1][1] = 0.5;
    a_cosine_matrices[0][1][2] = 0.5;
    a_cosine_matrices[0][2][0] = 0.5;
    a_cosine_matrices[0][2][1] = 0.5;
    a_cosine_matrices[0][2][2] = 0.5;

    a_cosine_matrices[1][0][0] = 0.1;
    a_cosine_matrices[1][0][1] = 0.2;
    a_cosine_matrices[1][0][2] = 0.3;
    a_cosine_matrices[1][1][0] = 0.4;
    a_cosine_matrices[1][1][1] = 0.5;
    a_cosine_matrices[1][1][2] = 0.6;
    a_cosine_matrices[1][2][0] = 0.7;
    a_cosine_matrices[1][2][1] = 0.8;
    a_cosine_matrices[1][2][2] = 0.9;

    a_cosine_matrices[2][0][0] = 0.1;
    a_cosine_matrices[2][0][1] = 0.2;
    a_cosine_matrices[2][0][2] = 0.3;
    a_cosine_matrices[2][1][0] = 0.4;
    a_cosine_matrices[2][1][1] = 0.5;
    a_cosine_matrices[2][1][2] = 0.6;
    a_cosine_matrices[2][2][0] = 0.7;
    a_cosine_matrices[2][2][1] = 0.8;
    a_cosine_matrices[2][2][2] = 0.9;

    a_cosine_matrices[3][0][0] = 0.1;
    a_cosine_matrices[3][0][1] = 0.2;
    a_cosine_matrices[3][0][2] = 0.3;
    a_cosine_matrices[3][1][0] = 0.4;
    a_cosine_matrices[3][1][1] = 0.5;
    a_cosine_matrices[3][1][2] = 0.6;
    a_cosine_matrices[3][2][0] = 0.7;
    a_cosine_matrices[3][2][1] = 0.8;
    a_cosine_matrices[3][2][2] = 0.9;

    a_cosine_matrices[4][0][0] = 0.1;
    a_cosine_matrices[4][0][1] = 0.2;
    a_cosine_matrices[4][0][2] = 0.3;
    a_cosine_matrices[4][1][0] = 0.4;
    a_cosine_matrices[4][1][1] = 0.5;
    a_cosine_matrices[4][1][2] = 0.6;
    a_cosine_matrices[4][2][0] = 0.7;
    a_cosine_matrices[4][2][1] = 0.8;
    a_cosine_matrices[4][2][2] = 0.9;


    lpo_3d.set_volume_fractions_grains(0,data,0,0,0.2);
    lpo_3d.set_volume_fractions_grains(0,data,0,1,0.2);
    lpo_3d.set_volume_fractions_grains(0,data,0,2,0.2);
    lpo_3d.set_volume_fractions_grains(0,data,0,3,0.2);
    lpo_3d.set_volume_fractions_grains(0,data,0,4,0.2);

    lpo_3d.set_rotation_matrix_grains(0,data,0,0,a_cosine_matrices[0]);
    lpo_3d.set_rotation_matrix_grains(0,data,0,1,a_cosine_matrices[1]);
    lpo_3d.set_rotation_matrix_grains(0,data,0,2,a_cosine_matrices[2]);
    lpo_3d.set_rotation_matrix_grains(0,data,0,3,a_cosine_matrices[3]);
    lpo_3d.set_rotation_matrix_grains(0,data,0,4,a_cosine_matrices[4]);

    SymmetricTensor<2,dim3> strain_rate_nondimensional; // e
    strain_rate_nondimensional[0][0] = 7.5;
    strain_rate_nondimensional[0][1] = 8;
    strain_rate_nondimensional[0][2] = 8.5;
    //strain_rate_nondimensional[1][0] = 8; // symmetry
    strain_rate_nondimensional[1][1] = 9.5;
    strain_rate_nondimensional[1][2] = 10;
    //strain_rate_nondimensional[2][0] = 8.5;// symmetry
    //strain_rate_nondimensional[2][1] = 10; // symmetry
    strain_rate_nondimensional[2][2] = 11.5;

    Tensor<2,dim3> velocity_gradient_tensor_nondimensional; // l
    velocity_gradient_tensor_nondimensional[0][0] = 2;
    velocity_gradient_tensor_nondimensional[0][1] = 2.5;
    velocity_gradient_tensor_nondimensional[0][2] = 3;
    velocity_gradient_tensor_nondimensional[1][0] = 3.5;
    velocity_gradient_tensor_nondimensional[1][1] = 4;
    velocity_gradient_tensor_nondimensional[1][2] = 4.5;
    velocity_gradient_tensor_nondimensional[2][0] = 5;
    velocity_gradient_tensor_nondimensional[2][1] = 5.5;
    velocity_gradient_tensor_nondimensional[2][2] = 6;


    std::array<double,4> ref_resolved_shear_stress;
    ref_resolved_shear_stress[0] = 1;
    ref_resolved_shear_stress[1] = 2;
    ref_resolved_shear_stress[2] = 3;
    ref_resolved_shear_stress[3] = 1e60; // can't really use nummerical limits max or infinite, because need to be able to square it without becoming infinite. This is the value fortran D-Rex uses.


    std::pair<std::vector<double>, std::vector<Tensor<2,3>>> derivatives;
    derivatives = lpo_3d.compute_derivatives_drex_2004(0,
                                                       data,
                                                       0,
                                                       strain_rate_nondimensional,
                                                       velocity_gradient_tensor_nondimensional,
                                                       ref_resolved_shear_stress,
                                                       true);

    // The correct analytical solution to check against
    double solution[5] = {2.5350823696, -0.6337705924, -0.6337705924, -0.6337705924 ,-0.6337705924};
    for (unsigned int i = 0; i < derivatives.first.size(); ++i)
      CHECK(derivatives.first[i] == Approx(solution[i]));


    Tensor<2,3> deriv_direction_solution_1;
    deriv_direction_solution_1[0][0] = 0;
    deriv_direction_solution_1[0][1] = 0.5;
    deriv_direction_solution_1[0][2] = 1.0;
    deriv_direction_solution_1[1][0] = -0.5;
    deriv_direction_solution_1[1][1] = 0;
    deriv_direction_solution_1[1][2] = 0.5;
    deriv_direction_solution_1[2][0] = -1.0;
    deriv_direction_solution_1[2][1] = -0.5;
    deriv_direction_solution_1[2][2] = 0;
    Tensor<2,3> deriv_direction_solution_2_to_5;
    deriv_direction_solution_2_to_5[0][0] = 0;
    deriv_direction_solution_2_to_5[0][1] = 0.5304624591;
    deriv_direction_solution_2_to_5[0][2] = 1.0609249182;
    deriv_direction_solution_2_to_5[1][0] = -0.5304624591;
    deriv_direction_solution_2_to_5[1][1] = 0;
    deriv_direction_solution_2_to_5[1][2] = 0.5304624591;
    deriv_direction_solution_2_to_5[2][0] = -1.0609249182;
    deriv_direction_solution_2_to_5[2][1] = -0.5304624591 ;
    deriv_direction_solution_2_to_5[2][2] = 0;

    Tensor<2,3> deriv_direction_solution[5] = {deriv_direction_solution_1,deriv_direction_solution_2_to_5,deriv_direction_solution_2_to_5,deriv_direction_solution_2_to_5,deriv_direction_solution_2_to_5};
    for (size_t index = 0; index < 5; index++)
      {
        for (size_t iii = 0; iii < 3; iii++)
          {
            for (size_t jjj = 0; jjj < 3; jjj++)
              {
                CHECK(deriv_direction_solution[index][iii][jjj] == Approx(derivatives.second[index][iii][jjj]));
              }
          }
      }
    // now check if the value is the same as the D-Rex output when contracted with the direction matrix
    Tensor<2,3> deriv_direction_full_solution_1;
    deriv_direction_full_solution_1[0][0] = -0.75;
    deriv_direction_full_solution_1[0][1] = 0.0;
    deriv_direction_full_solution_1[0][2] = 0.75;
    deriv_direction_full_solution_1[1][0] = -0.75;
    deriv_direction_full_solution_1[1][1] = 0;
    deriv_direction_full_solution_1[1][2] = 0.75;
    deriv_direction_full_solution_1[2][0] = -0.75;
    deriv_direction_full_solution_1[2][1] = 0.0;
    deriv_direction_full_solution_1[2][2] = 0.75;
    Tensor<2,3> deriv_direction_full_solution_2_to_5;
    deriv_direction_full_solution_2_to_5[0][0] = -0.4243699673 ;
    deriv_direction_full_solution_2_to_5[0][1] = -0.1060924918;
    deriv_direction_full_solution_2_to_5[0][2] = 0.2121849836;
    deriv_direction_full_solution_2_to_5[1][0] = -0.9017861805;
    deriv_direction_full_solution_2_to_5[1][1] = -0.1060924918;
    deriv_direction_full_solution_2_to_5[1][2] = 0.6896011968;
    deriv_direction_full_solution_2_to_5[2][0] = -1.3792023937;
    deriv_direction_full_solution_2_to_5[2][1] = -0.1060924918;
    deriv_direction_full_solution_2_to_5[2][2] = 1.16701741;
    Tensor<2,3> deriv_direction_full_solution[5] = {deriv_direction_full_solution_1,deriv_direction_full_solution_2_to_5,deriv_direction_full_solution_2_to_5,deriv_direction_full_solution_2_to_5,deriv_direction_full_solution_2_to_5};
    for (size_t index = 0; index < 5; index++)
      {
        auto full_solution = a_cosine_matrices[index] * deriv_direction_solution[index];
        for (size_t iii = 0; iii < 3; iii++)
          {
            for (size_t jjj = 0; jjj < 3; jjj++)
              {
                CHECK(full_solution[iii][jjj] == Approx(deriv_direction_full_solution[index][iii][jjj]));
              }
          }
      }
  }
}

TEST_CASE("CPO elastic tensor transform functions")
{
  dealii::SymmetricTensor<2,6> reference_elastic_tensor({1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21});

// first test whether the functions are invertable
  {
    dealii::SymmetricTensor<2,6> result_up_down = aspect::Particle::Property::CpoElasticTensor<3>::transform_4th_order_tensor_to_6x6_matrix(aspect::Particle::Property::CpoElasticTensor<3>::transform_6x6_matrix_to_4th_order_tensor(reference_elastic_tensor));

    for (size_t i = 0; i < 6; i++)
      {
        for (size_t j = 0; j < 6; j++)
          {
            REQUIRE(reference_elastic_tensor[i][j] == Approx(result_up_down[i][j]));
          }
      }
  }
  {
    dealii::SymmetricTensor<2,6> result_down_up = aspect::Particle::Property::CpoElasticTensor<3>::transform_21D_vector_to_6x6_matrix(aspect::Particle::Property::CpoElasticTensor<3>::transform_6x6_matrix_to_21D_vector(reference_elastic_tensor));

    for (size_t i = 0; i < 6; i++)
      {
        for (size_t j = 0; j < 6; j++)
          {
            REQUIRE(reference_elastic_tensor[i][j] == Approx(result_down_up[i][j]));
          }
      }
  }
  {
    dealii::SymmetricTensor<2,6> result_up_2down_up = aspect::Particle::Property::CpoElasticTensor<3>::transform_21D_vector_to_6x6_matrix(aspect::Particle::Property::CpoElasticTensor<3>::transform_4th_order_tensor_to_21D_vector(aspect::Particle::Property::CpoElasticTensor<3>::transform_6x6_matrix_to_4th_order_tensor(reference_elastic_tensor)));

    for (size_t i = 0; i < 6; i++)
      {
        for (size_t j = 0; j < 6; j++)
          {
            REQUIRE(reference_elastic_tensor[i][j] == Approx(result_up_2down_up[i][j]));
          }
      }
  }

  // test rotations
  // rotation matrix
  dealii::Tensor<2,3> rotation_tensor;

  {
    // fill the rotation matrix with a rotations in all directions
    {
      double radians = (dealii::numbers::PI/180.0)*(360/5); //0.35*dealii::numbers::PI; //(dealii::numbers::PI/180.0)*36;
      double alpha = radians;
      double beta = radians;
      double gamma = radians;
      rotation_tensor[0][0] = cos(alpha) * cos(beta);
      rotation_tensor[0][1] = sin(alpha) * cos(beta);
      rotation_tensor[0][2] = -sin(beta);
      rotation_tensor[1][0] = cos(alpha) * sin(beta) * sin(gamma) - sin(alpha)*cos(gamma);
      rotation_tensor[1][1] = sin(alpha) * sin(beta) * sin(gamma) + cos(alpha)*cos(gamma);
      rotation_tensor[1][2] = cos(beta) * sin(gamma);
      rotation_tensor[2][0] = cos(alpha) * sin(beta) * cos(gamma) + sin(alpha)*sin(gamma);
      rotation_tensor[2][1] = sin(alpha) * sin(beta) * cos(gamma) - cos(alpha)*sin(gamma);
      rotation_tensor[2][2] = cos(beta) * cos(gamma);
    }

    {
      dealii::SymmetricTensor<2,6> result_up_1_rotate_down = aspect::Particle::Property::CpoElasticTensor<3>::transform_4th_order_tensor_to_6x6_matrix(
                                                               aspect::Particle::Property::CpoElasticTensor<3>::rotate_4th_order_tensor(
                                                                 aspect::Particle::Property::CpoElasticTensor<3>::transform_6x6_matrix_to_4th_order_tensor(reference_elastic_tensor),rotation_tensor));
      dealii::SymmetricTensor<2,6> result_1_rotate = aspect::Particle::Property::CpoElasticTensor<3>::rotate_6x6_matrix(reference_elastic_tensor,rotation_tensor);

      for (size_t i = 0; i < 6; i++)
        {
          for (size_t j = 0; j < 6; j++)
            {
              REQUIRE(result_1_rotate[i][j] == Approx(result_up_1_rotate_down[i][j]));
            }
        }

      dealii::Tensor<4,3> result_up_10_rotate = aspect::Particle::Property::CpoElasticTensor<3>::transform_6x6_matrix_to_4th_order_tensor(result_up_1_rotate_down);

      dealii::SymmetricTensor<2,6> result_5_rotate = result_1_rotate;

      for (size_t i = 0; i < 4; i++)
        {
          result_up_10_rotate = aspect::Particle::Property::CpoElasticTensor<3>::rotate_4th_order_tensor(result_up_10_rotate, rotation_tensor);
          result_5_rotate = aspect::Particle::Property::CpoElasticTensor<3>::rotate_6x6_matrix(result_5_rotate, rotation_tensor);
        }

      dealii::SymmetricTensor<2,6> result_up_10_rotate_down = aspect::Particle::Property::CpoElasticTensor<3>::transform_4th_order_tensor_to_6x6_matrix(result_up_10_rotate);

      for (size_t i = 0; i < 6; i++)
        {
          for (size_t j = 0; j < 6; j++)
            {
              REQUIRE(result_5_rotate[i][j] == Approx(result_up_10_rotate_down[i][j]));
              // This test doesn't work when rotating in all rotations at the same time.
              //REQUIRE(result_1_rotate[i][j] == Approx(reference_elastic_tensor[i][j]));
            }
        }
    }
  }
  {
    // fill the rotation matrix with a rotations in the alpha direction
    {
      double radians = (dealii::numbers::PI/180.0)*(360/5); //0.35*dealii::numbers::PI; //(dealii::numbers::PI/180.0)*36;
      double alpha = radians;
      double beta = 0;
      double gamma = 0;
      rotation_tensor[0][0] = cos(alpha) * cos(beta);
      rotation_tensor[0][1] = sin(alpha) * cos(beta);
      rotation_tensor[0][2] = -sin(beta);
      rotation_tensor[1][0] = cos(alpha) * sin(beta) * sin(gamma) - sin(alpha)*cos(gamma);
      rotation_tensor[1][1] = sin(alpha) * sin(beta) * sin(gamma) + cos(alpha)*cos(gamma);
      rotation_tensor[1][2] = cos(beta) * sin(gamma);
      rotation_tensor[2][0] = cos(alpha) * sin(beta) * cos(gamma) + sin(alpha)*sin(gamma);
      rotation_tensor[2][1] = sin(alpha) * sin(beta) * cos(gamma) - cos(alpha)*sin(gamma);
      rotation_tensor[2][2] = cos(beta) * cos(gamma);
    }

    {
      dealii::SymmetricTensor<2,6> result_up_1_rotate_down = aspect::Particle::Property::CpoElasticTensor<3>::transform_4th_order_tensor_to_6x6_matrix(
                                                               aspect::Particle::Property::CpoElasticTensor<3>::rotate_4th_order_tensor(
                                                                 aspect::Particle::Property::CpoElasticTensor<3>::transform_6x6_matrix_to_4th_order_tensor(reference_elastic_tensor),rotation_tensor));
      dealii::SymmetricTensor<2,6> result_1_rotate = aspect::Particle::Property::CpoElasticTensor<3>::rotate_6x6_matrix(reference_elastic_tensor,rotation_tensor);

      for (size_t i = 0; i < 6; i++)
        {
          for (size_t j = 0; j < 6; j++)
            {
              REQUIRE(result_1_rotate[i][j] == Approx(result_up_1_rotate_down[i][j]));
            }
        }

      dealii::Tensor<4,3> result_up_10_rotate = aspect::Particle::Property::CpoElasticTensor<3>::transform_6x6_matrix_to_4th_order_tensor(result_up_1_rotate_down);

      dealii::SymmetricTensor<2,6> result_5_rotate = result_1_rotate;

      for (size_t i = 0; i < 4; i++)
        {
          result_up_10_rotate = aspect::Particle::Property::CpoElasticTensor<3>::rotate_4th_order_tensor(result_up_10_rotate, rotation_tensor);
          result_5_rotate = aspect::Particle::Property::CpoElasticTensor<3>::rotate_6x6_matrix(result_5_rotate, rotation_tensor);
        }

      dealii::SymmetricTensor<2,6> result_up_10_rotate_down = aspect::Particle::Property::CpoElasticTensor<3>::transform_4th_order_tensor_to_6x6_matrix(result_up_10_rotate);

      for (size_t i = 0; i < 6; i++)
        {
          for (size_t j = 0; j < 6; j++)
            {
              REQUIRE(result_5_rotate[i][j] == Approx(result_up_10_rotate_down[i][j]));
              REQUIRE(result_5_rotate[i][j] == Approx(reference_elastic_tensor[i][j]));
            }
        }
    }
  }
  {
    // fill the rotation matrix with a rotations in the beta direction
    {
      double radians = (dealii::numbers::PI/180.0)*(360/5); //0.35*dealii::numbers::PI; //(dealii::numbers::PI/180.0)*36;
      double alpha = 0;
      double beta = radians;
      double gamma = 0;
      rotation_tensor[0][0] = cos(alpha) * cos(beta);
      rotation_tensor[0][1] = sin(alpha) * cos(beta);
      rotation_tensor[0][2] = -sin(beta);
      rotation_tensor[1][0] = cos(alpha) * sin(beta) * sin(gamma) - sin(alpha)*cos(gamma);
      rotation_tensor[1][1] = sin(alpha) * sin(beta) * sin(gamma) + cos(alpha)*cos(gamma);
      rotation_tensor[1][2] = cos(beta) * sin(gamma);
      rotation_tensor[2][0] = cos(alpha) * sin(beta) * cos(gamma) + sin(alpha)*sin(gamma);
      rotation_tensor[2][1] = sin(alpha) * sin(beta) * cos(gamma) - cos(alpha)*sin(gamma);
      rotation_tensor[2][2] = cos(beta) * cos(gamma);
    }

    {
      dealii::SymmetricTensor<2,6> result_up_1_rotate_down = aspect::Particle::Property::CpoElasticTensor<3>::transform_4th_order_tensor_to_6x6_matrix(
                                                               aspect::Particle::Property::CpoElasticTensor<3>::rotate_4th_order_tensor(
                                                                 aspect::Particle::Property::CpoElasticTensor<3>::transform_6x6_matrix_to_4th_order_tensor(reference_elastic_tensor),rotation_tensor));
      dealii::SymmetricTensor<2,6> result_1_rotate = aspect::Particle::Property::CpoElasticTensor<3>::rotate_6x6_matrix(reference_elastic_tensor,rotation_tensor);

      for (size_t i = 0; i < 6; i++)
        {
          for (size_t j = 0; j < 6; j++)
            {
              REQUIRE(result_1_rotate[i][j] == Approx(result_up_1_rotate_down[i][j]));
            }
        }

      dealii::Tensor<4,3> result_up_10_rotate = aspect::Particle::Property::CpoElasticTensor<3>::transform_6x6_matrix_to_4th_order_tensor(result_up_1_rotate_down);

      dealii::SymmetricTensor<2,6> result_5_rotate = result_1_rotate;

      for (size_t i = 0; i < 4; i++)
        {
          result_up_10_rotate = aspect::Particle::Property::CpoElasticTensor<3>::rotate_4th_order_tensor(result_up_10_rotate, rotation_tensor);
          result_5_rotate = aspect::Particle::Property::CpoElasticTensor<3>::rotate_6x6_matrix(result_5_rotate, rotation_tensor);
        }

      dealii::SymmetricTensor<2,6> result_up_10_rotate_down = aspect::Particle::Property::CpoElasticTensor<3>::transform_4th_order_tensor_to_6x6_matrix(result_up_10_rotate);

      for (size_t i = 0; i < 6; i++)
        {
          for (size_t j = 0; j < 6; j++)
            {
              REQUIRE(result_5_rotate[i][j] == Approx(result_up_10_rotate_down[i][j]));
              REQUIRE(result_5_rotate[i][j] == Approx(reference_elastic_tensor[i][j]));
            }
        }
    }
  }

  {
    // fill the rotation matrix with a rotations in the gamma direction
    {
      double radians = (dealii::numbers::PI/180.0)*(360/5); //0.35*dealii::numbers::PI; //(dealii::numbers::PI/180.0)*36;
      double alpha = 0;
      double beta = 0;
      double gamma = radians;
      rotation_tensor[0][0] = cos(alpha) * cos(beta);
      rotation_tensor[0][1] = sin(alpha) * cos(beta);
      rotation_tensor[0][2] = -sin(beta);
      rotation_tensor[1][0] = cos(alpha) * sin(beta) * sin(gamma) - sin(alpha)*cos(gamma);
      rotation_tensor[1][1] = sin(alpha) * sin(beta) * sin(gamma) + cos(alpha)*cos(gamma);
      rotation_tensor[1][2] = cos(beta) * sin(gamma);
      rotation_tensor[2][0] = cos(alpha) * sin(beta) * cos(gamma) + sin(alpha)*sin(gamma);
      rotation_tensor[2][1] = sin(alpha) * sin(beta) * cos(gamma) - cos(alpha)*sin(gamma);
      rotation_tensor[2][2] = cos(beta) * cos(gamma);
    }

    {
      dealii::SymmetricTensor<2,6> result_up_1_rotate_down = aspect::Particle::Property::CpoElasticTensor<3>::transform_4th_order_tensor_to_6x6_matrix(
                                                               aspect::Particle::Property::CpoElasticTensor<3>::rotate_4th_order_tensor(
                                                                 aspect::Particle::Property::CpoElasticTensor<3>::transform_6x6_matrix_to_4th_order_tensor(reference_elastic_tensor),rotation_tensor));
      dealii::SymmetricTensor<2,6> result_1_rotate = aspect::Particle::Property::CpoElasticTensor<3>::rotate_6x6_matrix(reference_elastic_tensor,rotation_tensor);

      for (size_t i = 0; i < 6; i++)
        {
          for (size_t j = 0; j < 6; j++)
            {
              REQUIRE(result_1_rotate[i][j] == Approx(result_up_1_rotate_down[i][j]));
            }
        }

      dealii::Tensor<4,3> result_up_10_rotate = aspect::Particle::Property::CpoElasticTensor<3>::transform_6x6_matrix_to_4th_order_tensor(result_up_1_rotate_down);

      dealii::SymmetricTensor<2,6> result_5_rotate = result_1_rotate;

      for (size_t i = 0; i < 4; i++)
        {
          result_up_10_rotate = aspect::Particle::Property::CpoElasticTensor<3>::rotate_4th_order_tensor(result_up_10_rotate, rotation_tensor);
          result_5_rotate = aspect::Particle::Property::CpoElasticTensor<3>::rotate_6x6_matrix(result_5_rotate, rotation_tensor);
        }

      dealii::SymmetricTensor<2,6> result_up_10_rotate_down = aspect::Particle::Property::CpoElasticTensor<3>::transform_4th_order_tensor_to_6x6_matrix(result_up_10_rotate);

      for (size_t i = 0; i < 6; i++)
        {
          for (size_t j = 0; j < 6; j++)
            {
              REQUIRE(result_5_rotate[i][j] == Approx(result_up_10_rotate_down[i][j]));
              REQUIRE(result_5_rotate[i][j] == Approx(reference_elastic_tensor[i][j]));
            }
        }
    }
  }

}


TEST_CASE("CPO elastic tensor")
{

  std::vector<double> volume_fraction_mineral = {0.7,0.3};
  std::vector<std::vector<double>> volume_fractions_grains(2,std::vector<double>(8));
  std::vector<std::vector<dealii::Tensor<2,3>>> a_cosine_matrices_grains(2,std::vector<dealii::Tensor<2,3>>(8));

  dealii::Tensor<2,6> reference_elastic_tensor;
  dealii::Tensor<2,6> computed_elastic_tensor;


  unsigned int cpo_data_position = 0;
  std::vector<double> data_array(200,-1.);
  dealii::ArrayView<double> data_cpo(&data_array[0],200);


  // test the averaging function
  aspect::Particle::Property::CrystalPreferredOrientation<3> cpo;
  aspect::Particle::Property::CpoElasticTensor<3> cpo_elastic_tensor;
  aspect::ParameterHandler prm;
  cpo.declare_parameters(prm);
  cpo_elastic_tensor.declare_parameters(prm);
  prm.enter_subsection("Postprocess");
  {
    prm.enter_subsection("Particles");
    {
      prm.enter_subsection("Crystal Preferred Orientation");
      {
        prm.set("Random number seed","1");
        prm.set("Number of grains per particle","8");
        prm.set("CPO derivatives algorithm","Spin tensor");
        prm.set("Property advection method","Backward Euler");
        prm.enter_subsection("Initial grains");
        {
          prm.set("Model name","Uniform grains and random uniform rotations");
          prm.enter_subsection("Uniform grains and random uniform rotations");
          {
            // Let the minerals just passively rotate with the rotation of
            // the particle caused by the flow.
            prm.set("Minerals","Passive,Passive");
            prm.set("Volume fractions minerals","0.7,0.3");
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
      prm.leave_subsection ();
    }
    prm.leave_subsection ();
  }
  prm.leave_subsection ();

  cpo.parse_parameters(prm);
  cpo.initialize();


  cpo_elastic_tensor.parse_parameters(prm);
  //cpo_elastic_tensor.initialize();


  // All these numbers are directly from the Fortran D-Rex
  // Had to fix the random seed to get consistent results.
  // Fixed the random set to an array filled with zeros.
  using namespace dealii;
  Tensor<2,3> rotation_matrix;
  rotation_matrix[TableIndices<2>(0,0)] = -0.87492387659370430;
  rotation_matrix[TableIndices<2>(0,1)] = -0.47600252255715020;
  rotation_matrix[TableIndices<2>(0,2)] = -0.10151800968122601;
  rotation_matrix[TableIndices<2>(1,0)] = 0.13036031917262200;
  rotation_matrix[TableIndices<2>(1,1)] = -2.6406769698713056E-002;
  rotation_matrix[TableIndices<2>(1,2)] = -0.99232315682823224;
  rotation_matrix[TableIndices<2>(2,0)] = 0.46957683898444408;
  rotation_matrix[TableIndices<2>(2,1)] = -0.87974004016081919;
  rotation_matrix[TableIndices<2>(2,2)] = 8.6098835151007691E-002;
  cpo.set_rotation_matrix_grains(cpo_data_position,data_cpo,0,0,rotation_matrix);

  rotation_matrix[TableIndices<2>(0,0)] = -0.98046837857873570;
  rotation_matrix[TableIndices<2>(0,1)] = 0.19463893778994429;
  rotation_matrix[TableIndices<2>(0,2)] = -2.8239743415760400E-002;
  rotation_matrix[TableIndices<2>(1,0)] = 6.8942963409018593E-002;
  rotation_matrix[TableIndices<2>(1,1)] = 0.20565757913350766;
  rotation_matrix[TableIndices<2>(1,2)] = -0.97619251975954824;
  rotation_matrix[TableIndices<2>(2,0)] = -0.18419735979776022;
  rotation_matrix[TableIndices<2>(2,1)] = -0.95907282258851756;
  rotation_matrix[TableIndices<2>(2,2)] = -0.21505972600848655;
  cpo.set_rotation_matrix_grains(cpo_data_position,data_cpo,0,1,rotation_matrix);

  rotation_matrix[TableIndices<2>(0,0)] = -0.60475851993637786;
  rotation_matrix[TableIndices<2>(0,1)] = 0.70843624030907060;
  rotation_matrix[TableIndices<2>(0,2)] = -0.36543256811748415;
  rotation_matrix[TableIndices<2>(1,0)] = 0.14301138974031016;
  rotation_matrix[TableIndices<2>(1,1)] = -0.35406726088673490;
  rotation_matrix[TableIndices<2>(1,2)] = -0.92567249145068109;
  rotation_matrix[TableIndices<2>(2,0)] = -0.78533679283378455;
  rotation_matrix[TableIndices<2>(2,1)] = -0.61106385979914768;
  rotation_matrix[TableIndices<2>(2,2)] = 0.11279851062549377;
  cpo.set_rotation_matrix_grains(cpo_data_position,data_cpo,0,2,rotation_matrix);

  rotation_matrix[TableIndices<2>(0,0)] = 0.86791495614143876;
  rotation_matrix[TableIndices<2>(0,1)] = 0.11797028562530290;
  rotation_matrix[TableIndices<2>(0,2)] = -0.48441508819120616;
  rotation_matrix[TableIndices<2>(1,0)] = 0.29623368357270952;
  rotation_matrix[TableIndices<2>(1,1)] = 0.66075310029034506;
  rotation_matrix[TableIndices<2>(1,2)] = 0.69114249890201696;
  rotation_matrix[TableIndices<2>(2,0)] = 0.40054408016972737;
  rotation_matrix[TableIndices<2>(2,1)] = -0.74231595966649988;
  rotation_matrix[TableIndices<2>(2,2)] = 0.53869116329275069;
  cpo.set_rotation_matrix_grains(cpo_data_position,data_cpo,0,3,rotation_matrix);

  rotation_matrix[TableIndices<2>(0,0)] = 6.2142369991145308E-002;
  rotation_matrix[TableIndices<2>(0,1)] = 0.99798524450925208;
  rotation_matrix[TableIndices<2>(0,2)] = 1.4733749724082583E-002;
  rotation_matrix[TableIndices<2>(1,0)] = 0.16313304505497195;
  rotation_matrix[TableIndices<2>(1,1)] = 3.9946223138032123E-003;
  rotation_matrix[TableIndices<2>(1,2)] = -0.99141604427573704;
  rotation_matrix[TableIndices<2>(2,0)] = -0.98947772660983668;
  rotation_matrix[TableIndices<2>(2,1)] = 6.3633828841832010E-002;
  rotation_matrix[TableIndices<2>(2,2)] = -0.16253511002663851;
  cpo.set_rotation_matrix_grains(cpo_data_position,data_cpo,0,4,rotation_matrix);

  rotation_matrix[TableIndices<2>(0,0)] = 0.95811162445076037;
  rotation_matrix[TableIndices<2>(0,1)] = -0.24047036409864109;
  rotation_matrix[TableIndices<2>(0,2)] = -0.18394917461469307;
  rotation_matrix[TableIndices<2>(1,0)] = -0.26945849516984494;
  rotation_matrix[TableIndices<2>(1,1)] = -0.95216866624585605;
  rotation_matrix[TableIndices<2>(1,2)] = -0.14816174866727450;
  rotation_matrix[TableIndices<2>(2,0)] = -0.13948857525086517;
  rotation_matrix[TableIndices<2>(2,1)] = 0.18918653693665904;
  rotation_matrix[TableIndices<2>(2,2)] = -0.97685775746677384;
  cpo.set_rotation_matrix_grains(cpo_data_position,data_cpo,0,5,rotation_matrix);

  rotation_matrix[TableIndices<2>(0,0)] = 0.54599545795872684;
  rotation_matrix[TableIndices<2>(0,1)] = -0.79260950430410815;
  rotation_matrix[TableIndices<2>(0,2)] = -0.27534584625644454;
  rotation_matrix[TableIndices<2>(1,0)] = -0.83567561116582212;
  rotation_matrix[TableIndices<2>(1,1)] = -0.55085590166106368;
  rotation_matrix[TableIndices<2>(1,2)] = -6.8746495709015629E-002;
  rotation_matrix[TableIndices<2>(2,0)] = -9.6081601263635075E-002;
  rotation_matrix[TableIndices<2>(2,1)] = 0.26479508638287669;
  rotation_matrix[TableIndices<2>(2,2)] = -0.96221681521835400;
  cpo.set_rotation_matrix_grains(cpo_data_position,data_cpo,0,6,rotation_matrix);

  rotation_matrix[TableIndices<2>(0,0)] = 0.14540407652501869;
  rotation_matrix[TableIndices<2>(0,1)] = -0.61649901222701298;
  rotation_matrix[TableIndices<2>(0,2)] = -0.77881206212059828;
  rotation_matrix[TableIndices<2>(1,0)] = -0.59216315099851768;
  rotation_matrix[TableIndices<2>(1,1)] = -0.68282211925029168;
  rotation_matrix[TableIndices<2>(1,2)] = 0.43341920614609419;
  rotation_matrix[TableIndices<2>(2,0)] = -0.79819736735552915;
  rotation_matrix[TableIndices<2>(2,1)] = 0.39350172081601731;
  rotation_matrix[TableIndices<2>(2,2)] = -0.46322016633770996;
  cpo.set_rotation_matrix_grains(cpo_data_position,data_cpo,0,7,rotation_matrix);


  cpo.set_volume_fractions_grains(cpo_data_position,data_cpo,0,0,2.5128593570287589E-002);
  cpo.set_volume_fractions_grains(cpo_data_position,data_cpo,0,1,0.83128842847575013);
  cpo.set_volume_fractions_grains(cpo_data_position,data_cpo,0,2,2.4387041141724769E-002);
  cpo.set_volume_fractions_grains(cpo_data_position,data_cpo,0,3,2.4763275182773107E-002);
  cpo.set_volume_fractions_grains(cpo_data_position,data_cpo,0,4,2.4801714431754770E-002);
  cpo.set_volume_fractions_grains(cpo_data_position,data_cpo,0,5,2.3943562805875843E-002);
  cpo.set_volume_fractions_grains(cpo_data_position,data_cpo,0,6,2.1493810045792379E-002);
  cpo.set_volume_fractions_grains(cpo_data_position,data_cpo,0,7,2.4193574346041427E-002);

  rotation_matrix[TableIndices<2>(0,0)] = -0.66168933252008499;
  rotation_matrix[TableIndices<2>(0,1)] = -0.27722421136423192;
  rotation_matrix[TableIndices<2>(0,2)] = 0.70016104334335305;
  rotation_matrix[TableIndices<2>(1,0)] = -0.45052346117292291;
  rotation_matrix[TableIndices<2>(1,1)] = -0.59946225647123619;
  rotation_matrix[TableIndices<2>(1,2)] = -0.66370395454608444;
  rotation_matrix[TableIndices<2>(2,0)] = 0.60350910238307343;
  rotation_matrix[TableIndices<2>(2,1)] = -0.75116099269655345;
  rotation_matrix[TableIndices<2>(2,2)] = 0.27207473610935284;
  cpo.set_rotation_matrix_grains(cpo_data_position,data_cpo,1,0,rotation_matrix);

  rotation_matrix[TableIndices<2>(0,0)] = 0.70309122563849258;
  rotation_matrix[TableIndices<2>(0,1)] = -0.23393734574397834;
  rotation_matrix[TableIndices<2>(0,2)] = -0.67775637808568567;
  rotation_matrix[TableIndices<2>(1,0)] = 0.68298185906364617;
  rotation_matrix[TableIndices<2>(1,1)] = -6.6406501211459870E-002;
  rotation_matrix[TableIndices<2>(1,2)] = 0.73168002139436472;
  rotation_matrix[TableIndices<2>(2,0)] = -0.21654097292603913;
  rotation_matrix[TableIndices<2>(2,1)] = -0.97001835897279087;
  rotation_matrix[TableIndices<2>(2,2)] = 0.11341644508992850;
  cpo.set_rotation_matrix_grains(cpo_data_position,data_cpo,1,1,rotation_matrix);

  rotation_matrix[TableIndices<2>(0,0)] = -0.37892641930874921;
  rotation_matrix[TableIndices<2>(0,1)] = 0.92610670487847191;
  rotation_matrix[TableIndices<2>(0,2)] = 3.9338024409460159E-002;
  rotation_matrix[TableIndices<2>(1,0)] = -0.10014673867324062;
  rotation_matrix[TableIndices<2>(1,1)] = 2.1922389732184896E-004;
  rotation_matrix[TableIndices<2>(1,2)] = -0.99514251026853373;
  rotation_matrix[TableIndices<2>(2,0)] = -0.92466780075091537;
  rotation_matrix[TableIndices<2>(2,1)] = -0.37833500123197994;
  rotation_matrix[TableIndices<2>(2,2)] = 9.2539397053637271E-002;
  cpo.set_rotation_matrix_grains(cpo_data_position,data_cpo,1,2,rotation_matrix);

  rotation_matrix[TableIndices<2>(0,0)] = 0.73973803569795438;
  rotation_matrix[TableIndices<2>(0,1)] = -0.58084801380447959;
  rotation_matrix[TableIndices<2>(0,2)] = -0.35134329241205425;
  rotation_matrix[TableIndices<2>(1,0)] = 0.30883050284698427;
  rotation_matrix[TableIndices<2>(1,1)] = -0.17349072757172229;
  rotation_matrix[TableIndices<2>(1,2)] = 0.94047567596335968;
  rotation_matrix[TableIndices<2>(2,0)] = -0.60686583881969203;
  rotation_matrix[TableIndices<2>(2,1)] = -0.79564553570809793;
  rotation_matrix[TableIndices<2>(2,2)] = 5.0387432541084881E-002;
  cpo.set_rotation_matrix_grains(cpo_data_position,data_cpo,1,3,rotation_matrix);

  rotation_matrix[TableIndices<2>(0,0)] = 4.2730545437086563E-002;
  rotation_matrix[TableIndices<2>(0,1)] = 0.99790446393350407;
  rotation_matrix[TableIndices<2>(0,2)] = 4.9677348699965519E-002;
  rotation_matrix[TableIndices<2>(1,0)] = 0.24607319829629554;
  rotation_matrix[TableIndices<2>(1,1)] = 3.7594987859276820E-002;
  rotation_matrix[TableIndices<2>(1,2)] = -0.97243353311600111;
  rotation_matrix[TableIndices<2>(2,0)] = -0.97326576151708344;
  rotation_matrix[TableIndices<2>(2,1)] = 5.2990062123162360E-002;
  rotation_matrix[TableIndices<2>(2,2)] = -0.24420040893728351;
  cpo.set_rotation_matrix_grains(cpo_data_position,data_cpo,1,4,rotation_matrix);

  rotation_matrix[TableIndices<2>(0,0)] = -0.92112012021624434;
  rotation_matrix[TableIndices<2>(0,1)] = 0.17590913101287936;
  rotation_matrix[TableIndices<2>(0,2)] = 0.34726602737040008;
  rotation_matrix[TableIndices<2>(1,0)] = -0.27292179963530816;
  rotation_matrix[TableIndices<2>(1,1)] = -0.92793421776529450;
  rotation_matrix[TableIndices<2>(1,2)] = -0.25387354780599874;
  rotation_matrix[TableIndices<2>(2,0)] = 0.27758135269283285;
  rotation_matrix[TableIndices<2>(2,1)] = -0.32862450412044819;
  rotation_matrix[TableIndices<2>(2,2)] = 0.90274831467569783;
  cpo.set_rotation_matrix_grains(cpo_data_position,data_cpo,1,5,rotation_matrix);

  rotation_matrix[TableIndices<2>(0,0)] = -4.6028553384029648E-002;
  rotation_matrix[TableIndices<2>(0,1)] = -0.82828827663204008;
  rotation_matrix[TableIndices<2>(0,2)] = -0.56150509327989218;
  rotation_matrix[TableIndices<2>(1,0)] = -0.62780239844730912;
  rotation_matrix[TableIndices<2>(1,1)] = -0.41264414122407550;
  rotation_matrix[TableIndices<2>(1,2)] = 0.66416683570654200;
  rotation_matrix[TableIndices<2>(2,0)] = -0.78150657001547463;
  rotation_matrix[TableIndices<2>(2,1)] = 0.37975977805647526;
  rotation_matrix[TableIndices<2>(2,2)] = -0.50060220610535922;
  cpo.set_rotation_matrix_grains(cpo_data_position,data_cpo,1,6,rotation_matrix);

  rotation_matrix[TableIndices<2>(0,0)] = -0.69974343277254114;
  rotation_matrix[TableIndices<2>(0,1)] = -0.60514581197791628;
  rotation_matrix[TableIndices<2>(0,2)] = 0.38074455421574593;
  rotation_matrix[TableIndices<2>(1,0)] = 0.27181144918830025;
  rotation_matrix[TableIndices<2>(1,1)] = -0.71756260285634454;
  rotation_matrix[TableIndices<2>(1,2)] = -0.64160793166977359;
  rotation_matrix[TableIndices<2>(2,0)] = 0.66130789657394939;
  rotation_matrix[TableIndices<2>(2,1)] = -0.34496383827187421;
  rotation_matrix[TableIndices<2>(2,2)] = 0.66656470768691123;
  cpo.set_rotation_matrix_grains(cpo_data_position,data_cpo,1,7,rotation_matrix);


  cpo.set_volume_fractions_grains(cpo_data_position,data_cpo,1,0,2.4802805652404735E-002);
  cpo.set_volume_fractions_grains(cpo_data_position,data_cpo,1,1,2.4917064186342413E-002);
  cpo.set_volume_fractions_grains(cpo_data_position,data_cpo,1,2,2.4790722064907112E-002);
  cpo.set_volume_fractions_grains(cpo_data_position,data_cpo,1,3,2.4932751194567736E-002);
  cpo.set_volume_fractions_grains(cpo_data_position,data_cpo,1,4,2.4896744967217745E-002);
  cpo.set_volume_fractions_grains(cpo_data_position,data_cpo,1,5,0.79902139535472505);
  cpo.set_volume_fractions_grains(cpo_data_position,data_cpo,1,6,2.4861121542485400E-002);
  cpo.set_volume_fractions_grains(cpo_data_position,data_cpo,1,7,5.1777395037349808E-002);

  reference_elastic_tensor[0][0] = 282.99195951271281;
  reference_elastic_tensor[0][1] = 74.161110997372660;
  reference_elastic_tensor[0][2] = 69.528000044099443;
  reference_elastic_tensor[0][3] = 0.85449958913948032;
  reference_elastic_tensor[0][4] = 0.15156631865980030;
  reference_elastic_tensor[0][5] = -10.295196344728696;
  reference_elastic_tensor[1][0] = 74.161110997372674;
  reference_elastic_tensor[1][1] = 223.44404040361212;
  reference_elastic_tensor[1][2] = 70.938305212304968;
  reference_elastic_tensor[1][3] = 1.4935368335052783;
  reference_elastic_tensor[1][4] = -1.5555954844298270;
  reference_elastic_tensor[1][5] = -3.5461136157235025;
  reference_elastic_tensor[2][0] = 69.528000044099443;
  reference_elastic_tensor[2][1] = 70.938305212304940;
  reference_elastic_tensor[2][2] = 208.89404139751849;
  reference_elastic_tensor[2][3] = 0.48656291118743428;
  reference_elastic_tensor[2][4] = 0.49424401802786699;
  reference_elastic_tensor[2][5] = -0.15651560991351129;
  reference_elastic_tensor[3][0] = 0.85449958913948365;
  reference_elastic_tensor[3][1] = 1.4935368335052732;
  reference_elastic_tensor[3][2] = 0.48656291118742612;
  reference_elastic_tensor[3][3] = 71.502776245041289;
  reference_elastic_tensor[3][4] = -2.3939425644793477;
  reference_elastic_tensor[3][5] = 0.62714033620769472;
  reference_elastic_tensor[4][0] = 0.15156631865980935;
  reference_elastic_tensor[4][1] = -1.5555954844298292;
  reference_elastic_tensor[4][2] = 0.49424401802786488;
  reference_elastic_tensor[4][3] = -2.3939425644793459;
  reference_elastic_tensor[4][4] = 78.565442407530099;
  reference_elastic_tensor[4][5] = -0.69890743507815323;
  reference_elastic_tensor[5][0] = -10.295196344728710;
  reference_elastic_tensor[5][1] = -3.5461136157235131;
  reference_elastic_tensor[5][2] = -0.15651560991350758;
  reference_elastic_tensor[5][3] = 0.62714033620769460;
  reference_elastic_tensor[5][4] = -0.69890743507815523;
  reference_elastic_tensor[5][5] = 80.599981331604567;


  cpo.set_volume_fraction_mineral(cpo_data_position,data_cpo,0,0.7);
  cpo.set_volume_fraction_mineral(cpo_data_position,data_cpo,1,0.3);

  cpo.set_deformation_type(cpo_data_position,data_cpo,0,(double)aspect::Particle::Property::DeformationType::olivine_a_fabric);
  cpo.set_deformation_type(cpo_data_position,data_cpo,1,(double)aspect::Particle::Property::DeformationType::enstatite);

  computed_elastic_tensor = cpo_elastic_tensor.voigt_average_elastic_tensor(cpo,
                                                                            cpo_data_position,
                                                                            data_cpo);


  for (size_t i = 0; i < 6; i++)
    {
      for (size_t j = 0; j < 6; j++)
        {
          CHECK(computed_elastic_tensor[i][j] == Approx(reference_elastic_tensor[i][j]));
        }
    }

  // test store and load functions
  // the first and last element should not be changed
  // by these functions.
  std::vector<double> array_ref = {0.0,
                                   1.,2.,3.,4.,5,6,
                                   7.,8.,9.,10,11,12,
                                   13,14,15,16,17,18,
                                   19,20,21,22,23,24,
                                   25,26,27,28,29,30,
                                   31,32,33,34,35,36,
                                   37
                                  };

  std::vector<double> array = {0.0,
                               1.,2.,3.,4.,5.,6.,
                               7.,8.,9.,10,11,12,
                               13,14,15,16,17,18,
                               19,20,21,22,23,24,
                               25,26,27,28,29,30,
                               31,32,33,34,35,36,
                               37
                              };

  // There used be be 36 unique entries, but now because we are using the
  // symmetric tensor, there are only 21 unique entries.
  std::vector<double> array_plus_100 = {0.0,
                                        101.,102.,103.,104.,105.,106.,
                                        107.,108.,109.,110.,111.,112.,
                                        113.,114.,115.,116.,117.,118.,
                                        119.,120.,121.,22.,23.,24.,
                                        25.,26.,27.,28.,29.,30.,
                                        31.,32.,33.,34.,35.,36.,
                                        37.
                                       };

  cpo_data_position = 1;
  dealii::ArrayView<double> data(&array[0],38);
  dealii::SymmetricTensor<2,6> tensor = cpo_elastic_tensor.get_elastic_tensor(cpo_data_position,data);

  for (unsigned int i = 0; i < dealii::SymmetricTensor<2,6>::n_independent_components ; ++i)
    {
      CHECK(data[cpo_data_position + i] == tensor[dealii::SymmetricTensor<2,6>::unrolled_to_component_indices(i)]);
    }

  cpo_elastic_tensor.set_elastic_tensor(cpo_data_position,data,tensor);

  for (unsigned int i = 0; i < array.size() ; ++i)
    CHECK(data[i] == array_ref[i]);

  for (unsigned int i = 0; i < dealii::SymmetricTensor<2,6>::n_independent_components ; ++i)
    {
      tensor[dealii::SymmetricTensor<2,6>::unrolled_to_component_indices(i)] += 100;
    }

  cpo_elastic_tensor.set_elastic_tensor(cpo_data_position,data,tensor);

  for (unsigned int i = 0; i < array.size() ; ++i)
    CHECK(data[i] == array_plus_100[i]);

  tensor = cpo_elastic_tensor.get_elastic_tensor(cpo_data_position,data);

  for (unsigned int i = 0; i < dealii::SymmetricTensor<2,6>::n_independent_components ; ++i)
    {
      CHECK(data[cpo_data_position + i] == tensor[dealii::SymmetricTensor<2,6>::unrolled_to_component_indices(i)]);
    }

  for (unsigned int i = 0; i < dealii::SymmetricTensor<2,6>::n_independent_components ; ++i)
    {
      CHECK(array_plus_100[cpo_data_position + i] == tensor[dealii::SymmetricTensor<2,6>::unrolled_to_component_indices(i)]);
    }
}
