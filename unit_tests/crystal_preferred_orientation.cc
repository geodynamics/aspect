/*
  Copyright (C) 2022 by the authors of the ASPECT code.

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

  std::vector<unsigned int> deformation_types_ref = {(unsigned int)aspect::Particle::Property::DeformationType::passive,
                                                     (unsigned int)aspect::Particle::Property::DeformationType::passive
                                                    };

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
