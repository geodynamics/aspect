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
          ((computed[0][0] == Approx(expected[0][0]) || std::fabs(computed[0][0] < tol))
           && (computed[0][1] == Approx(expected[0][1]) || std::fabs(computed[0][1] < tol))
           && (computed[0][2] == Approx(expected[0][2]) || std::fabs(computed[0][2] < tol))
           && (computed[1][0] == Approx(expected[1][0]) || std::fabs(computed[1][0] < tol))
           && (computed[1][1] == Approx(expected[1][1]) || std::fabs(computed[1][1] < tol))
           && (computed[1][2] == Approx(expected[1][2]) || std::fabs(computed[1][2] < tol))
           && (computed[2][0] == Approx(expected[2][0]) || std::fabs(computed[2][0] < tol))
           && (computed[2][1] == Approx(expected[2][1]) || std::fabs(computed[2][1] < tol))
           && (computed[2][2] == Approx(expected[2][2]) || std::fabs(computed[2][2] < tol)))
          ||
          ((computed[0][0] == Approx(-expected[0][0]) || std::fabs(computed[0][0] < tol))
           && (computed[0][1] == Approx(-expected[0][1]) || std::fabs(computed[0][1] < tol))
           && (computed[0][2] == Approx(-expected[0][2]) || std::fabs(computed[0][2] < tol))
           && (computed[1][0] == Approx(-expected[1][0]) || std::fabs(computed[1][0] < tol))
           && (computed[1][1] == Approx(-expected[1][1]) || std::fabs(computed[1][1] < tol))
           && (computed[1][2] == Approx(-expected[1][2]) || std::fabs(computed[1][2] < tol))
           && (computed[2][0] == Approx(-expected[2][0]) || std::fabs(computed[2][0] < tol))
           && (computed[2][1] == Approx(-expected[2][1]) || std::fabs(computed[2][1] < tol))
           && (computed[2][2] == Approx(-expected[2][2]) || std::fabs(computed[2][2] < tol)))));
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

  unsigned int cpo_data_position = 1;
  std::vector<double> data_array(70,-1.);
  dealii::ArrayView<double> data(&data_array[0],70);

  cpo.ref_volume_fraction_mineral(cpo_data_position,data,0) = 0.7;
  cpo.ref_volume_fraction_mineral(cpo_data_position,data,1) = 0.3;

  cpo.ref_rotation_matrix_grains_indices(cpo_data_position,data,0,0,0,0) = 0;
  cpo.ref_rotation_matrix_grains_indices(cpo_data_position,data,0,0,0,1) = 1./1000.;
  cpo.ref_rotation_matrix_grains_indices(cpo_data_position,data,0,0,0,2) = 2./1000.;
  cpo.ref_rotation_matrix_grains_indices(cpo_data_position,data,0,0,1,0) = 3./1000.;
  cpo.ref_rotation_matrix_grains_indices(cpo_data_position,data,0,0,1,1) = 4./1000.;
  cpo.ref_rotation_matrix_grains_indices(cpo_data_position,data,0,0,1,2) = 5./1000.;
  cpo.ref_rotation_matrix_grains_indices(cpo_data_position,data,0,0,2,0) = 6./1000.;
  cpo.ref_rotation_matrix_grains_indices(cpo_data_position,data,0,0,2,1) = 7./1000.;
  cpo.ref_rotation_matrix_grains_indices(cpo_data_position,data,0,0,2,2) = 8./1000.;
  cpo.ref_rotation_matrix_grains_indices(cpo_data_position,data,0,1,0,0) = 9./1000.;
  cpo.ref_rotation_matrix_grains_indices(cpo_data_position,data,0,1,0,1) = 10./1000.;
  cpo.ref_rotation_matrix_grains_indices(cpo_data_position,data,0,1,0,2) = 11./1000.;
  cpo.ref_rotation_matrix_grains_indices(cpo_data_position,data,0,1,1,0) = 12./1000.;
  cpo.ref_rotation_matrix_grains_indices(cpo_data_position,data,0,1,1,1) = 13./1000.;
  cpo.ref_rotation_matrix_grains_indices(cpo_data_position,data,0,1,1,2) = 14./1000.;
  cpo.ref_rotation_matrix_grains_indices(cpo_data_position,data,0,1,2,0) = 15./1000.;
  cpo.ref_rotation_matrix_grains_indices(cpo_data_position,data,0,1,2,1) = 16./1000.;
  cpo.ref_rotation_matrix_grains_indices(cpo_data_position,data,0,1,2,2) = 17./1000.;
  cpo.ref_rotation_matrix_grains_indices(cpo_data_position,data,0,2,0,0) = 18./1000.;
  cpo.ref_rotation_matrix_grains_indices(cpo_data_position,data,0,2,0,1) = 19./1000.;
  cpo.ref_rotation_matrix_grains_indices(cpo_data_position,data,0,2,0,2) = 20./1000.;
  cpo.ref_rotation_matrix_grains_indices(cpo_data_position,data,0,2,1,0) = 21./1000.;
  cpo.ref_rotation_matrix_grains_indices(cpo_data_position,data,0,2,1,1) = 22./1000.;
  cpo.ref_rotation_matrix_grains_indices(cpo_data_position,data,0,2,1,2) = 23./1000.;
  cpo.ref_rotation_matrix_grains_indices(cpo_data_position,data,0,2,2,0) = 24./1000.;
  cpo.ref_rotation_matrix_grains_indices(cpo_data_position,data,0,2,2,1) = 25./1000.;
  cpo.ref_rotation_matrix_grains_indices(cpo_data_position,data,0,2,2,2) = 26./1000.;

  cpo.ref_volume_fractions_grains(cpo_data_position,data,0,0) = 0.1;
  cpo.ref_volume_fractions_grains(cpo_data_position,data,0,1) = 0.2;
  cpo.ref_volume_fractions_grains(cpo_data_position,data,0,2) = 0.3;

  cpo.ref_rotation_matrix_grains_indices(cpo_data_position,data,1,0,0,0) = 27./1000.;
  cpo.ref_rotation_matrix_grains_indices(cpo_data_position,data,1,0,0,1) = 28./1000.;
  cpo.ref_rotation_matrix_grains_indices(cpo_data_position,data,1,0,0,2) = 29./1000.;
  cpo.ref_rotation_matrix_grains_indices(cpo_data_position,data,1,0,1,0) = 30./1000.;
  cpo.ref_rotation_matrix_grains_indices(cpo_data_position,data,1,0,1,1) = 31./1000.;
  cpo.ref_rotation_matrix_grains_indices(cpo_data_position,data,1,0,1,2) = 32./1000.;
  cpo.ref_rotation_matrix_grains_indices(cpo_data_position,data,1,0,2,0) = 33./1000.;
  cpo.ref_rotation_matrix_grains_indices(cpo_data_position,data,1,0,2,1) = 34./1000.;
  cpo.ref_rotation_matrix_grains_indices(cpo_data_position,data,1,0,2,2) = 35./1000.;
  cpo.ref_rotation_matrix_grains_indices(cpo_data_position,data,1,1,0,0) = 36./1000.;
  cpo.ref_rotation_matrix_grains_indices(cpo_data_position,data,1,1,0,1) = 37./1000.;
  cpo.ref_rotation_matrix_grains_indices(cpo_data_position,data,1,1,0,2) = 38./1000.;
  cpo.ref_rotation_matrix_grains_indices(cpo_data_position,data,1,1,1,0) = 39./1000.;
  cpo.ref_rotation_matrix_grains_indices(cpo_data_position,data,1,1,1,1) = 40./1000.;
  cpo.ref_rotation_matrix_grains_indices(cpo_data_position,data,1,1,1,2) = 41./1000.;
  cpo.ref_rotation_matrix_grains_indices(cpo_data_position,data,1,1,2,0) = 42./1000.;
  cpo.ref_rotation_matrix_grains_indices(cpo_data_position,data,1,1,2,1) = 43./1000.;
  cpo.ref_rotation_matrix_grains_indices(cpo_data_position,data,1,1,2,2) = 44./1000.;
  cpo.ref_rotation_matrix_grains_indices(cpo_data_position,data,1,2,0,0) = 45./1000.;
  cpo.ref_rotation_matrix_grains_indices(cpo_data_position,data,1,2,0,1) = 46./1000.;
  cpo.ref_rotation_matrix_grains_indices(cpo_data_position,data,1,2,0,2) = 47./1000.;
  cpo.ref_rotation_matrix_grains_indices(cpo_data_position,data,1,2,1,0) = 48./1000.;
  cpo.ref_rotation_matrix_grains_indices(cpo_data_position,data,1,2,1,1) = 49./1000.;
  cpo.ref_rotation_matrix_grains_indices(cpo_data_position,data,1,2,1,2) = 50./1000.;
  cpo.ref_rotation_matrix_grains_indices(cpo_data_position,data,1,2,2,0) = 51./1000.;
  cpo.ref_rotation_matrix_grains_indices(cpo_data_position,data,1,2,2,1) = 52./1000.;
  cpo.ref_rotation_matrix_grains_indices(cpo_data_position,data,1,2,2,2) = 53./1000.;

  cpo.ref_volume_fractions_grains(cpo_data_position,data,1,0) = 0.4;
  cpo.ref_volume_fractions_grains(cpo_data_position,data,1,1) = 0.5;
  cpo.ref_volume_fractions_grains(cpo_data_position,data,1,2) = 0.6;

  std::vector<unsigned int> deformation_types_ref = {(unsigned int)aspect::Particle::Property::DeformationType::passive,
                                                     (unsigned int)aspect::Particle::Property::DeformationType::passive
                                                    };

  data[0] = 20847932.2;
  data[65] = 6541684.3;

  cpo.ref_deformation_type(cpo_data_position,data,0) = (double)aspect::Particle::Property::DeformationType::passive;
  cpo.ref_deformation_type(cpo_data_position,data,1) = (double)aspect::Particle::Property::DeformationType::passive;


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
            prm.enter_subsection("Uniform grains and random uniform rotations");
            {
              // Let the minerals just passively rotate with the rotation of
              // the particle caused by the flow.
              prm.set("Minerals","Passive,Passive");
              prm.set("Volume fractions minerals","0.5,0.5");
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
    CHECK(data[3] == Approx(0.159063));
    CHECK(data[4] == Approx(-0.11941));
    CHECK(data[5] == Approx(0.9800204275));
    CHECK(data[6] == Approx(-0.0888556));
    CHECK(data[7] == Approx(-0.990362));
    CHECK(data[8] == Approx(-0.1062486256));
    CHECK(data[9] == Approx(0.983261702));
    CHECK(data[10] == Approx(-0.0701800114));
    CHECK(data[11] == Approx(-0.1681403917));
    CHECK(data[12] == Approx(0.2));
    CHECK(data[13] == Approx(0.4095335744));
    CHECK(data[14] == Approx(-0.3401753011));
    CHECK(data[15] == Approx(0.8465004524));
    CHECK(data[16] == Approx(0.7605716382));
    CHECK(data[17] == Approx(0.639714977));
    CHECK(data[18] == Approx(-0.1108852174));
    CHECK(data[19] == Approx(-0.5037986052));
    CHECK(data[20] == Approx(0.6892354553));
    CHECK(data[21] == Approx(0.5207124471));
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
    CHECK(data[103] == Approx(-0.816855));

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

    Tensor<2,3> velocity_gradient_tensor_nondimensional;
    velocity_gradient_tensor_nondimensional[0][1] = 2.0* 0.5959;
    velocity_gradient_tensor_nondimensional[1][0] = 2.0* 0.5959;

    std::array<double,4> ref_resolved_shear_stress;
    ref_resolved_shear_stress[0] = 1;
    ref_resolved_shear_stress[1] = 2;
    ref_resolved_shear_stress[2] = 3;
    ref_resolved_shear_stress[3] = 1e60; // can't really use numerical limits max or infinite, because need to be able to square it without becoming infinite. This is the value fortran D-Rex uses.

    std::pair<std::vector<double>, std::vector<Tensor<2,3>>> derivatives;
    derivatives = cpo_3d.compute_derivatives(0, data,0,
                                             strain_rate_nondimensional, velocity_gradient_tensor_nondimensional,
                                             ref_resolved_shear_stress);

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
