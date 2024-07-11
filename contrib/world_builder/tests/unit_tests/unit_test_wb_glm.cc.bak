/*
  Copyright (C) 2020 by the authors of the World Builder code.

  This file is part of the World Builder.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published
   by the Free Software Foundation, either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#define DOCTEST_CONFIG_SUPER_FAST_ASSERTS

#include "doctest/doctest.h"

#include "glm/glm.h"

#include "world_builder/coordinate_system.h"
#include "world_builder/utilities.h"

#include <array>

using namespace WorldBuilder;
using doctest::Approx;


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
 * Compare two rotation matrices
 */
inline void compare_rotation_matrices_approx(
  const std::array<std::array<double,3>,3> &computed,
  const std::array<std::array<double,3>,3> &expected)
{
  // sign of eigenvector is not important
  INFO("rotation matrices are not the same: \n" <<
       "expected = " << expected[0][0] << ' ' << expected[0][1] << ' ' << expected[0][2] << "\n" <<
       "           " << expected[1][0] << ' ' << expected[1][1] << ' ' << expected[1][2] << "\n" <<
       "           " << expected[2][0] << ' ' << expected[2][1] << ' ' << expected[2][2] << "\n" <<
       "computed = " << computed[0][0] << ' ' << computed[0][1] << ' ' << computed[0][2] << "\n" <<
       "           " << computed[1][0] << ' ' << computed[1][1] << ' ' << computed[1][2] << "\n" <<
       "           " << computed[2][0] << ' ' << computed[2][1] << ' ' << computed[2][2] << "\n" );
  CHECK((
          (computed[0][0] == Approx(expected[0][0]) && computed[0][1] == Approx(expected[0][1]) && computed[0][2] == Approx(expected[0][2]) &&
           computed[1][0] == Approx(expected[1][0]) && computed[1][1] == Approx(expected[1][1]) && computed[1][2] == Approx(expected[1][2]) &&
           computed[2][0] == Approx(expected[2][0]) && computed[2][1] == Approx(expected[2][1]) && computed[2][2] == Approx(expected[2][2]))
          ||
          (computed[0][0] == Approx(-expected[0][0]) && computed[0][1] == Approx(-expected[0][1]) && computed[0][2] == Approx(-expected[0][2]) &&
           computed[1][0] == Approx(-expected[1][0]) && computed[1][1] == Approx(-expected[1][1]) && computed[1][2] == Approx(-expected[1][2]) &&
           computed[2][0] == Approx(-expected[2][0]) && computed[2][1] == Approx(-expected[2][1]) && computed[2][2] == Approx(-expected[2][2]))));
}

inline
void
compare_quaternions(
  const glm::quaternion::quat &computed,
  const glm::quaternion::quat &expected
)
{
  CHECK(computed.w == Approx(expected.w));
  CHECK(computed.x == Approx(expected.x));
  CHECK(computed.y == Approx(expected.y));
  CHECK(computed.z == Approx(expected.z));
}

TEST_CASE("glm quat basic functions")
{
  typedef glm::quaternion::quat QT;
  QT q1(2,3,4,5);
  QT const q2(6,7,8,9);
  QT const q3(9,8,7,6);


  QT q_res = -q1;
  compare_quaternions(q_res, QT(-2,-3,-4,-5));

  q_res = q1 + q2;
  compare_quaternions(q_res, QT(8,10,12,14));

  q_res = q2 - q1;
  compare_quaternions(q_res, QT(-4,-4,-4,-4));
  q_res = q3 - q1;
  compare_quaternions(q_res, QT(-7,-5,-3,-1));


  q_res = q1 * 2.0;
  compare_quaternions(q_res, QT(4,6,8,10));
  q_res = 2.0 * q1;
  compare_quaternions(q_res, QT(4,6,8,10));


  q_res = (2.0 * q1) / 2.0;
  compare_quaternions(q_res, q1);

  CHECK(glm::quaternion::dot(q1,q2) ==  Approx(110.0));
  CHECK(glm::quaternion::dot(q2,q3) ==  Approx(220.0));
  CHECK(glm::quaternion::dot(q2,-q3) ==  Approx(-220.0));
  CHECK(glm::quaternion::dot(-q2,-q3) ==  Approx(220.0));

  double res = glm::quaternion::mix(0.0,1.0,0.0);
  CHECK(res == Approx(0.0));

  res = glm::quaternion::mix(0.0,1.0,0.25);
  CHECK(res == Approx(Approx(0.25)));

  res = glm::quaternion::mix(0.0,1.0,0.5);
  CHECK(res == Approx(0.5));

  res = glm::quaternion::mix(0.0,1.0,0.75);
  CHECK(res == Approx(0.75));

  res = glm::quaternion::mix(0.0,1.0,1.0);
  CHECK(res == Approx(1.0));


  res = glm::quaternion::mix(1.0,0.0,0.0);
  CHECK(res == Approx(1.0));

  res = glm::quaternion::mix(1.0,0.0,0.25);
  CHECK(res == Approx(0.75));

  res = glm::quaternion::mix(1.0,0.0,0.5);
  CHECK(res == Approx(0.5));

  res = glm::quaternion::mix(1.0,0.0,0.75);
  CHECK(res == Approx(0.25));

  res = glm::quaternion::mix(1.0,0.0,1.0);
  CHECK(res == Approx(0.0));


  res = glm::quaternion::mix(-3.0,1.0,0.0);
  CHECK(res == Approx(-3.0));

  res = glm::quaternion::mix(-3.0,1.0,0.25);
  CHECK(res == Approx(-2.0));

  res = glm::quaternion::mix(-3.0,1.0,0.5);
  CHECK(res == Approx(-1.0));

  res = glm::quaternion::mix(-3.0,1.0,0.75);
  CHECK(res == Approx(0.0));

  res = glm::quaternion::mix(-3.0,1.0,1.0);
  CHECK(res == Approx(1.0));

}

TEST_CASE("glm quat slerp functions")
{
  {
    // This shows that two equal rotation matrices average in the same rotation matrix
    const std::array<std::array<double,3>,3> rot1 = {{{{0.36,0.48,-0.8}},{{-0.8,0.6,0}}, {{0.48,0.64, 0.6}}}}; //{{normalize({{1,2,3}}),normalize({{4,5,6}}),normalize({{7,8,9}})}};
    const glm::quaternion::quat quat_current = glm::quaternion::quat_cast(rot1);
    const glm::quaternion::quat quat_next = glm::quaternion::quat_cast(rot1);

    const glm::quaternion::quat quat_average = glm::quaternion::slerp(quat_current,quat_next,0.5);

    const std::array<std::array<double,3>,3> result1 = glm::quaternion::mat3_cast(quat_average);

    auto ea1 = Utilities::euler_angles_from_rotation_matrix(rot1);
    auto ear = Utilities::euler_angles_from_rotation_matrix(result1);

    compare_rotation_matrices_approx(result1, rot1);
    compare_3d_arrays_approx(ea1,ear);
  }

  {
    // This shows that a very small difference in the rotation matrix causes a very small difference in the resulting rotation matrix.
    const std::array<std::array<double,3>,3> rot1 = {{{{0.360000,0.480000,-0.800000}},{{-0.800000,0.600000,0.000000000}}, {{0.480000,0.640000, 0.600000}}}};
    const std::array<std::array<double,3>,3> rot2 = {{{{0.350000,0.480000,-0.800000}},{{-0.800000,0.600000,0.000000000}}, {{0.480000,0.640000, 0.600000}}}};
    const std::array<std::array<double,3>,3> expt = {{{{0.358571,0.479818,-0.800533}},{{-0.800533,0.599107,0.000627009}}, {{0.479818,0.640802, 0.599107}}}};

    const glm::quaternion::quat quat_current = glm::quaternion::quat_cast(rot1);
    const glm::quaternion::quat quat_next = glm::quaternion::quat_cast(rot2);

    const glm::quaternion::quat quat_average = glm::quaternion::slerp(quat_current,quat_next,0.5);

    const std::array<std::array<double,3>,3> result1 = glm::quaternion::mat3_cast(quat_average);

    compare_rotation_matrices_approx(result1, expt);
  }

  {
    // This shows that a small difference in the rotation matrix causes a small difference in the resulting rotation matrix.
    const std::array<std::array<double,3>,3> rot1 = {{{{0.360000,0.480000,-0.800000}},{{-0.800000,0.600000,0.000000000}}, {{0.480000,0.640000, 0.600000}}}};
    const std::array<std::array<double,3>,3> rot2 = {{{{0.300000,0.480000,-0.800000}},{{-0.800000,0.600000,0.000000000}}, {{0.480000,0.640000, 0.600000}}}};
    const std::array<std::array<double,3>,3> expt = {{{{0.351289,0.478886,-0.803242}},{{-0.803242,0.594555,0.00382358}}, {{0.478886,0.644888, 0.594555}}}};

    const glm::quaternion::quat quat_current = glm::quaternion::quat_cast(rot1);
    const glm::quaternion::quat quat_next = glm::quaternion::quat_cast(rot2);

    const glm::quaternion::quat quat_average = glm::quaternion::slerp(quat_current,quat_next,0.5);

    const std::array<std::array<double,3>,3> result1 = glm::quaternion::mat3_cast(quat_average);

    compare_rotation_matrices_approx(result1, expt);
  }

  {
    // This shows that a small difference in the rotation matrix causes a small difference in the resulting rotation matrix.
    const std::array<std::array<double,3>,3> rot1 = {{{{0.360000,0.480000,-0.800000}},{{-0.800000,0.600000,0.000000000}}, {{0.480000,0.640000, 0.600000}}}};
    const std::array<std::array<double,3>,3> rot2 = {{{{0.300000,0.450000,-0.800000}},{{-0.800000,0.600000,0.000000000}}, {{0.480000,0.640000, 0.600000}}}};
    const std::array<std::array<double,3>,3> expt = {{{{0.35767,0.472227,-0.802856}},{{-0.7972,0.6014,0.000000000}}, {{0.481714,0.642285, 0.593783}}}};

    const glm::quaternion::quat quat_current = glm::quaternion::quat_cast(rot1);
    const glm::quaternion::quat quat_next = glm::quaternion::quat_cast(rot2);

    const glm::quaternion::quat quat_average = glm::quaternion::slerp(quat_current,quat_next,0.5);

    const std::array<std::array<double,3>,3> result1 = glm::quaternion::mat3_cast(quat_average);

    compare_rotation_matrices_approx(result1, expt);
  }

  {
    // This shows that a significant difference in the rotation matrix causes a significant difference in the resulting rotation matrix.
    const std::array<std::array<double,3>,3> rot1 = {{{{0.360000,0.480000,-0.800000}},{{-0.800000,0.600000,0.000000000}}, {{0.480000,0.640000, 0.600000}}}};
    const std::array<std::array<double,3>,3> rot2 = {{{{0.200000,0.250000,-0.400000}},{{-0.500000,0.200000,0.300000000}}, {{0.380000,0.840000, 0.500000}}}};
    const std::array<std::array<double,3>,3> expt = {{{{0.493495,0.403007,-0.707894}},{{-0.70133,0.664351,-0.0720395}}, {{0.414325,0.578479, 0.656215}}}};

    const glm::quaternion::quat quat_current = glm::quaternion::quat_cast(rot1);
    const glm::quaternion::quat quat_next = glm::quaternion::quat_cast(rot2);

    const glm::quaternion::quat quat_average = glm::quaternion::slerp(quat_current,quat_next,0.5);

    const std::array<std::array<double,3>,3> result1 = glm::quaternion::mat3_cast(quat_average);

    compare_rotation_matrices_approx(result1, expt);
  }
}