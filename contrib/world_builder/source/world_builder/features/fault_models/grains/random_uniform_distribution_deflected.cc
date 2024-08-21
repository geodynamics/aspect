/*
  Copyright (C) 2018-2024 by the authors of the World Builder code.

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

#include "world_builder/features/fault_models/grains/random_uniform_distribution_deflected.h"

#include "world_builder/nan.h"
#include "world_builder/types/array.h"
#include "world_builder/types/bool.h"
#include "world_builder/types/double.h"
#include "world_builder/types/object.h"
#include "world_builder/types/unsigned_int.h"
#include "world_builder/utilities.h"
#include "world_builder/world.h"

namespace WorldBuilder
{
  using namespace Utilities;

  namespace Features
  {
    namespace FaultModels
    {
      namespace Grains
      {
        RandomUniformDistributionDeflected::RandomUniformDistributionDeflected(WorldBuilder::World *world_)
          :
          min_depth(NaN::DSNAN),
          max_depth(NaN::DSNAN)

        {
          this->world = world_;
          this->name = "random uniform distribution deflected";
        }

        RandomUniformDistributionDeflected::~RandomUniformDistributionDeflected()
          = default;

        void
        RandomUniformDistributionDeflected::declare_entries(Parameters &prm, const std::string & /*unused*/)
        {
          // Document plugin and require entries if needed.
          // Add compositions to the required parameters.
          prm.declare_entry("", Types::Object({"compositions"}),
                            "Random uniform distribution grains model. The size of the grains can be independently set "
                            "to a single value or to a random distribution.");

          // Declare entries of this plugin
          prm.declare_entry("min distance fault center", Types::Double(0),
                            "The distance from the fault center in meters from which the composition of this feature is present.");
          prm.declare_entry("max distance fault center", Types::Double(std::numeric_limits<double>::max()),
                            "The distance from the fault in meters to which the composition of this feature is present.");

          prm.declare_entry("compositions", Types::Array(Types::UnsignedInt(),0),
                            "A list with the integer labels of the composition which are present there.");

          prm.declare_entry("orientation operation", Types::String("replace", std::vector<std::string> {"replace"}),
                            "Whether the value should replace any value previously defined at this location (replace) or "
                            "add the value to the previously define value (add, not implemented). Replacing implies that all values not "
                            "explicitly defined are set to zero.");

          prm.declare_entry("grain sizes",
                            Types::Array(Types::Double(1),0),
                            "A list of the size of all of the grains in each composition. If set to <0, the size will be randomized between 0 and 1.");

          prm.declare_entry("normalize grain sizes",
                            Types::Array(Types::Bool(true),0),
                            "A list of whether the sizes of the grains should be normalized or not. If normalized, the total of the grains of a composition will be equal to 1.");

          prm.declare_entry("deflections",
                            Types::Array(Types::Double(1),0),
                            "A list of the deflections of all of the grains in each composition between 0 and 1.");

          prm.declare_entry("basis rotation matrices", Types::Array(Types::Array(Types::Array(Types::Double(0),3,3),3,3),0),
                            "A list with the rotation matrices of the grains which are present there for each compositions.");

          prm.declare_entry("basis Euler angles z-x-z", Types::Array(Types::Array(Types::Double(0),3,3),0),
                            "A list with the z-x-z Euler angles of the grains which are present there for each compositions.");



        }

        void
        RandomUniformDistributionDeflected::parse_entries(Parameters &prm)
        {
          min_depth = prm.get<double>("min distance fault center");
          max_depth = prm.get<double>("max distance fault center");
          compositions = prm.get_vector<unsigned int>("compositions");

          const bool set_euler_angles = prm.check_entry("basis Euler angles z-x-z");
          const bool set_rotation_matrices = prm.check_entry("basis rotation matrices");

          WBAssertThrow(!(set_euler_angles == true && set_rotation_matrices == true),
                        "Only Euler angles or Rotation matrices may be set, but both are set for " << prm.get_full_json_path());


          WBAssertThrow(!(set_euler_angles == false && set_rotation_matrices == false),
                        "Euler angles or Rotation matrices have to be set, but neither are set for " << prm.get_full_json_path());

          if (set_euler_angles)
            {
              std::vector<std::array<double,3> > basis_euler_angles_vector = prm.get_vector<std::array<double,3> >("basis Euler angles z-x-z");
              basis_rotation_matrices.resize(basis_euler_angles_vector.size());
              for (size_t i = 0; i<basis_euler_angles_vector.size(); ++i)
                {
                  basis_rotation_matrices[i] = Utilities::euler_angles_to_rotation_matrix(basis_euler_angles_vector[i][0],basis_euler_angles_vector[i][1],basis_euler_angles_vector[i][2]);
                }

            }
          else
            {
              basis_rotation_matrices = prm.get_vector<std::array<std::array<double,3>,3> >("basis rotation matrices");
            }

          operation = prm.get<std::string>("orientation operation");
          grain_sizes = prm.get_vector<double>("grain sizes");
          normalize_grain_sizes = prm.get_vector<bool>("normalize grain sizes");
          deflections = prm.get_vector<double>("deflections");

          WBAssertThrow(compositions.size() == grain_sizes.size(),
                        "There are not the same amount of compositions (" << compositions.size()
                        << ") and grain_sizes (" << grain_sizes.size() << ").");
          WBAssertThrow(compositions.size() == normalize_grain_sizes.size(),
                        "There are not the same amount of compositions (" << compositions.size()
                        << ") and normalize_grain_sizes (" << normalize_grain_sizes.size() << ").");
          WBAssertThrow(compositions.size() == deflections.size(),
                        "There are not the same amount of compositions (" << compositions.size()
                        << ") and deflections (" << deflections.size() << ").");
          WBAssertThrow(compositions.size() == basis_rotation_matrices.size(),
                        "There are not the same amount of compositions (" << compositions.size()
                        << ") and rotation_matrices (" << basis_rotation_matrices.size() << ").");
        }

        WorldBuilder::grains
        RandomUniformDistributionDeflected::get_grains(const Point<3> & /*position_in_cartesian_coordinates*/,
                                                       const double /*depth*/,
                                                       const unsigned int composition_number,
                                                       WorldBuilder::grains grains_,
                                                       const double /*feature_min_depth*/,
                                                       const double /*feature_max_depth*/,
                                                       const WorldBuilder::Utilities::PointDistanceFromCurvedPlanes &distance_from_planes,
                                                       const AdditionalParameters & /*additional_parameters*/) const
        {
          WorldBuilder::grains  grains_local = grains_;
          if (std::fabs(distance_from_planes.distance_from_plane) <= max_depth && std::fabs(distance_from_planes.distance_from_plane) >= min_depth)
            {
              for (unsigned int i =0; i < compositions.size(); ++i)
                {
                  if (compositions[i] == composition_number)
                    {
                      std::uniform_real_distribution<> dist(0.0,1.0);
                      for (auto &&it_rotation_matrices : grains_local.rotation_matrices)
                        {
                          // set a uniform random a_cosine_matrix per grain
                          // This function is based on an article in Graphic Gems III, written by James Arvo, Cornell University (p 116-120).
                          // The original code can be found on  http://www.realtimerendering.com/resources/GraphicsGems/gemsiii/rand_rotation.c
                          // and is licenend accourding to this website with the following licence:
                          //
                          // "The Graphics Gems code is copyright-protected. In other words, you cannot claim the text of the code as your own and
                          // resell it. Using the code is permitted in any program, product, or library, non-commercial or commercial. Giving credit
                          // is not required, though is a nice gesture. The code comes as-is, and if there are any flaws or problems with any Gems
                          // code, nobody involved with Gems - authors, editors, publishers, or webmasters - are to be held responsible. Basically,
                          // don't be a jerk, and remember that anything free comes with no guarantee.""
                          //
                          // The book saids in the preface the following: "As in the first two volumes, all of the C and C++ code in this book is in
                          // the public domain, and is yours to study, modify, and use."

                          // first generate three random numbers between 0 and 1 and multiply them with 2 PI or 2 for z. Note that these are not the same as phi_1, theta and phi_2.
                          const double one = dist(world->get_random_number_engine());
                          const double two = dist(world->get_random_number_engine());
                          const double three = dist(world->get_random_number_engine());

                          const double theta = 2.0 * Consts::PI * one * deflections[i]; // Rotation about the pole (Z).
                          const double phi = 2.0 * Consts::PI * two; // For direction of pole deflection.
                          const double z = 2.0* three * deflections[i]; //For magnitude of pole deflection.

                          // Compute a vector V used for distributing points over the sphere
                          // via the reflection I - V Transpose(V).  This formulation of V
                          // will guarantee that if x[1] and x[2] are uniformly distributed,
                          // the reflected points will be uniform on the sphere.  Note that V
                          // has length sqrt(2) to eliminate the 2 in the Householder matrix.

                          const double r  = std::sqrt( z );
                          const double Vx = std::sin( phi ) * r;
                          const double Vy = std::cos( phi ) * r;
                          const double Vz = std::sqrt( 2.F - z );

                          // Compute the row vector S = Transpose(V) * R, where R is a simple
                          // rotation by theta about the z-axis.  No need to compute Sz since
                          // it's just Vz.

                          const double st = std::sin( theta );
                          const double ct = std::cos( theta );
                          const double Sx = Vx * ct - Vy * st;
                          const double Sy = Vx * st + Vy * ct;

                          // Construct the rotation matrix  ( V Transpose(V) - I ) R, which
                          // is equivalent to V S - R.

                          std::array<std::array<double,3>,3> rotation_matrices;
                          rotation_matrices[0][0] = (Vx * Sx - ct);
                          rotation_matrices[0][1] = (Vx * Sy - st);
                          rotation_matrices[0][2] = Vx * Vz;

                          rotation_matrices[1][0] = (Vy * Sx + st);
                          rotation_matrices[1][1] = (Vy * Sy - ct);
                          rotation_matrices[1][2] = Vy * Vz;

                          rotation_matrices[2][0] = Vz * Sx;
                          rotation_matrices[2][1] = Vz * Sy;
                          rotation_matrices[2][2] = 1.0 - z;   // This equals Vz * Vz - 1.0

                          // Rotate the basis rotation matrix with the random uniform distribution rotation matrix
                          // First get the transpose of the rotation matrix
                          std::array<std::array<double, 3>, 3> rot_T= rotation_matrices;
                          rot_T[0][1] = rotation_matrices[1][0];
                          rot_T[1][0] = rotation_matrices[0][1];
                          rot_T[1][2] = rotation_matrices[2][1];
                          rot_T[2][1] = rotation_matrices[1][2];
                          rot_T[0][2] = rotation_matrices[2][0];
                          rot_T[2][0] = rotation_matrices[0][2];

                          // Then U' = R * U * R^T
                          std::array<std::array<double,3>,3> result1 = multiply_3x3_matrices(rotation_matrices, basis_rotation_matrices[i]);

                          it_rotation_matrices = result1;
                        }

                      double total_size = 0;
                      for (auto &&it_sizes : grains_local.sizes)
                        {
                          it_sizes = grain_sizes[i] < 0 ? dist(world->get_random_number_engine()) : grain_sizes[i];
                          total_size += it_sizes;
                        }

                      if (normalize_grain_sizes[i])
                        {
                          const double one_over_total_size = 1/total_size;
                          // std::transform is a c++17 feature, while current GWB is c++14
                          // update this after switching to c+=17
                          // std::transform(grains_local.sizes.begin(), grains_local.sizes.end(), grains_local.sizes.begin(),
                          //                [one_over_total_size](double sizes) -> double { return sizes *one_over_total_size; });
                          // Apply the transformation using a loop
                          for (auto &&size : grains_local.sizes)
                            {
                              size = size*one_over_total_size;
                            }
                        }

                      return grains_local;
                    }
                }
            }
          return grains_local;
        }
        WB_REGISTER_FEATURE_FAULT_GRAINS_MODEL(RandomUniformDistributionDeflected, random uniform distribution deflected)
      } // namespace Grains
    } // namespace FaultModels
  } // namespace Features
} // namespace WorldBuilder

