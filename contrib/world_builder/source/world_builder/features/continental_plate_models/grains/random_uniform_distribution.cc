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

#include "world_builder/features/continental_plate_models/grains/random_uniform_distribution.h"

#include <algorithm>

#include "world_builder/nan.h"
#include "world_builder/types/array.h"
#include "world_builder/types/bool.h"
#include "world_builder/types/double.h"
#include "world_builder/types/object.h"
#include "world_builder/types/one_of.h"
#include "world_builder/types/unsigned_int.h"
#include "world_builder/types/value_at_points.h"
#include "world_builder/world.h"

namespace WorldBuilder
{
  using namespace Utilities;

  namespace Features
  {
    namespace ContinentalPlateModels
    {
      namespace Grains
      {
        RandomUniformDistribution::RandomUniformDistribution(WorldBuilder::World *world_)
          :
          min_depth(NaN::DSNAN),
          max_depth(NaN::DSNAN)

        {
          this->world = world_;
          this->name = "random uniform distribution";
        }

        RandomUniformDistribution::~RandomUniformDistribution()
          = default;

        void
        RandomUniformDistribution::declare_entries(Parameters &prm, const std::string & /*unused*/)
        {
          // Document plugin and require entries if needed.
          // Add compositions the required parameters.
          prm.declare_entry("", Types::Object({"compositions"}),
                            "Random uniform distribution grains model. The size of the grains can be independently set "
                            "to a single value or to a random distribution.");

          // Declare entries of this plugin
          prm.declare_entry("min depth", Types::OneOf(Types::Double(0),Types::Array(Types::ValueAtPoints(0., 2.))),
                            "The depth in meters from which the composition of this feature is present.");
          prm.declare_entry("max depth", Types::OneOf(Types::Double(std::numeric_limits<double>::max()),Types::Array(Types::ValueAtPoints(std::numeric_limits<double>::max(), 2.))),
                            "The depth in meters to which the composition of this feature is present.");

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



        }

        void
        RandomUniformDistribution::parse_entries(Parameters &prm, const std::vector<Point<2>> &coordinates)
        {
          min_depth_surface = Objects::Surface(prm.get("min depth",coordinates));
          min_depth = min_depth_surface.minimum;
          max_depth_surface = Objects::Surface(prm.get("max depth",coordinates));
          max_depth = max_depth_surface.maximum;
          compositions = prm.get_vector<unsigned int>("compositions");

          operation = prm.get<std::string>("orientation operation");
          grain_sizes = prm.get_vector<double>("grain sizes");
          normalize_grain_sizes = prm.get_vector<bool>("normalize grain sizes");


          WBAssertThrow(compositions.size() == grain_sizes.size(),
                        "There are not the same amount of compositions (" << compositions.size()
                        << ") and grain_sizes (" << grain_sizes.size() << ").");
          WBAssertThrow(compositions.size() == normalize_grain_sizes.size(),
                        "There are not the same amount of compositions (" << compositions.size()
                        << ") and normalize_grain_sizes (" << normalize_grain_sizes.size() << ").");
        }


        WorldBuilder::grains
        RandomUniformDistribution::get_grains(const Point<3> & /*position_in_cartesian_coordinates*/,
                                              const Objects::NaturalCoordinate &position_in_natural_coordinates,
                                              const double depth,
                                              const unsigned int composition_number,
                                              WorldBuilder::grains grains_,
                                              const double  /*feature_min_depth*/,
                                              const double /*feature_max_depth*/) const
        {
          WorldBuilder::grains  grains_local = grains_;
          if (depth <= max_depth && depth >= min_depth)
            {
              const double min_depth_local = min_depth_surface.constant_value ? min_depth : min_depth_surface.local_value(position_in_natural_coordinates.get_surface_point()).interpolated_value;
              const double max_depth_local = max_depth_surface.constant_value ? max_depth : max_depth_surface.local_value(position_in_natural_coordinates.get_surface_point()).interpolated_value;
              if (depth <= max_depth_local &&  depth >= min_depth_local)
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

                              const double theta = 2.0 * Consts::PI * one; // Rotation about the pole (Z)
                              const double phi = 2.0 * Consts::PI * two; // For direction of pole deflection.
                              const double z = 2.0* three; //For magnitude of pole deflection.

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

                              it_rotation_matrices[0][0] = Vx * Sx - ct;
                              it_rotation_matrices[0][1] = Vx * Sy - st;
                              it_rotation_matrices[0][2] = Vx * Vz;

                              it_rotation_matrices[1][0] = Vy * Sx + st;
                              it_rotation_matrices[1][1] = Vy * Sy - ct;
                              it_rotation_matrices[1][2] = Vy * Vz;

                              it_rotation_matrices[2][0] = Vz * Sx;
                              it_rotation_matrices[2][1] = Vz * Sy;
                              it_rotation_matrices[2][2] = 1.0 - z;   // This equals Vz * Vz - 1.0
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
                              std::transform(grains_local.sizes.begin(), grains_local.sizes.end(), grains_local.sizes.begin(),
                                             [one_over_total_size](double sizes) -> double { return sizes *one_over_total_size; });
                            }


                          return grains_local;
                        }
                    }
                }
            }
          return grains_local;
        }
        WB_REGISTER_FEATURE_CONTINENTAL_PLATE_GRAINS_MODEL(RandomUniformDistribution, random uniform distribution)
      } // namespace Grains
    } // namespace ContinentalPlateModels
  } // namespace Features
} // namespace WorldBuilder


