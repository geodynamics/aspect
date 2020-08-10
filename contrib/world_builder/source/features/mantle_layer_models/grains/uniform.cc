/*
  Copyright (C) 2018 - 2020 by the authors of the World Builder code.

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

#include <world_builder/utilities.h>
#include <world_builder/assert.h>
#include <world_builder/nan.h>
#include <world_builder/parameters.h>

#include <world_builder/types/array.h>
#include <world_builder/types/double.h>
#include <world_builder/types/string.h>
#include <world_builder/types/object.h>
#include <world_builder/types/unsigned_int.h>
#include <world_builder/types/plugin_system.h>
#include <world_builder/features/mantle_layer_models/grains/uniform.h>


namespace WorldBuilder
{
  using namespace Utilities;

  namespace Features
  {
    namespace MantleLayerModels
    {
      namespace Grains
      {
        Uniform::Uniform(WorldBuilder::World *world_)
          :
          min_depth(NaN::DSNAN),
          max_depth(NaN::DSNAN),
          operation("")
        {
          this->world = world_;
          this->name = "uniform";
        }

        Uniform::~Uniform()
        { }

        void
        Uniform::declare_entries(Parameters &prm, const std::string &)
        {
          // Add compositions, rotation matrices and grain size models to the required parameters.
          prm.declare_entry("", Types::Object({"compositions"}), "Uniform grains model object");

          prm.declare_entry("min depth", Types::Double(0),
                            "The depth in meters from which the composition of this feature is present.");
          prm.declare_entry("max depth", Types::Double(std::numeric_limits<double>::max()),
                            "The depth in meters to which the composition of this feature is present.");

          prm.declare_entry("compositions", Types::Array(Types::UnsignedInt(),0),
                            "A list with the integer labels of the composition which are present there.");

          prm.declare_entry("rotation matrices", Types::Array(Types::Array(Types::Array(Types::Double(0),3,3),3,3),0),
                            "A list with the labels of the grains which are present there for each compositions.");

          prm.declare_entry("Euler angles z-x-z", Types::Array(Types::Array(Types::Double(0),3,3),0),
                            "A list with the z-x-z Euler angles of the grains which are present there for each compositions.");

          prm.declare_entry("orientation operation", Types::String("replace", std::vector<std::string> {"replace"}),
                            "Whether the value should replace any value previously defined at this location (replace) or "
                            "add the value to the previously define value (add, not implemented). Replacing implies that all values not "
                            "explicitly defined are set to zero.");

          prm.declare_entry("grain sizes",
                            Types::Array(Types::Double(-1),0),
                            "A list of the size of all of the grains in each composition. If set to <0, the size will be set so that the total is equal to 1.");


        }

        void
        Uniform::parse_entries(Parameters &prm)
        {
          min_depth = prm.get<double>("min depth");
          max_depth = prm.get<double>("max depth");
          compositions = prm.get_vector<unsigned int>("compositions");

          const bool set_euler_angles = prm.check_entry("Euler angles z-x-z");
          const bool set_rotation_matrices = prm.check_entry("rotation matrices");

          WBAssertThrow(!(set_euler_angles == true && set_rotation_matrices == true),
                        "Only Euler angles or Rotation matrices may be set, but both are set for " << prm.get_full_json_path());


          WBAssertThrow(!(set_euler_angles == false && set_rotation_matrices == false),
                        "Euler angles or Rotation matrices have to be set, but neither are set for " << prm.get_full_json_path());

          if (set_euler_angles)
            {
              std::vector<std::array<double,3> > euler_angles_vector = prm.get_vector<std::array<double,3> >("Euler angles z-x-z");
              rotation_matrices.resize(euler_angles_vector.size());
              for (size_t i = 0; i<euler_angles_vector.size(); ++i)
                {
                  rotation_matrices[i] = WorldBuilder::Utilities::euler_angles_to_rotation_matrix(euler_angles_vector[i][0],euler_angles_vector[i][1],euler_angles_vector[i][2]);
                }

            }
          else
            {
              rotation_matrices = prm.get_vector<std::array<std::array<double,3>,3> >("rotation matrices");
            }
          operation = prm.get<std::string>("orientation operation");
          grain_sizes = prm.get_vector<double>("grain sizes");


          WBAssertThrow(compositions.size() == rotation_matrices.size(),
                        "There are not the same amount of compositions (" << compositions.size()
                        << ") and rotation_matrices (" << rotation_matrices.size() << ").");
          WBAssertThrow(compositions.size() == grain_sizes.size(),
                        "There are not the same amount of compositions (" << compositions.size()
                        << ") and grain_sizes (" << grain_sizes.size() << ").");
        }


        WorldBuilder::grains
        Uniform::get_grains(const Point<3> &,
                            const double depth,
                            const unsigned int composition_number,
                            WorldBuilder::grains grains_,
                            const double ,
                            const double) const
        {
          WorldBuilder::grains  grains_local = grains_;
          if (depth <= max_depth && depth >= min_depth)
            {
              for (unsigned int i =0; i < compositions.size(); ++i)
                {
                  if (compositions[i] == composition_number)
                    {
                      for (auto &&it_rotation_matrices : grains_local.rotation_matrices)
                        {
                          it_rotation_matrices = rotation_matrices[i];
                        }

                      const double size = grain_sizes[i] < 0 ? 1.0/static_cast<double>(grains_local.sizes.size()) :  grain_sizes[i];

                      for (auto &&it_sizes : grains_local.sizes)
                        {
                          it_sizes = size;
                        }


                      return grains_local;
                    }
                }
            }
          return grains_local;
        }
        WB_REGISTER_FEATURE_MANTLE_LAYER_GRAINS_MODEL(Uniform, uniform)
      }
    }
  }
}


