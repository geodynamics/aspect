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

#include "world_builder/features/continental_plate_models/composition/random.h"


#include "world_builder/nan.h"
#include "world_builder/types/array.h"
#include "world_builder/types/double.h"
#include "world_builder/types/object.h"
#include "world_builder/types/one_of.h"
#include "world_builder/types/unsigned_int.h"
#include "world_builder/types/value_at_points.h"
#include "world_builder/world.h"
#include <random>



namespace WorldBuilder
{

  using namespace Utilities;

  namespace Features
  {
    using namespace FeatureUtilities;
    namespace ContinentalPlateModels
    {
      namespace Composition
      {
        Random::Random(WorldBuilder::World *world_)
          :
          min_depth(NaN::DSNAN),
          max_depth(NaN::DSNAN)
        {
          this->world = world_;
          this->name = "random";
        }

        Random::~Random()
          = default;

        void
        Random::declare_entries(Parameters &prm, const std::string & /*unused*/)
        {
          // Document plugin and require entries if needed.
          // Add compositions to the required parameters.
          prm.declare_entry("", Types::Object({"compositions"}),
                            "Uniform compositional model. Sets constant compositional field.");

          // Declare entries of this plugin
          prm.declare_entry("min depth", Types::OneOf(Types::Double(0),Types::Array(Types::ValueAtPoints(0., 2.))),
                            "The depth in meters from which the composition of this feature is present.");

          prm.declare_entry("max depth", Types::OneOf(Types::Double(std::numeric_limits<double>::max()),Types::Array(Types::ValueAtPoints(std::numeric_limits<double>::max(), 2.))),
                            "The depth in meters to which the composition of this feature is present.");

          prm.declare_entry("compositions", Types::Array(Types::UnsignedInt(),0),
                            "A list with the labels of the composition which are present there.");

          prm.declare_entry("min value", Types::Array(Types::Double(0.0),1),
                            "Minimum value of the range within which we want to generate a random compositional value "
                            "corresponding to the compositional field.");

          prm.declare_entry("max value", Types::Array(Types::Double(1.0),1),
                            "Maximum value of the range within which we want to generate a random compositional value "
                            "corresponding to the compositional field.");

          prm.declare_entry("operation", Types::String("replace", std::vector<std::string> {"replace", "replace defined only", "add", "subtract"}),
                            "Whether the value should replace any value previously defined at this location (replace) or "
                            "add the value to the previously define value. Replacing implies that all compositions not "
                            "explicitly defined are set to zero. To only replace the defined compositions use the replace only defined option.");

        }

        void
        Random::parse_entries(Parameters &prm, const std::vector<Point<2>> &coordinates)
        {
          min_depth_surface = Objects::Surface(prm.get("min depth",coordinates));
          min_depth = min_depth_surface.minimum;
          max_depth_surface = Objects::Surface(prm.get("max depth",coordinates));
          max_depth = max_depth_surface.maximum;

          compositions = prm.get_vector<unsigned int>("compositions");
          min_value    = prm.get_vector<double>("min value");
          max_value    = prm.get_vector<double>("max value");
          operation    = string_operations_to_enum(prm.get<std::string>("operation"));
        }


        double
        Random::get_composition(const Point<3> & /*position_in_cartesian_coordinates*/,
                                const Objects::NaturalCoordinate &position_in_natural_coordinates,
                                const double depth,
                                const unsigned int composition_number,
                                double composition,
                                const double  /*feature_min_depth*/,
                                const double  /*feature_max_depth*/) const
        {
          if (depth <= max_depth && depth >= min_depth)
            {
              const double min_depth_local = min_depth_surface.constant_value ? min_depth : min_depth_surface.local_value(position_in_natural_coordinates.get_surface_point()).interpolated_value;
              const double max_depth_local = max_depth_surface.constant_value ? max_depth : max_depth_surface.local_value(position_in_natural_coordinates.get_surface_point()).interpolated_value;
              if (depth <= max_depth_local &&  depth >= min_depth_local)
                {
                  for (unsigned int i : compositions)
                    {
                      if (i == composition_number)
                        {
                          std::uniform_real_distribution<> dist(min_value[0], max_value[0]);
                          const double composition_ = dist(world->get_random_number_engine());
                          return apply_operation(operation,composition,composition_);
                        }
                    }

                  if (operation == Operations::REPLACE)
                    return 0.0;
                }
            }
          return composition;
        }
        WB_REGISTER_FEATURE_CONTINENTAL_PLATE_COMPOSITION_MODEL(Random, random)
      } // namespace Composition
    } // namespace ContinentalPlateModels
  } // namespace Features
} // namespace WorldBuilder


