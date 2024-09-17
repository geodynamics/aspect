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

#include "world_builder/features/continental_plate_models/temperature/linear.h"


#include "world_builder/nan.h"
#include "world_builder/types/array.h"
#include "world_builder/types/double.h"
#include "world_builder/types/object.h"
#include "world_builder/types/one_of.h"
#include "world_builder/types/value_at_points.h"
#include "world_builder/world.h"


namespace WorldBuilder
{
  using namespace Utilities;

  namespace Features
  {
    using namespace FeatureUtilities;
    namespace ContinentalPlateModels
    {
      namespace Temperature
      {
        Linear::Linear(WorldBuilder::World *world_)
          :
          min_depth(NaN::DSNAN),
          max_depth(NaN::DSNAN),
          top_temperature(NaN::DSNAN),
          bottom_temperature(NaN::DSNAN),
          operation(Operations::REPLACE)
        {
          this->world = world_;
          this->name = "linear";
        }

        Linear::~Linear()
          = default;

        void
        Linear::declare_entries(Parameters &prm, const std::string & /*unused*/)
        {
          // Document plugin and require entries if needed.
          // Add `max depth` to the required parameters.
          prm.declare_entry("", Types::Object({"max depth"}),
                            "Linear temperature model. Can be set to use an adiabatic temperature at the boundaries.");

          // Declare entries of this plugin
          prm.declare_entry("min depth", Types::OneOf(Types::Double(0),Types::Array(Types::ValueAtPoints(0., 2.))),
                            "The depth in meters from which the composition of this feature is present.");

          prm.declare_entry("max depth", Types::OneOf(Types::Double(std::numeric_limits<double>::max()),Types::Array(Types::ValueAtPoints(std::numeric_limits<double>::max(), 2.))),
                            "The depth in meters to which the composition of this feature is present.");

          prm.declare_entry("top temperature", Types::Double(293.15),
                            "The temperature at the top in degree Kelvin of this feature."
                            "If the value is below zero, the an adiabatic temperature is used.");

          prm.declare_entry("bottom temperature", Types::Double(-1),
                            "The temperature at the top in degree Kelvin of this feature. "
                            "If the value is below zero, an adiabatic temperature is used.");

        }

        void
        Linear::parse_entries(Parameters &prm, const std::vector<Point<2>> &coordinates)
        {
          min_depth_surface = Objects::Surface(prm.get("min depth",coordinates));
          min_depth = min_depth_surface.minimum;
          max_depth_surface = Objects::Surface(prm.get("max depth",coordinates));
          max_depth = max_depth_surface.maximum;
          WBAssert(max_depth >= min_depth, "max depth needs to be larger or equal to min depth.");
          operation = string_operations_to_enum(prm.get<std::string>("operation"));
          top_temperature = prm.get<double>("top temperature");
          bottom_temperature = prm.get<double>("bottom temperature");
        }


        double
        Linear::get_temperature(const Point<3> & /*position_in_cartesian_coordinates*/,
                                const Objects::NaturalCoordinate &position_in_natural_coordinates,
                                const double depth,
                                const double gravity_norm,
                                double temperature_,
                                const double feature_min_depth,
                                const double feature_max_depth) const
        {
          if (depth <= max_depth && depth >= min_depth)
            {
              const double min_depth_local = min_depth_surface.constant_value ? min_depth : min_depth_surface.local_value(position_in_natural_coordinates.get_surface_point()).interpolated_value;
              const double max_depth_local = max_depth_surface.constant_value ? max_depth : max_depth_surface.local_value(position_in_natural_coordinates.get_surface_point()).interpolated_value;
              if (depth <= max_depth_local &&  depth >= min_depth_local)
                {
                  const double min_depth_local_local = std::max(feature_min_depth, min_depth_local);
                  const double max_depth_local_local = std::min(feature_max_depth, max_depth_local);

                  double top_temperature_local = top_temperature;
                  if (top_temperature_local < 0)
                    {
                      top_temperature_local =  this->world->potential_mantle_temperature *
                                               std::exp(((this->world->thermal_expansion_coefficient * gravity_norm) /
                                                         this->world->specific_heat) * min_depth_local_local);
                    }

                  double bottom_temperature_local = bottom_temperature;
                  if (bottom_temperature_local < 0)
                    {
                      bottom_temperature_local =  this->world->potential_mantle_temperature *
                                                  std::exp(((this->world->thermal_expansion_coefficient * gravity_norm) /
                                                            this->world->specific_heat) * max_depth_local_local);
                    }

                  const double new_temperature = top_temperature_local + (max_depth_local_local - min_depth_local_local < 10.0*std::numeric_limits<double>::epsilon() ? 0.0 :
                                                                          (depth - min_depth_local) *
                                                                          ((bottom_temperature_local - top_temperature_local) / (max_depth_local_local - min_depth_local_local)));

                  WBAssert(!std::isnan(new_temperature), "Temperature is not a number: " << new_temperature
                           << ", based on a temperature model with the name " << this->name);
                  WBAssert(std::isfinite(new_temperature), "Temperature is not a finite: " << new_temperature
                           << ", based on a temperature model with the name " << this->name
                           << ", top_temperature_local = " << top_temperature_local << ", depth = " << depth << ", min_depth_local = " << min_depth_local
                           << ", bottom_temperature_local= " << bottom_temperature_local << ", top_temperature_local=" << top_temperature_local
                           << ",max_depth_local_local=" << max_depth_local_local << ", min_depth_local_local =" << min_depth_local_local
                           << ", feature_max_depth = " << feature_max_depth << ", feature_max_depth = " << feature_max_depth);

                  return apply_operation(operation,temperature_,new_temperature);
                }
            }


          return temperature_;
        }

        WB_REGISTER_FEATURE_CONTINENTAL_PLATE_TEMPERATURE_MODEL(Linear, linear)
      } // namespace Temperature
    } // namespace ContinentalPlateModels
  } // namespace Features
} // namespace WorldBuilder

