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

#include "world_builder/features/oceanic_plate_models/temperature/plate_model_constant_age.h"

#include "world_builder/nan.h"
#include "world_builder/types/array.h"
#include "world_builder/types/double.h"
#include "world_builder/types/object.h"
#include "world_builder/types/one_of.h"
#include "world_builder/types/point.h"
#include "world_builder/types/value_at_points.h"
#include "world_builder/world.h"


namespace WorldBuilder
{
  using namespace Utilities;

  namespace Features
  {
    namespace OceanicPlateModels
    {
      namespace Temperature
      {
        PlateModelConstantAge::PlateModelConstantAge(WorldBuilder::World *world_)
          :
          min_depth(NaN::DSNAN),
          max_depth(NaN::DSNAN),
          top_temperature(NaN::DSNAN),
          bottom_temperature(NaN::DSNAN),
          operation(Operations::REPLACE)
        {
          this->world = world_;
          this->name = "plate model constant age";
        }

        PlateModelConstantAge::~PlateModelConstantAge()
          = default;

        void
        PlateModelConstantAge::declare_entries(Parameters &prm, const std::string & /*unused*/)
        {
          // Document plugin and require entries if needed.
          // Add temperature to the required parameters.
          prm.declare_entry("", Types::Object({"max depth"}),
                            "Plate model, but with a fixed age.");

          // Declare entries of this plugin
          prm.declare_entry("min depth", Types::OneOf(Types::Double(0),Types::Array(Types::ValueAtPoints(0., 2.))),
                            "The depth in meters from which the temperature of this feature is present.");

          prm.declare_entry("max depth", Types::OneOf(Types::Double(std::numeric_limits<double>::max()),Types::Array(Types::ValueAtPoints(std::numeric_limits<double>::max(), 2.))),
                            "The depth in meters to which the temperature of this feature is present.");

          prm.declare_entry("top temperature", Types::Double(293.15),
                            "The temperature in degree Kelvin which this feature should have");

          prm.declare_entry("bottom temperature", Types::Double(-1),
                            "The temperature in degree Kelvin which this feature should have");

          prm.declare_entry("plate age", Types::Double(80e3),
                            "The age of the plate in year. "
                            "This age is assigned to the whole plate. ");
        }

        void
        PlateModelConstantAge::parse_entries(Parameters &prm, const std::vector<Point<2>> &coordinates)
        {

          min_depth_surface = Objects::Surface(prm.get("min depth",coordinates));
          min_depth = min_depth_surface.minimum;
          max_depth_surface = Objects::Surface(prm.get("max depth",coordinates));
          max_depth = max_depth_surface.maximum;
          operation = string_operations_to_enum(prm.get<std::string>("operation"));
          top_temperature = prm.get<double>("top temperature");
          bottom_temperature = prm.get<double>("bottom temperature");
          plate_age = prm.get<double>("plate age")*31557600;
        }


        double
        PlateModelConstantAge::get_temperature(const Point<3> & /*position*/,
                                               const Objects::NaturalCoordinate &position_in_natural_coordinates,
                                               const double depth,
                                               const double gravity_norm,
                                               double temperature_,
                                               const double /*feature_min_depth*/,
                                               const double /*feature_max_depth*/) const
        {
          if (depth <= max_depth && depth >= min_depth)
            {
              const double min_depth_local = min_depth_surface.constant_value ? min_depth : min_depth_surface.local_value(position_in_natural_coordinates.get_surface_point()).interpolated_value;
              const double max_depth_local = max_depth_surface.constant_value ? max_depth : max_depth_surface.local_value(position_in_natural_coordinates.get_surface_point()).interpolated_value;
              if (depth <= max_depth_local &&  depth >= min_depth_local)
                {
                  double bottom_temperature_local = bottom_temperature;

                  if (bottom_temperature_local < 0)
                    {
                      bottom_temperature_local =  this->world->potential_mantle_temperature *
                                                  std::exp(((this->world->thermal_expansion_coefficient* gravity_norm) /
                                                            this->world->specific_heat) * depth);
                    }
                  const int sommation_number = 100;

                  // some aliases
                  //const double top_temperature = top_temperature;
                  const double thermal_diffusivity = this->world->thermal_diffusivity;
                  double temperature = top_temperature + (bottom_temperature_local - top_temperature) * (depth / max_depth);

                  // This formula ignores the horizontal heat transfer by just having the age of the plate in it.
                  // (Chapter 7 Heat, Fowler M. The solid earth: an introduction to global geophysics[M]. Cambridge University Press, 1990)
                  for (int i = 1; i<sommation_number+1; ++i)
                    {
                      temperature = temperature + (bottom_temperature_local - top_temperature) *
                                    ((2 / (double(i) * Consts::PI)) * std::sin((double(i) * Consts::PI * depth) / max_depth) *
                                     std::exp(-1.0 * i * i * Consts::PI * Consts::PI * thermal_diffusivity * plate_age / (max_depth * max_depth)));
                    }

                  WBAssert(!std::isnan(temperature), "Temperature inside plate model constant age is not a number: " << temperature
                           << ". Relevant variables: bottom_temperature_local = " << bottom_temperature_local
                           << ", top_temperature = " << top_temperature
                           << ", max_depth = " << max_depth
                           << ", thermal_diffusivity = " << thermal_diffusivity
                           << ", age = " << plate_age << '.');
                  WBAssert(std::isfinite(temperature), "Temperature inside plate model constant age is not a finite: " << temperature                           << ". Relevant variables: bottom_temperature_local = " << bottom_temperature_local
                           << ", top_temperature = " << top_temperature
                           << ", thermal_diffusivity = " << thermal_diffusivity
                           << ", age = " << plate_age << '.');


                  return apply_operation(operation,temperature_,temperature);

                }
            }
          return temperature_;
        }

        WB_REGISTER_FEATURE_OCEANIC_PLATE_TEMPERATURE_MODEL(PlateModelConstantAge, plate model constant age)
      } // namespace Temperature
    } // namespace OceanicPlateModels
  } // namespace Features
} // namespace WorldBuilder

