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

#include "world_builder/features/subducting_plate_models/temperature/linear.h"


#include "world_builder/nan.h"
#include "world_builder/types/double.h"
#include "world_builder/types/object.h"
#include "world_builder/utilities.h"
#include "world_builder/world.h"


namespace WorldBuilder
{
  using namespace Utilities;

  namespace Features
  {
    namespace SubductingPlateModels
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
          // Add `max distance slab top` center to the required parameters.
          prm.declare_entry("", Types::Object({"max distance slab top"}),
                            "Linear temperature model. Can be set to use an adiabatic temperature at the boundaries.");

          // Declare entries of this plugin
          prm.declare_entry("min distance slab top", Types::Double(0),
                            "todo The depth in meters from which the composition of this feature is present.");

          prm.declare_entry("max distance slab top", Types::Double(std::numeric_limits<double>::max()),
                            "todo The depth in meters to which the composition of this feature is present.");

          prm.declare_entry("top temperature", Types::Double(293.15),
                            "The temperature at the top in degree Kelvin of this feature."
                            "If the value is below zero, the an adiabatic temperature is used.");

          prm.declare_entry("bottom temperature", Types::Double(-1),
                            "The temperature at the bottom in degree Kelvin of this feature. "
                            "If the value is below zero, an adiabatic temperature is used.");
        }

        void
        Linear::parse_entries(Parameters &prm)
        {
          min_depth = prm.get<double>("min distance slab top");
          max_depth = prm.get<double>("max distance slab top");
          WBAssert(max_depth >= min_depth, "max depth needs to be larger or equal to min depth.");
          operation = string_operations_to_enum(prm.get<std::string>("operation"));
          top_temperature = prm.get<double>("top temperature");
          bottom_temperature = prm.get<double>("bottom temperature");
        }


        double
        Linear::get_temperature(const Point<3> & /*position_in_cartesian_coordinates*/,
                                const double /*depth*/,
                                const double gravity_norm,
                                double temperature_,
                                const double /*feature_min_depth*/,
                                const double /*feature_max_depth*/,
                                const WorldBuilder::Utilities::PointDistanceFromCurvedPlanes &distance_from_plane,
                                const AdditionalParameters & /*additional_parameters*/) const
        {
          if (distance_from_plane.distance_from_plane <= max_depth && distance_from_plane.distance_from_plane >= min_depth)
            {
              const double min_depth_local = min_depth;
              const double max_depth_local = max_depth;


              double top_temperature_local = top_temperature;
              if (top_temperature_local < 0)
                {
                  top_temperature_local =  this->world->potential_mantle_temperature *
                                           std::exp(((this->world->thermal_expansion_coefficient * gravity_norm) /
                                                     this->world->specific_heat) * min_depth_local);
                }

              double bottom_temperature_local = bottom_temperature;
              if (bottom_temperature_local < 0)
                {
                  bottom_temperature_local =  this->world->potential_mantle_temperature *
                                              std::exp(((this->world->thermal_expansion_coefficient * gravity_norm) /
                                                        this->world->specific_heat) * max_depth_local);
                }

              const double new_temperature = top_temperature_local +
                                             (distance_from_plane.distance_from_plane - min_depth_local) *
                                             ((bottom_temperature_local - top_temperature_local) / (max_depth_local - min_depth_local));
              return apply_operation(operation,temperature_,new_temperature);


            }
          return temperature_;
        }

        WB_REGISTER_FEATURE_SUBDUCTING_PLATE_TEMPERATURE_MODEL(Linear, linear)
      } // namespace Temperature
    } // namespace SubductingPlateModels
  } // namespace Features
} // namespace WorldBuilder

