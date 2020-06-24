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

#include <world_builder/types/double.h>
#include <world_builder/types/string.h>
#include <world_builder/types/object.h>
#include <world_builder/features/subducting_plate_models/temperature/linear.h>


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
          operation(Utilities::Operations::REPLACE)
        {
          this->world = world_;
          this->name = "linear";
        }

        Linear::~Linear()
        { }

        void
        Linear::declare_entries(Parameters &prm, const std::string &)
        {
          // Add max depth to the required parameters.
          prm.declare_entry("", Types::Object({"max distance slab top"}), "Temperature model object");

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
          operation = Utilities::string_operations_to_enum(prm.get<std::string>("operation"));
          top_temperature = prm.get<double>("top temperature");
          bottom_temperature = prm.get<double>("bottom temperature");
        }


        double
        Linear::get_temperature(const Point<3> &,
                                const double,
                                const double gravity_norm,
                                double temperature_,
                                const double,
                                const double,
                                const std::map<std::string,double> &distance_from_plane) const
        {
          if (distance_from_plane.at("distanceFromPlane") <= max_depth && distance_from_plane.at("distanceFromPlane") >= min_depth)
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

              const double new_temperature = top_temperature +
                                             (distance_from_plane.at("distanceFromPlane") - min_depth_local) *
                                             ((bottom_temperature_local - top_temperature_local) / (max_depth_local - min_depth_local));
              return Utilities::apply_operation(operation,temperature_,new_temperature);


            }
          return temperature_;
        }

        WB_REGISTER_FEATURE_SUBDUCTING_PLATE_TEMPERATURE_MODEL(Linear, linear)
      }
    }
  }
}

