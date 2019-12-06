/*
  Copyright (C) 2018 by the authors of the World Builder code.

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
#include <world_builder/features/fault_models/temperature/uniform.h>


namespace WorldBuilder
{
  using namespace Utilities;

  namespace Features
  {
    namespace FaultModels
    {
      namespace Temperature
      {
        Uniform::Uniform(WorldBuilder::World *world_)
          :
          min_depth(NaN::DSNAN),
          max_depth(NaN::DSNAN),
          temperature(NaN::DSNAN),
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

          // Add temperature to the required parameters.
          prm.declare_entry("", Types::Object({"temperature"}), "Temperature model object");


          prm.declare_entry("min distance fault center", Types::Double(0),
                            "todo The depth in meters from which the composition of this feature is present.");

          prm.declare_entry("max distance fault center", Types::Double(std::numeric_limits<double>::max()),
                            "todo The depth in meters to which the composition of this feature is present.");

          prm.declare_entry("temperature", Types::Double(293.15),
                            "The temperature in degree Kelvin which this feature should have");

          prm.declare_entry("operation", Types::String("replace", std::vector<std::string> {"replace", "add", "substract"}),
                            "Whether the value should replace any value previously defined at this location (replace), "
                            "add the value to the previously define value (add) or substract the value to the previously "
                            "define value (substract).");

        }

        void
        Uniform::parse_entries(Parameters &prm)
        {
          min_depth = prm.get<double>("min distance fault center");
          max_depth = prm.get<double>("max distance fault center");
          operation = prm.get<std::string>("operation");
          temperature = prm.get<double>("temperature");
        }


        double
        Uniform::get_temperature(const Point<3> &,
                                 const double ,
                                 const double ,
                                 double temperature_,
                                 const double ,
                                 const double ,
                                 const std::map<std::string,double> &distance_from_plane) const
        {

          if (std::fabs(distance_from_plane.at("distanceFromPlane")) <= max_depth && std::fabs(distance_from_plane.at("distanceFromPlane")) >= min_depth)
            if (operation == "replace")
              {
                return temperature;
              }
            else if ("add")
              return temperature_ + temperature;
            else if ("substract")
              return temperature_ - temperature;

          return temperature_;
        }

        WB_REGISTER_FEATURE_FAULT_TEMPERATURE_MODEL(Uniform, uniform)
      }
    }
  }
}

