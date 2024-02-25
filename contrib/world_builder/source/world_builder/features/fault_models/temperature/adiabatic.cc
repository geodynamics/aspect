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

#include "world_builder/features/fault_models/temperature/adiabatic.h"


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
    namespace FaultModels
    {
      namespace Temperature
      {
        Adiabatic::Adiabatic(WorldBuilder::World *world_)
          :
          min_depth(NaN::DSNAN),
          max_depth(NaN::DSNAN),
          potential_mantle_temperature(NaN::DSNAN),
          thermal_expansion_coefficient(NaN::DSNAN),
          specific_heat(NaN::DSNAN),
          operation(Operations::REPLACE)
        {
          this->world = world_;
          this->name = "adiabatic";
        }

        Adiabatic::~Adiabatic()
          = default;

        void
        Adiabatic::declare_entries(Parameters &prm, const std::string & /*unused*/)
        {
          // Document plugin and require entries if needed.
          prm.declare_entry("", Types::Object(),
                            "Adiabatic temperature model. Uses global values by default.");

          // Declare entries of this plugin
          prm.declare_entry("min distance fault center", Types::Double(0),
                            "todo The depth in meters from which the composition of this feature is present.");

          prm.declare_entry("max distance fault center", Types::Double(std::numeric_limits<double>::max()),
                            "todo The depth in meters to which the composition of this feature is present.");

          prm.declare_entry("potential mantle temperature", Types::Double(-1),
                            "The potential temperature of the mantle at the surface in Kelvin. "
                            "If the value is lower then zero, the global value is used.");

          prm.declare_entry("thermal expansion coefficient", Types::Double(-1),
                            "The thermal expansion coefficient in $K^{-1}$. "
                            "If the value is lower then zero, the global value is used.");

          prm.declare_entry("specific heat", Types::Double(-1),
                            "The specific heat in $J kg^{-1} K^{-1}$. "
                            "If the value is lower then zero, the global value is used.");

        }

        void
        Adiabatic::parse_entries(Parameters &prm)
        {

          min_depth = prm.get<double>("min distance fault center");
          max_depth = prm.get<double>("max distance fault center");
          operation = string_operations_to_enum(prm.get<std::string>("operation"));

          potential_mantle_temperature = prm.get<double>("potential mantle temperature");
          if (potential_mantle_temperature < 0)
            potential_mantle_temperature =  this->world->potential_mantle_temperature;


          thermal_expansion_coefficient = prm.get<double>("thermal expansion coefficient");
          if (thermal_expansion_coefficient < 0)
            thermal_expansion_coefficient =  this->world->thermal_expansion_coefficient;

          specific_heat = prm.get<double>("specific heat");
          if (specific_heat < 0)
            specific_heat =  this->world->specific_heat;

          // Some assertions in debug mode can't hurt and have helped before.
          WBAssert(!std::isnan(potential_mantle_temperature),
                   "potential_mantle_temperature is not a number: " << potential_mantle_temperature << '.');
          WBAssert(std::isfinite(potential_mantle_temperature),
                   "potential_mantle_temperature is not a finite: " << potential_mantle_temperature << '.');

          WBAssert(!std::isnan(thermal_expansion_coefficient),
                   "thermal_expansion_coefficient is not a number: " << thermal_expansion_coefficient << '.');
          WBAssert(std::isfinite(thermal_expansion_coefficient),
                   "thermal_expansion_coefficient is not a finite: " << thermal_expansion_coefficient << '.');

          WBAssert(!std::isnan(specific_heat),
                   "specific_heat is not a number: " << specific_heat << '.');
          WBAssert(std::isfinite(thermal_expansion_coefficient),
                   "specific_heat is not a finite: " << specific_heat << '.');

        }


        double
        Adiabatic::get_temperature(const Point<3> & /*position_in_cartesian_coordinates*/,
                                   const double depth,
                                   const double gravity_norm,
                                   double temperature_,
                                   const double  /*feature_min_depth*/,
                                   const double  /*feature_max_depth*/,
                                   const WorldBuilder::Utilities::PointDistanceFromCurvedPlanes & /*distance_from_planes*/,
                                   const AdditionalParameters & /*additional_parameters*/) const
        {

          if (depth <= max_depth && depth >= min_depth)
            {
              const double adabatic_temperature = potential_mantle_temperature *
                                                  std::exp(((thermal_expansion_coefficient * gravity_norm) /
                                                            specific_heat) * depth);


              WBAssert(!std::isnan(adabatic_temperature),
                       "adabatic_temperature is not a number: " << adabatic_temperature << ". "
                       <<"Parameters: potential_mantle_temperature = " << potential_mantle_temperature
                       <<", thermal_expansion_coefficient = " << thermal_expansion_coefficient
                       << ", gravity_norm = " << gravity_norm
                       << ", specific_heat = " << specific_heat
                       << ", depth = " << depth);

              WBAssert(std::isfinite(adabatic_temperature),
                       "adabatic_temperature is not a finite: " << adabatic_temperature << '.');

              return apply_operation(operation,temperature_,adabatic_temperature);
            }

          return temperature_;
        }

        WB_REGISTER_FEATURE_FAULT_TEMPERATURE_MODEL(Adiabatic, adiabatic)
      } // namespace Temperature
    } // namespace FaultModels
  } // namespace Features
} // namespace WorldBuilder

