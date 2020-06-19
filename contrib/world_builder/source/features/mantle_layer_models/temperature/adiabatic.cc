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
#include <world_builder/features/mantle_layer_models/temperature/adiabatic.h>


namespace WorldBuilder
{
  using namespace Utilities;

  namespace Features
  {
    namespace MantleLayerModels
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
          operation(Utilities::Operations::REPLACE)
        {
          this->world = world_;
          this->name = "adiabatic";
        }

        Adiabatic::~Adiabatic()
        { }

        void
        Adiabatic::declare_entries(Parameters &prm, const std::string &)
        {
          prm.declare_entry("min depth", Types::Double(0),
                            "The depth in meters from which the temperature of this feature is present.");

          prm.declare_entry("max depth", Types::Double(std::numeric_limits<double>::max()),
                            "The depth in meters to which the temperature of this feature is present.");

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

          min_depth = prm.get<double>("min depth");
          max_depth = prm.get<double>("max depth");
          operation = Utilities::string_operations_to_enum(prm.get<std::string>("operation"));

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
                   "potential_mantle_temperature is not a number: " << potential_mantle_temperature << ".");
          WBAssert(std::isfinite(potential_mantle_temperature),
                   "potential_mantle_temperature is not a finite: " << potential_mantle_temperature << ".");

          WBAssert(!std::isnan(thermal_expansion_coefficient),
                   "thermal_expansion_coefficient is not a number: " << thermal_expansion_coefficient << ".");
          WBAssert(std::isfinite(thermal_expansion_coefficient),
                   "thermal_expansion_coefficient is not a finite: " << thermal_expansion_coefficient << ".");

          WBAssert(!std::isnan(specific_heat),
                   "specific_heat is not a number: " << specific_heat << ".");
          WBAssert(std::isfinite(thermal_expansion_coefficient),
                   "specific_heat is not a finite: " << specific_heat << ".");

        }


        double
        Adiabatic::get_temperature(const Point<3> &,
                                   const double depth,
                                   const double gravity_norm,
                                   double temperature_,
                                   const double ,
                                   const double ) const
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
                       "adabatic_temperature is not a finite: " << adabatic_temperature << ".");

              return Utilities::apply_operation(operation,temperature_,adabatic_temperature);
            }

          return temperature_;
        }

        WB_REGISTER_FEATURE_MANTLE_LAYER_TEMPERATURE_MODEL(Adiabatic, adiabatic)
      }
    }
  }
}

