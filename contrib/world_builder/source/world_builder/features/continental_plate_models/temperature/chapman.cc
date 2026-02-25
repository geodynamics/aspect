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

#include "world_builder/features/continental_plate_models/temperature/chapman.h"


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
        Chapman::Chapman(WorldBuilder::World *world_)
          :
          top_temperature(NaN::DSNAN),
          top_heat_flux(NaN::DSNAN),
          thermal_conductivity(NaN::DSNAN),
          heat_production_per_unit_volume(NaN::DSNAN),
          min_depth(NaN::DSNAN),
          max_depth(NaN::DSNAN),
          operation(Operations::REPLACE)
        {
          this->world = world_;
          this->name = "chapman";
        }

        Chapman::~Chapman()
          = default;

        void
        Chapman::declare_entries(Parameters &prm, const std::string & /*unused*/)
        {
          // Declare entries of this plugin
          prm.declare_entry("", Types::Object(),
                            "Continental geotherm using the steady-state 1-D heat conduction equation from Chapman (1986).");

          prm.declare_entry("min depth", Types::OneOf(Types::Double(0),
                                                      Types::Array(Types::ValueAtPoints(0.,2)),
                                                      Types::String("")),
                            "The depth in meters from which the composition of this feature is present.");

          prm.declare_entry("max depth", Types::OneOf(Types::Double(std::numeric_limits<double>::max()),
                                                      Types::Array(Types::ValueAtPoints(std::numeric_limits<double>::max(),2)),
                                                      Types::String("")),
                            "The depth in meters to which the composition of this feature is present.");

          prm.declare_entry("top temperature", Types::Double(293.15),
                            "The temperature at the top surface in K of this feature."
                            "If the value is below zero, then an adiabatic temperature is used.");

          prm.declare_entry("top heat flux", Types::Double(0.055),
                            "The heat flux at the top surface in W m^(-2) of this feature."
                            "The default value is 0.055.");

          prm.declare_entry("thermal conductivity", Types::Double(2.5),
                            "The thermal conductivity in W m^(-1) K^(-1) of this feature."
                            "The default value is 2.5.");

          prm.declare_entry("heat generation per unit volume", Types::Double(1.e-6),
                            "The heat generation per unit volume in W m^(-3) of this feature."
                            "The default value is 1e-6.");

        }

        void
        Chapman::parse_entries(Parameters &prm, const std::vector<Point<2>> &coordinates)
        {
          min_depth_surface = Objects::Surface(prm.get("min depth",coordinates));
          min_depth = min_depth_surface.minimum;
          max_depth_surface = Objects::Surface(prm.get("max depth",coordinates));
          max_depth = max_depth_surface.maximum;
          WBAssert(max_depth >= min_depth, "max depth needs to be larger or equal to min depth.");
          operation = string_operations_to_enum(prm.get<std::string>("operation"));
          thermal_conductivity = prm.get<double>("thermal conductivity");
          heat_production_per_unit_volume = prm.get<double>("heat generation per unit volume");
          top_heat_flux = prm.get<double>("top heat flux");
          top_temperature = prm.get<double>("top temperature");
          WBAssert(!std::isnan(top_temperature), "Top surface temperature is not a number: " << top_temperature
                   << ", based on a temperature model with the name " << this->name);
        }


        double
        Chapman::get_temperature(const Point<3> & /*position_in_cartesian_coordinates*/,
                                 const Objects::NaturalCoordinate &position_in_natural_coordinates,
                                 const double depth,
                                 const double gravity_norm,
                                 double temperature_,
                                 const double feature_min_depth,
                                 const double /*feature_max_depth*/) const
        {
          if (depth <= max_depth && depth >= min_depth)
            {
              const double min_depth_local = min_depth_surface.constant_value ? min_depth : min_depth_surface.local_value(position_in_natural_coordinates.get_surface_point()).interpolated_value;
              const double max_depth_local = max_depth_surface.constant_value ? max_depth : max_depth_surface.local_value(position_in_natural_coordinates.get_surface_point()).interpolated_value;
              if (depth <= max_depth_local && depth >= min_depth_local)
                {
                  const double min_depth_local_local = std::max(feature_min_depth, min_depth_local);

                  double top_temperature_local = top_temperature;

                  // use adiabatic temperature if top temperature is below zero
                  if (top_temperature_local < 0)
                    {
                      top_temperature_local =  this->world->potential_mantle_temperature *
                                               std::exp(((this->world->thermal_expansion_coefficient * gravity_norm) /
                                                         this->world->specific_heat) * min_depth_local_local);
                    }

                  // calculate the temperature at depth z using steady-state 1-D heat conduction equation:
                  // T(z) = t_top + (q_top / k) * (Δz) - (A / (2 * k)) * (Δz)^2
                  const double new_temperature = top_temperature + (top_heat_flux / thermal_conductivity) * (depth-min_depth_local_local) - heat_production_per_unit_volume / (2. * thermal_conductivity) * (depth-min_depth_local_local)*(depth-min_depth_local_local);

                  WBAssert(!std::isnan(new_temperature), "Temperature is not a number: " << new_temperature
                           << ", based on a temperature model with the name " << this->name);
                  WBAssert(std::isfinite(new_temperature), "Temperature is not a finite: " << new_temperature
                           << ", based on a temperature model with the name " << this->name
                           << ", top_temperature_local = " << top_temperature_local
                           << ", depth = " << depth
                           << ", min_depth_local = " << min_depth_local
                           << ", top_heat_flux = " << top_heat_flux
                           << ", thermal_conductivity=" << thermal_conductivity
                           << ", min_depth_local_local =" << min_depth_local_local);

                  return apply_operation(operation,temperature_,new_temperature);
                }
            }


          return temperature_;
        }

        WB_REGISTER_FEATURE_CONTINENTAL_PLATE_TEMPERATURE_MODEL(Chapman, chapman)
      } // namespace Temperature
    } // namespace ContinentalPlateModels
  } // namespace Features
} // namespace WorldBuilder

