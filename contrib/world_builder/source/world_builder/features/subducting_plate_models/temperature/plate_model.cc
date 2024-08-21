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

#include "world_builder/features/subducting_plate_models/temperature/plate_model.h"

#include "world_builder/nan.h"
#include "world_builder/types/bool.h"
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
        PlateModel::PlateModel(WorldBuilder::World *world_)
          :
          min_depth(NaN::DSNAN),
          max_depth(NaN::DSNAN),
          density(NaN::DSNAN),
          plate_velocity(NaN::DSNAN),
          thermal_conductivity(NaN::DSNAN),
          thermal_expansion_coefficient(NaN::DSNAN),
          specific_heat(NaN::DSNAN),
          potential_mantle_temperature(NaN::DSNAN),
          surface_temperature(NaN::DSNAN),
          adiabatic_heating(true),
          operation(Operations::REPLACE)
        {
          this->world = world_;
          this->name = "plate model";
        }

        PlateModel::~PlateModel()
          = default;

        void
        PlateModel::declare_entries(Parameters &prm, const std::string & /*unused*/)
        {
          // Document plugin and require entries if needed.
          // Add `plate velocity` to the required parameters.
          prm.declare_entry("", Types::Object({"plate velocity"}),
                            "Plate model (based on McKenzie, 1970).");

          // Declare entries of this plugin
          prm.declare_entry("min distance slab top", Types::Double(0),
                            "todo The depth in meters from which the composition of this feature is present.");

          prm.declare_entry("max distance slab top", Types::Double(std::numeric_limits<double>::max()),
                            "todo The depth in meters to which the composition of this feature is present.");

          prm.declare_entry("density", Types::Double(3300),
                            "The reference density of the subducting plate in $kg/m^3$");

          prm.declare_entry("plate velocity", Types::Double(NaN::DQNAN),
                            "The velocity in meters per year with which the plate subducts in meters per year.");

          prm.declare_entry("thermal conductivity", Types::Double(2.0),
                            "The thermal conductivity of the subducting plate material in $W m^{-1} K^{-1}$.");

          prm.declare_entry("thermal expansion coefficient", Types::Double(-1),
                            "The thermal expansivity of the subducting plate material in $K^{-1}$. If smaller than zero, the global value is used.");

          prm.declare_entry("specific heat", Types::Double(-1),
                            "The specific heat of the subducting plate material in $J kg^{-1} K^{-1}$. If smaller than zero, the global value is used.");

          prm.declare_entry("adiabatic heating", Types::Bool(true),
                            "Whether adiabatic heating should be used for the slab. Setting the parameter to false leads to equation 26 from McKenzie (1970),"
                            "which is the result obtained from McKenzie 1969.");

          prm.declare_entry("potential mantle temperature", Types::Double(-1),
                            "The potential temperature of the mantle at the surface in Kelvin. If smaller than zero, the global value is used.");
        }

        void
        PlateModel::parse_entries(Parameters &prm)
        {

          min_depth = prm.get<double>("min distance slab top");
          max_depth = prm.get<double>("max distance slab top");
          operation = string_operations_to_enum(prm.get<std::string>("operation"));

          density = prm.get<double>("density");
          plate_velocity = prm.get<double>("plate velocity");
          thermal_conductivity = prm.get<double>("thermal conductivity");
          thermal_expansion_coefficient = prm.get<double>("thermal expansion coefficient");

          if (thermal_expansion_coefficient < 0 )
            thermal_expansion_coefficient = this->world->thermal_expansion_coefficient;

          specific_heat = prm.get<double>("specific heat");

          if (specific_heat < 0)
            specific_heat = this->world->specific_heat;

          adiabatic_heating = prm.get<bool>("adiabatic heating");

          potential_mantle_temperature = this->world->potential_mantle_temperature >= 0
                                         ?
                                         this->world->potential_mantle_temperature
                                         :
                                         prm.get<double>("potential mantle temperature");
          surface_temperature = this->world->surface_temperature;
        }


        double
        PlateModel::get_temperature(const Point<3> & /*position_in_cartesian_coordinates*/,
                                    const double depth,
                                    const double gravity_norm,
                                    double temperature_,
                                    const double /*feature_min_depth*/,
                                    const double /*feature_max_depth*/,
                                    const WorldBuilder::Utilities::PointDistanceFromCurvedPlanes &distance_from_planes,
                                    const AdditionalParameters &additional_parameters) const
        {
          const double thickness_local = std::min(additional_parameters.local_thickness, max_depth);
          const double distance_from_plane = distance_from_planes.distance_from_plane;
          const double distance_along_plane = distance_from_planes.distance_along_plane;

          if (distance_from_plane <= max_depth && distance_from_plane >= min_depth)
            {
              /*
               * We now use the McKenzie (1970) equation to determine the
               * temperature inside the slab. The McKenzie equation was
               * designed for a straight slab, but we have a potentially
               * curved slab. Since the angle is used to compute the depth
               * of the point, we directly use the depth.
               */
              const double R = (density * specific_heat
                                * (plate_velocity /(365.25 * 24.0 * 60.0 * 60.0))
                                * thickness_local) / (2.0 * thermal_conductivity);

              WBAssert(!std::isnan(R), "Internal error: R is not a number: " << R << '.');

              const int n_sum = 500;
              // distance_from_plane can be zero, so protect division.
              const double z_scaled = 1 - (std::fabs(distance_from_plane) < 2.0 * std::numeric_limits<double>::epsilon() ?
                                           2.0 * std::numeric_limits<double>::epsilon()
                                           :
                                           distance_from_plane/thickness_local);

              // distance_along_plane can be zero, so protect division.
              const double x_scaled = (std::fabs(distance_along_plane) < 2.0 * std::numeric_limits<double>::epsilon() ?
                                       2.0 *std::numeric_limits<double>::epsilon()
                                       :
                                       distance_along_plane/thickness_local);

              // the paper uses `(x_scaled * sin(average_angle) - z_scaled * cos(average_angle))` to compute the
              // depth (except that you do not use average angles since they only have on angle). On recomputing
              // their result it seems to me (Menno) that it should have been `(1-z_scaled)` instead of `z_scaled`.
              // To avoid this whole problem we just use the depth directly since we have that.
              // todo: get the local thickniss out of H, that prevents an other division.
              // If we want to specify the bottom temperature, because we have defined a linear temperature increase in the
              // mantle and/or oceanic plate, we have to switch off adiabatic heating for now.
              // Todo: there may be a better way to deal with this.
              ;
              const double temp = adiabatic_heating ? std::exp(((thermal_expansion_coefficient * gravity_norm * depth) / specific_heat)) : 1;

              WBAssert(!std::isnan(temp), "Internal error: temp is not a number: " << temp << ". In exponent: "
                       << std::exp(((thermal_expansion_coefficient * gravity_norm) / specific_heat) * depth)
                       << ", thermal_expansion_coefficient = " << thermal_expansion_coefficient << ", gravity_norm = " << gravity_norm
                       << ", specific_heat = "<< specific_heat << ", depth = " << depth );

              double sum=0;
              for (int i=1; i<=n_sum; i++)
                {
                  sum += (std::pow((-1.0),i)/(i*Consts::PI)) *
                         (exp((R - std::pow(R * R + i * i * Consts::PI * Consts::PI, 0.5)) * x_scaled))
                         * (sin(i * Consts::PI * z_scaled));
                }
              // todo: investigate whether this 273.15 should just be the surface temperature.
              const double temperature = temp * (potential_mantle_temperature
                                                 + 2.0 * (potential_mantle_temperature - 273.15) * sum);

              WBAssert(!std::isnan(temperature), "Internal error: temperature is not a number: " << temperature << '.');
              WBAssert(std::isfinite(temperature), "Internal error: temperature is not finite: " << temperature << '.');


              return apply_operation(operation,temperature_,temperature);
            }

          return temperature_;
        }

        WB_REGISTER_FEATURE_SUBDUCTING_PLATE_TEMPERATURE_MODEL(PlateModel, plate model)
      } // namespace Temperature
    } // namespace SubductingPlateModels
  } // namespace Features
} // namespace WorldBuilder

