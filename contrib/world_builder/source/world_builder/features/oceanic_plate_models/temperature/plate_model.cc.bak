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

#include "world_builder/features/oceanic_plate_models/temperature/plate_model.h"

#include "world_builder/nan.h"
#include "world_builder/types/array.h"
#include "world_builder/types/double.h"
#include "world_builder/types/object.h"
#include "world_builder/types/one_of.h"
#include "world_builder/types/point.h"
#include "world_builder/types/value_at_points.h"
#include "world_builder/utilities.h"
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
        PlateModel::PlateModel(WorldBuilder::World *world_)
          :
          min_depth(NaN::DSNAN),
          max_depth(NaN::DSNAN),
          top_temperature(NaN::DSNAN),
          bottom_temperature(NaN::DSNAN),
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
          // Add temperature to the required parameters.
          prm.declare_entry("", Types::Object({"max depth"}),
                            "Plate model.");

          // Declare entries of this plugin
          prm.declare_entry("min depth", Types::OneOf(Types::Double(0),Types::Array(Types::ValueAtPoints(0., 2.))),
                            "The depth in meters from which the temperature of this feature is present.");

          prm.declare_entry("max depth", Types::OneOf(Types::Double(std::numeric_limits<double>::max()),Types::Array(Types::ValueAtPoints(std::numeric_limits<double>::max(), 2.))),
                            "The depth in meters to which the temperature of this feature is present.");

          prm.declare_entry("top temperature", Types::Double(293.15),
                            "The temperature in degree Kelvin which this feature should have");

          prm.declare_entry("bottom temperature", Types::Double(-1),
                            "The temperature in degree Kelvin which this feature should have");

          prm.declare_entry("spreading velocity", Types::OneOf(Types::Double(0.05),Types::Array(Types::ValueAtPoints(0.05, std::numeric_limits<uint64_t>::max()))),
                            "The spreading velocity of the plate in meter per year. "
                            "This is the velocity with which one side moves away from the ridge.");

          prm.declare_entry("ridge coordinates", Types::Array(Types::Array(Types::Point<2>(), 2),1),
                            "An list of ridges. Each ridge is a lists of at least 2 2d points which "
                            "define the location of the ridge. You need to define at least one ridge."
                            "So the an example with two ridges is "
                            "[[[10,20],[20,30],[10,40]],[[50,10],[60,10]]].");

        }

        void
        PlateModel::parse_entries(Parameters &prm, const std::vector<Point<2>> &coordinates)
        {

          min_depth_surface = Objects::Surface(prm.get("min depth",coordinates));
          min_depth = min_depth_surface.minimum;
          max_depth_surface = Objects::Surface(prm.get("max depth",coordinates));
          max_depth = max_depth_surface.maximum;
          operation = string_operations_to_enum(prm.get<std::string>("operation"));
          top_temperature = prm.get<double>("top temperature");
          bottom_temperature = prm.get<double>("bottom temperature");
          spreading_velocities = prm.get_value_at_array("spreading velocity");

          mid_oceanic_ridges = prm.get_vector<std::vector<Point<2>>>("ridge coordinates");
          const double dtr = prm.coordinate_system->natural_coordinate_system() == spherical ? Consts::PI / 180.0 : 1.0;
          for (auto &ridge_coordinates : mid_oceanic_ridges)
            for (auto &ridge_coordinate : ridge_coordinates)
              {
                ridge_coordinate *= dtr;
              }

          unsigned int ridge_point_index = 0;
          for (const auto &mid_oceanic_ridge : mid_oceanic_ridges)
            {
              std::vector<double> spreading_rates_for_ridge;
              for (unsigned int index_y = 0; index_y < mid_oceanic_ridge.size(); index_y++)
                {
                  if (spreading_velocities.second.size() <= 1)
                    {
                      spreading_rates_for_ridge.push_back(spreading_velocities.second[0]);
                    }
                  else
                    {
                      spreading_rates_for_ridge.push_back(spreading_velocities.second[ridge_point_index]);
                    }
                  ridge_point_index += 1;
                }
              spreading_velocities_at_each_ridge_point.push_back(spreading_rates_for_ridge);
            }
        }


        double
        PlateModel::get_temperature(const Point<3> &position,
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
                  Objects::NaturalCoordinate position_in_natural_coordinates_at_min_depth = Objects::NaturalCoordinate(position,
                                                                                            *(world->parameters.coordinate_system));
                  position_in_natural_coordinates_at_min_depth.get_ref_depth_coordinate() += depth-min_depth;
                  std::vector<std::vector<double>> subducting_plate_velocities = {{0}};
                  std::vector<double> ridge_migration_times = {0.0};
                  double bottom_temperature_local = bottom_temperature;

                  if (bottom_temperature_local < 0)
                    {
                      bottom_temperature_local =  this->world->potential_mantle_temperature *
                                                  std::exp(((this->world->thermal_expansion_coefficient* gravity_norm) /
                                                            this->world->specific_heat) * depth);
                    }

                  const int summation_number = 100;

                  std::vector<double> ridge_parameters = Utilities::calculate_ridge_distance_and_spreading(mid_oceanic_ridges,
                                                         spreading_velocities_at_each_ridge_point,
                                                         world->parameters.coordinate_system,
                                                         position_in_natural_coordinates_at_min_depth,
                                                         subducting_plate_velocities,
                                                         ridge_migration_times);

                  const double thermal_diffusivity = this->world->thermal_diffusivity;
                  const double age = ridge_parameters[1] / ridge_parameters[0];
                  double temperature = top_temperature + (bottom_temperature_local - top_temperature) * (depth / max_depth);

                  // This formula addresses the horizontal heat transfer by having the spreading velocity and distance to the ridge in it.
                  // (Chapter 7 Heat, Fowler M. The solid earth: an introduction to global geophysics[M]. Cambridge University Press, 1990)
                  for (int i = 1; i<summation_number+1; ++i)
                    {
                      temperature = temperature + (bottom_temperature_local - top_temperature) *
                                    ((2 / (double(i) * Consts::PI)) * std::sin((double(i) * Consts::PI * depth) / max_depth) *
                                     std::exp((((ridge_parameters[0] * max_depth)/(2 * thermal_diffusivity)) -
                                               std::sqrt(((ridge_parameters[0]*ridge_parameters[0]*max_depth*max_depth) /
                                                          (4*thermal_diffusivity*thermal_diffusivity)) + double(i) * double(i) * Consts::PI * Consts::PI)) *
                                              ((ridge_parameters[0] * age) / max_depth)));

                    }

                  WBAssert(!std::isnan(temperature), "Temperature inside plate model is not a number: " << temperature
                           << ". Relevant variables: bottom_temperature_local = " << bottom_temperature_local
                           << ", top_temperature = " << top_temperature
                           << ", max_depth = " << max_depth
                           << ", spreading_velocity = " << ridge_parameters[0]
                           << ", thermal_diffusivity = " << thermal_diffusivity
                           << ", age = " << age << '.');
                  WBAssert(std::isfinite(temperature), "Temperature inside plate model is not a finite: " << temperature                           << ". Relevant variables: bottom_temperature_local = " << bottom_temperature_local
                           << ", top_temperature = " << top_temperature
                           << ", spreading_velocity = " << ridge_parameters[0]
                           << ", thermal_diffusivity = " << thermal_diffusivity
                           << ", age = " << age << '.');


                  return apply_operation(operation,temperature_,temperature);

                }
            }
          return temperature_;
        }

        WB_REGISTER_FEATURE_OCEANIC_PLATE_TEMPERATURE_MODEL(PlateModel, plate model)
      } // namespace Temperature
    } // namespace OceanicPlateModels
  } // namespace Features
} // namespace WorldBuilder

