/*
  Copyright (C) 2018 - 2021 by the authors of the World Builder code.

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

#include "world_builder/features/oceanic_plate_models/temperature/half_space_model.h"


#include "world_builder/kd_tree.h"
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
        HalfSpaceModel::HalfSpaceModel(WorldBuilder::World *world_)
          :
          min_depth(NaN::DSNAN),
          max_depth(NaN::DSNAN),
          top_temperature(NaN::DSNAN),
          bottom_temperature(NaN::DSNAN),
          spreading_velocity(NaN::DSNAN),
          operation(Operations::REPLACE)
        {
          this->world = world_;
          this->name = "half space model";
        }

        HalfSpaceModel::~HalfSpaceModel()
          = default;

        void
        HalfSpaceModel::declare_entries(Parameters &prm, const std::string & /*unused*/)
        {
          // Document plugin and require entries if needed.
          // Add `ridge coordinates`, `spreading velocity`, `max depth` to the required parameters.
          prm.declare_entry("", Types::Object({"ridge coordinates", "spreading velocity", "max depth"}),
                            "Half space cooling mode");

          // Declare entries of this plugin
          prm.declare_entry("min depth", Types::OneOf(Types::Double(0),Types::Array(Types::ValueAtPoints(0.))),
                            "The depth in meters from which the temperature of this feature is present.");

          prm.declare_entry("max depth", Types::OneOf(Types::Double(std::numeric_limits<double>::max()),Types::Array(Types::ValueAtPoints(std::numeric_limits<double>::max()))),
                            "The depth in meters to which the temperature of this feature is present."
                            "Because half-space reaches background temperature asymptotically, "
                            "this value should be ~2 times the nominal plate thickness of 100 km" );

          prm.declare_entry("top temperature", Types::Double(293.15),
                            "The actual surface temperature in degree Kelvin for this feature.");

          prm.declare_entry("bottom temperature", Types::Double(-1),
                            "The mantle temperature for the half-space cooling model"
                            "in degree Kelvin for this feature. If the model has an adiabatic gradient"
                            "this should be the mantle potential temperature, and T = Tad + Thalf. ");

          prm.declare_entry("spreading velocity", Types::Double(-1),
                            "The spreading velocity of the plate in meter per year. "
                            "This is the velocity with which one side moves away from the ridge.");

          prm.declare_entry("ridge coordinates", Types::Array(Types::Array(Types::Point<2>(), 2),1),
                            "An list of ridges. Each ridge is a lists of at least 2 2d points which "
                            "define the location of the ridge. You need to define at least one ridge."
                            "So the an example with two ridges is "
                            "[[[10,20],[20,30],[10,40]],[[50,10],[60,10]]].");
        }

        void
        HalfSpaceModel::parse_entries(Parameters &prm, const std::vector<Point<2>> &coordinates)
        {

          min_depth_surface = Objects::Surface(prm.get("min depth",coordinates));
          min_depth = min_depth_surface.minimum;
          max_depth_surface = Objects::Surface(prm.get("max depth",coordinates));
          max_depth = max_depth_surface.maximum;
          operation = string_operations_to_enum(prm.get<std::string>("operation"));
          top_temperature = prm.get<double>("top temperature");
          bottom_temperature = prm.get<double>("bottom temperature");
          spreading_velocity = prm.get<double>("spreading velocity")/31557600;  // m/seconds

          mid_oceanic_ridges = prm.get_vector<std::vector<Point<2>>>("ridge coordinates");
          const double dtr = prm.coordinate_system->natural_coordinate_system() == spherical ? Consts::PI / 180.0 : 1.0;
          for (auto &ridge_coordinates : mid_oceanic_ridges)
            for (auto &ridge_coordinate : ridge_coordinates)
              {
                ridge_coordinate *= dtr;
              }
        }


        double
        HalfSpaceModel::get_temperature(const Point<3> &position,
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


                  double bottom_temperature_local = bottom_temperature;

                  if (bottom_temperature_local < 0)
                    {
                      bottom_temperature_local =  this->world->potential_mantle_temperature *
                                                  std::exp(((this->world->thermal_expansion_coefficient* gravity_norm) /
                                                            this->world->specific_heat) * depth);
                    }

                  double distance_ridge = std::numeric_limits<double>::max();

                  const CoordinateSystem coordinate_system = world->parameters.coordinate_system->natural_coordinate_system();


                  // first find if the coordinate is on this side of a ridge
                  unsigned int relevant_ridge = 0;
                  const Point<2> check_point(position_in_natural_coordinates_at_min_depth.get_surface_coordinates(),
                                             position_in_natural_coordinates_at_min_depth.get_coordinate_system());

                  // if there is only one ridge, there is no transform
                  if (mid_oceanic_ridges.size() > 1)
                    {
                      // There are more than one ridge, so there are transform faults
                      // Find the first which is on the same side
                      for (relevant_ridge = 0; relevant_ridge < mid_oceanic_ridges.size()-1; relevant_ridge++)
                        {
                          const Point<2> transform_point_0 = mid_oceanic_ridges[relevant_ridge+1][0];
                          const Point<2> transform_point_1 = mid_oceanic_ridges[relevant_ridge][mid_oceanic_ridges[relevant_ridge].size()-1];
                          const Point<2> reference_point   = mid_oceanic_ridges[relevant_ridge][0];

                          const bool reference_on_side_of_line = (transform_point_1[0] - transform_point_0[0])
                                                                 * (reference_point[1] - transform_point_0[1])
                                                                 - (transform_point_1[1] - transform_point_0[1])
                                                                 * (reference_point[0] - transform_point_0[0])
                                                                 < 0;
                          const bool checkpoint_on_side_of_line = (transform_point_1[0] - transform_point_0[0])
                                                                  * (check_point[1] - transform_point_0[1])
                                                                  - (transform_point_1[1] - transform_point_0[1])
                                                                  * (check_point[0] - transform_point_0[0])
                                                                  < 0;


                          if (reference_on_side_of_line == checkpoint_on_side_of_line)
                            {
                              break;
                            }

                        }
                    }

                  for (unsigned int i_coordinate = 0; i_coordinate < mid_oceanic_ridges[relevant_ridge].size() - 1; i_coordinate++)
                    {
                      const Point<2> segment_point0 = mid_oceanic_ridges[relevant_ridge][i_coordinate];
                      const Point<2> segment_point1 = mid_oceanic_ridges[relevant_ridge][i_coordinate + 1];

                      // based on http://geomalgorithms.com/a02-_lines.html
                      const Point<2> v = segment_point1 - segment_point0;
                      const Point<2> w = check_point - segment_point0;

                      const double c1 = (w[0] * v[0] + w[1] * v[1]);
                      const double c2 = (v[0] * v[0] + v[1] * v[1]);

                      Point<2> Pb(coordinate_system);
                      // This part is needed when we want to consider segments instead of lines
                      // If you want to have infinite lines, use only the else statement.

                      if (c1 <= 0)
                        Pb=segment_point0;
                      else if (c2 <= c1)
                        Pb=segment_point1;
                      else
                        Pb = segment_point0 + (c1 / c2) * v;

                      Point<3> compare_point(coordinate_system);

                      compare_point[0] = coordinate_system == cartesian ? Pb[0] :  position_in_natural_coordinates_at_min_depth.get_depth_coordinate();
                      compare_point[1] = coordinate_system == cartesian ? Pb[1] : Pb[0];
                      compare_point[2] = coordinate_system == cartesian ? position_in_natural_coordinates_at_min_depth.get_depth_coordinate() : Pb[1];

                      distance_ridge = std::min(distance_ridge,
                                                this->world->parameters.coordinate_system->distance_between_points_at_same_depth(Point<3>(position_in_natural_coordinates_at_min_depth.get_coordinates(),
                                                    position_in_natural_coordinates_at_min_depth.get_coordinate_system()),
                                                    compare_point));

                    }

                  const double thermal_diffusivity = this->world->thermal_diffusivity;
                  const double age = distance_ridge / spreading_velocity;

                  double  temperature = bottom_temperature_local;

                  temperature = temperature + (age > 0 ? (top_temperature - bottom_temperature_local)*std::erfc(depth/(2*std::sqrt(thermal_diffusivity*age))) : 0.);

                  WBAssert(!std::isnan(temperature), "Temperature inside half-space cooling model is not a number: " << temperature
                           << ". Relevant variables: bottom_temperature_local = " << bottom_temperature_local
                           << ", top_temperature = " << top_temperature
                           << ", max_depth = " << max_depth
                           << ", spreading_velocity = " << spreading_velocity
                           << ", thermal_diffusivity = " << thermal_diffusivity
                           << ", age = " << age << '.');
                  WBAssert(std::isfinite(temperature), "Temperature inside half-space cooling model is not a finite: " << temperature                           << ". Relevant variables: bottom_temperature_local = " << bottom_temperature_local
                           << ", top_temperature = " << top_temperature
                           << ", spreading_velocity = " << spreading_velocity
                           << ", thermal_diffusivity = " << thermal_diffusivity
                           << ", age = " << age << '.');


                  return apply_operation(operation,temperature_,temperature);

                }
            }
          return temperature_;
        }

        WB_REGISTER_FEATURE_OCEANIC_PLATE_TEMPERATURE_MODEL(HalfSpaceModel, half space model)
      } // namespace Temperature
    } // namespace OceanicPlateModels
  } // namespace Features
} // namespace WorldBuilder

