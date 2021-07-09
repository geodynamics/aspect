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

#include <world_builder/utilities.h>
#include <world_builder/assert.h>
#include <world_builder/nan.h>
#include <world_builder/parameters.h>

#include <world_builder/types/array.h>
#include <world_builder/types/double.h>
#include <world_builder/types/point.h>
#include <world_builder/types/string.h>
#include <world_builder/types/object.h>
#include <world_builder/features/oceanic_plate_models/temperature/plate_model.h>


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
          spreading_velocity(NaN::DSNAN),
          operation(Utilities::Operations::REPLACE)
        {
          this->world = world_;
          this->name = "plate model";
        }

        PlateModel::~PlateModel()
          = default;

        void
        PlateModel::declare_entries(Parameters &prm, const std::string & /*unused*/)
        {

          // Add temperature to the required parameters.
          prm.declare_entry("", Types::Object({"ridge coordinates", "spreading velocity", "max depth"}), "Temperature model object");

          prm.declare_entry("min depth", Types::Double(0),
                            "The depth in meters from which the temperature of this feature is present.");

          prm.declare_entry("max depth", Types::Double(std::numeric_limits<double>::max()),
                            "The depth in meters to which the temperature of this feature is present.");

          prm.declare_entry("top temperature", Types::Double(293.15),
                            "The temperature in degree Kelvin which this feature should have");

          prm.declare_entry("bottom temperature", Types::Double(-1),
                            "The temperature in degree Kelvin which this feature should have");

          prm.declare_entry("spreading velocity", Types::Double(-1),
                            "The spreading velocity of the plate in meter per year. "
                            "This is the velocity with which one side moves away from the ridge.");

          prm.declare_entry("ridge coordinates", Types::Array(Types::Point<2>(),2),
                            "A list of 2d points which define the location of the ridge.");

        }

        void
        PlateModel::parse_entries(Parameters &prm)
        {

          min_depth = prm.get<double>("min depth");
          max_depth = prm.get<double>("max depth");
          operation = Utilities::string_operations_to_enum(prm.get<std::string>("operation"));
          top_temperature = prm.get<double>("top temperature");
          bottom_temperature = prm.get<double>("bottom temperature");
          spreading_velocity = prm.get<double>("spreading velocity")/31557600;
          ridge_coordinates = prm.get_vector<Point<2> >("ridge coordinates");

          const double dtr = prm.coordinate_system->natural_coordinate_system() == spherical ? const_pi / 180.0 : 1.0;
          for (auto &ridge_coordinate : ridge_coordinates)
            {
              ridge_coordinate *= dtr;
            }
        }


        double
        PlateModel::get_temperature(const Point<3> &position,
                                    const double depth,
                                    const double gravity_norm,
                                    double temperature_,
                                    const double /*feature_min_depth*/,
                                    const double /*feature_max_depth*/) const
        {
          if (depth <= max_depth && depth >= min_depth)
            {
              WorldBuilder::Utilities::NaturalCoordinate natural_coordinate = WorldBuilder::Utilities::NaturalCoordinate(position,
                                                                              *(world->parameters.coordinate_system));


              double bottom_temperature_local = bottom_temperature;

              if (bottom_temperature_local < 0)
                {
                  bottom_temperature_local =  this->world->potential_mantle_temperature *
                                              std::exp(((this->world->thermal_expansion_coefficient* gravity_norm) /
                                                        this->world->specific_heat) * depth);
                }

              const int sommation_number = 100;
              double distance_ridge = std::numeric_limits<double>::max();

              const CoordinateSystem coordinate_system = world->parameters.coordinate_system->natural_coordinate_system();

              for (unsigned int i_ridge = 0; i_ridge < ridge_coordinates.size()-1; i_ridge++)
                {
                  const Point<2> segment_point0 = ridge_coordinates[i_ridge];
                  const Point<2> segment_point1 = ridge_coordinates[i_ridge+1];

                  const Point<2> check_point(natural_coordinate.get_surface_coordinates(),natural_coordinate.get_coordinate_system());
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

                  compare_point[0] = coordinate_system == cartesian ? Pb[0] :  natural_coordinate.get_depth_coordinate();
                  compare_point[1] = coordinate_system == cartesian ? Pb[1] : Pb[0];
                  compare_point[2] = coordinate_system == cartesian ? natural_coordinate.get_depth_coordinate() : Pb[1];

                  distance_ridge = std::min(distance_ridge,this->world->parameters.coordinate_system->distance_between_points_at_same_depth(Point<3>(natural_coordinate.get_coordinates(),natural_coordinate.get_coordinate_system()),compare_point));

                }

              // some aliases
              //const double top_temperature = top_temperature;
              //const double spreading_velocity = spreading_velocity;
              const double thermal_diffusivity = this->world->thermal_diffusivity;
              const double age = distance_ridge / spreading_velocity;
              double temperature = top_temperature + (bottom_temperature_local - top_temperature) * (depth / max_depth);

              for (int i = 1; i<sommation_number+1; ++i)
                {
                  temperature = temperature + (bottom_temperature_local - top_temperature) *
                                ((2 / (double(i) * const_pi)) * std::sin((double(i) * const_pi * depth) / max_depth) *
                                 std::exp((((spreading_velocity * max_depth)/(2 * thermal_diffusivity)) -
                                           std::sqrt(((spreading_velocity*spreading_velocity*max_depth*max_depth) /
                                                      (4*thermal_diffusivity*thermal_diffusivity)) + double(i) * double(i) * const_pi * const_pi)) *
                                          ((spreading_velocity * age) / max_depth)));

                }

              WBAssert(!std::isnan(temperature), "Temparture inside plate model is not a number: " << temperature
                       << ". Relevant variables: bottom_temperature_local = " << bottom_temperature_local
                       << ", top_temperature = " << top_temperature
                       << ", max_depth = " << max_depth
                       << ", spreading_velocity = " << spreading_velocity
                       << ", thermal_diffusivity = " << thermal_diffusivity
                       << ", age = " << age << ".");
              WBAssert(std::isfinite(temperature), "Temparture inside plate model is not a finite: " << temperature                           << ". Relevant variables: bottom_temperature_local = " << bottom_temperature_local
                       << ", top_temperature = " << top_temperature
                       << ", spreading_velocity = " << spreading_velocity
                       << ", thermal_diffusivity = " << thermal_diffusivity
                       << ", age = " << age << ".");


              return Utilities::apply_operation(operation,temperature_,temperature);

            }
          return temperature_;
        }

        WB_REGISTER_FEATURE_OCEANIC_PLATE_TEMPERATURE_MODEL(PlateModel, plate model)
      }
    }
  }
}
