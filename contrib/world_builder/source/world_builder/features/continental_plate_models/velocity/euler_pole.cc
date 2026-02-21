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

#include "world_builder/features/continental_plate_models/velocity/euler_pole.h"

#include "world_builder/assert.h"
#include "world_builder/consts.h"
#include "world_builder/coordinate_system.h"
#include "world_builder/nan.h"
#include "world_builder/types/array.h"
#include "world_builder/types/double.h"
#include "world_builder/types/object.h"
#include "world_builder/types/one_of.h"
#include "world_builder/types/value_at_points.h"
#include "world_builder/utilities.h"


namespace WorldBuilder
{

  using namespace Utilities;

  namespace Features
  {
    namespace ContinentalPlateModels
    {
      namespace Velocity
      {
        EulerPole::EulerPole(WorldBuilder::World *world_)
          :
          min_depth(NaN::DSNAN),
          max_depth(NaN::DSNAN),
          euler_pole(NaN::DSNAN,NaN::DSNAN,NaN::DSNAN,CoordinateSystem::invalid),
          operation(Operations::REPLACE)
        {
          this->world = world_;
          this->name = "euler pole";
        }

        EulerPole::~EulerPole()
          = default;

        void
        EulerPole::declare_entries(Parameters &prm, const std::string & /*unused*/)
        {
          // Document plugin and require entries if needed.
          // Add `velocity` and to the required parameters.
          prm.declare_entry("", Types::Object({"euler pole"}),
                            "Uniform velocity model. Set the velocity to a constant value.");

          // Declare entries of this plugin
          prm.declare_entry("min depth", Types::OneOf(Types::Double(0),
                                                      Types::Array(Types::ValueAtPoints(0.,2)),
                                                      Types::String("")),
                            "The depth in meters from which the composition of this feature is present.");

          prm.declare_entry("max depth", Types::OneOf(Types::Double(std::numeric_limits<double>::max()),
                                                      Types::Array(Types::ValueAtPoints(std::numeric_limits<double>::max(),2)),
                                                      Types::String("")),
                            "The depth in meters to which the composition of this feature is present.");


          prm.declare_entry("euler pole", Types::Array(Types::Double(0.0),2,2),
                            "The euler pole for the plate (longitude, latitude) in degree.");

          prm.declare_entry("angular velocity", Types::Double(0.0), "The angular velocity of the Euler pole in degree/Myr.");

        }

        void
        EulerPole::parse_entries(Parameters &prm, const std::vector<Point<2>> &coordinates)
        {
          constexpr double year_in_seconds = 60*60*24*365.2425;
          WBAssertThrow(prm.coordinate_system->natural_coordinate_system() == CoordinateSystem::spherical,
                        "The Euler pole velocity model can only be used in a spherical coordinate system.");
          min_depth_surface = Objects::Surface(prm.get("min depth",coordinates));
          min_depth = min_depth_surface.minimum;
          max_depth_surface = Objects::Surface(prm.get("max depth",coordinates));
          max_depth = max_depth_surface.maximum;
          operation = string_operations_to_enum(prm.get<std::string>("operation"));
          std::vector<double> euler_pole_degree = prm.get_vector<double>("euler pole");
          const double angular_speed = (1e-6/year_in_seconds)*sin(prm.get<double>("angular velocity")*Consts::PI/180.);
          const double longitude = euler_pole_degree[0]*Consts::PI/180.;
          const double latitude = euler_pole_degree[1]*Consts::PI/180.;
          euler_pole = Point<3>( angular_speed*cos(latitude)*cos(longitude),angular_speed*cos(latitude)*sin(longitude),angular_speed*sin(latitude) ,CoordinateSystem::cartesian);
        }


        std::array<double,3>
        EulerPole::get_velocity(const Point<3> &position_in_cartesian_coordinates,
                                const Objects::NaturalCoordinate &position_in_natural_coordinates,
                                const double depth,
                                const double  /*gravity*/,
                                std::array<double,3> velocity_,
                                const double /*feature_min_depth*/,
                                const double /*feature_max_depth*/) const
        {

          if (depth <= max_depth && depth >= min_depth)
            {
              const double min_depth_local = min_depth_surface.constant_value ? min_depth : min_depth_surface.local_value(position_in_natural_coordinates.get_surface_point()).interpolated_value;
              const double max_depth_local = max_depth_surface.constant_value ? max_depth : max_depth_surface.local_value(position_in_natural_coordinates.get_surface_point()).interpolated_value;
              if (depth <= max_depth_local &&  depth >= min_depth_local)
                {
                  const double position_norm = position_in_cartesian_coordinates.norm();
                  const Point<3> position_normalized = position_in_cartesian_coordinates/position_norm;
                  const Point<3> euler_velocity = position_norm * Utilities::cross_product(euler_pole, position_normalized);
                  //std::terminate();
                  return {{
                      apply_operation(operation,velocity_[0],euler_velocity[0]),
                      apply_operation(operation,velocity_[1],euler_velocity[1]),
                      apply_operation(operation,velocity_[2],euler_velocity[2])
                    }
                  };
                }
            }

          return velocity_;
        }

        WB_REGISTER_FEATURE_CONTINENTAL_PLATE_VELOCITY_MODEL(EulerPole, euler pole)
      } // namespace Velocity
    } // namespace ContinentalPlateModels
  } // namespace Features
} // namespace WorldBuilder

