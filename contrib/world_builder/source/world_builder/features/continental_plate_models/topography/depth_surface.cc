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

#include "world_builder/features/continental_plate_models/topography/depth_surface.h"


#include "world_builder/features/continental_plate_models/topography/interface.h"
#include "world_builder/nan.h"
#include "world_builder/types/array.h"
#include "world_builder/types/double.h"
#include "world_builder/types/object.h"
#include "world_builder/types/one_of.h"
#include "world_builder/types/value_at_points.h"

namespace WorldBuilder
{

  using namespace Utilities;

  namespace Features
  {
    namespace ContinentalPlateModels
    {
      namespace Topography
      {
        DepthSurface::DepthSurface(WorldBuilder::World *world_)
          :
          min_depth(NaN::DSNAN),
          max_depth(NaN::DSNAN),
          operation(Operations::REPLACE)
        {
          this->world = world_;
          this->name = "depth surface";
        }

        DepthSurface::~DepthSurface()
          = default;

        void
        DepthSurface::declare_entries(Parameters &prm, const std::string & /*unused*/)
        {
          // Document plugin and require entries if needed.
          // Add `topography` to the required parameters.
          prm.declare_entry("", Types::Object({"topography"}),
                            "DepthSurface topography model. Set the topography to a constant value.");

          // Declare entries of this plugin
          prm.declare_entry("min depth", Types::OneOf(Types::Double(0),
                                                      Types::Array(Types::ValueAtPoints(0.,2)),
                                                      Types::String("")),
                            "The depth in meters from which the composition of this feature is present.");

          prm.declare_entry("max depth", Types::OneOf(Types::Double(std::numeric_limits<double>::max()),
                                                      Types::Array(Types::ValueAtPoints(std::numeric_limits<double>::max(),2)),
                                                      Types::String("")),
                            "The depth in meters to which the composition of this feature is present.");


          prm.declare_entry("topography", Types::OneOf(Types::Double(0),
                                                       Types::Array(Types::ValueAtPoints(0.,2)),
                                                       Types::String("")),
                            "The topography in meters.");

        }

        void
        DepthSurface::parse_entries(Parameters &prm, const std::vector<Point<2>> &coordinates)
        {

          min_depth_surface = Objects::Surface(prm.get("min depth",coordinates));
          min_depth = min_depth_surface.minimum;
          max_depth_surface = Objects::Surface(prm.get("max depth",coordinates));
          max_depth = max_depth_surface.maximum;
          operation = string_operations_to_enum(prm.get<std::string>("operation"));
          topography_surface = Objects::Surface(prm.get("topography",coordinates));
        }


        double
        DepthSurface::get_topography(const Point<3> & /*position_in_cartesian_coordinates*/,
                                     const Objects::NaturalCoordinate &position_in_natural_coordinates,
                                     double /*topography*/) const
        {
          return -topography_surface.local_value(position_in_natural_coordinates.get_surface_point()).interpolated_value;
        }

        WB_REGISTER_FEATURE_CONTINENTAL_PLATE_TOPOGRAPHY_MODEL(DepthSurface, depth surface)
      } // namespace Topography
    } // namespace ContinentalPlateModels
  } // namespace Features
} // namespace WorldBuilder

