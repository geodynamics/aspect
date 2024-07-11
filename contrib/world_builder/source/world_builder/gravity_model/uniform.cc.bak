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

#include "world_builder/gravity_model/uniform.h"
#include "world_builder/types/double.h"
#include "world_builder/types/object.h"
#include "world_builder/world.h"


namespace WorldBuilder
{
  namespace GravityModel
  {
    Uniform::Uniform(WorldBuilder::World *world_)
    {
      this->world = world_;
    }

    Uniform::~Uniform()
      = default;

    void
    Uniform::declare_entries(Parameters &prm, const std::string & /*unused*/)
    {
      // Add depth method to the required parameters.
      prm.declare_entry("", Types::Object(),
                        "Uniform gravity model. It returns the gravity vector in a Cartesian coordinate system at "
                        "a given position, which has a constant magitude for the whole domain. The vector points "
                        "down in cartesian coordinates and to the center of the sphere in spherical coordinates.");

      prm.declare_entry("magnitude",
                        Types::Double(9.81), "The magnitude of the gravity.");
    }

    void
    Uniform::parse_entries(Parameters &prm)
    {
      gravity_magnitude = prm.get<double>("magnitude");
    }


    Point<3>
    Uniform::gravity_vector(Point<3> point) const
    {
      const CoordinateSystem coordinate_system = world->parameters.coordinate_system->natural_coordinate_system();

      switch (coordinate_system)
        {
          case CoordinateSystem::cartesian:
            return Point<3>(0.,0.,-gravity_magnitude,CoordinateSystem::cartesian);
            break;

          case CoordinateSystem::spherical:
            return (point/point.norm())*-gravity_magnitude;
            break;

          default:
            WBAssertThrow(false, "Invalid coordinate system when using the gravity vector function.");
        }

      return Point<3>(NaN::DSNAN,NaN::DSNAN,NaN::DSNAN,CoordinateSystem::invalid);
    }


    double
    Uniform::gravity_norm(Point<3> /*point*/) const
    {
      return gravity_magnitude;
    }


    /**
     * Register plugin
     */
    WB_REGISTER_GRAVITY_MODEL(Uniform, uniform)
  } // namespace GravityModel
} // namespace WorldBuilder

