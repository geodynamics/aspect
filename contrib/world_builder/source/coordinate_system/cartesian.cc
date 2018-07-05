/*
  Copyright (C) 2018 by the authors of the World Builder code.

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

#include <world_builder/coordinate_systems/cartesian.h>

namespace WorldBuilder
{
  namespace CoordinateSystems
  {
    Cartesian::Cartesian()
    {}

    Cartesian::~Cartesian()
    {}

    void
    Cartesian::read(const ptree &/*tree*/, std::string & /*path*/)
    {}


    CoordinateSystem
    Cartesian::natural_coordinate_system() const
    {
      return CoordinateSystem::cartesian;
    }


    std::array<double,3>
    Cartesian::cartesian_to_natural_coordinates(const std::array<double,3> &position) const
    {
      return position;
    }


    std::array<double,3>
    Cartesian::natural_to_cartesian_coordinates(const std::array<double,3> &position) const
    {
      return position;
    }
  }
}

