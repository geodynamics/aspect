/*
  Copyright (C) 2018 by the authors of the World Builder code.

  This file is part of the World Builder.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published
   by the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/


#include <boost/algorithm/string.hpp>

#include <world_builder/coordinate_systems/interface.h>
#include <world_builder/coordinate_systems/cartesian.h>
#include <world_builder/assert.h>


namespace WorldBuilder
{
  namespace CoordinateSystems
  {
    Interface::Interface()
    {

    }

    Interface *
    create_coordinate_system(const std::string name)
    {
      std::string feature_name = boost::algorithm::to_lower_copy(name);
      boost::algorithm::trim(feature_name);
      if (feature_name == "cartesian")
        return new CoordinateSystems::Cartesian();
      else
        AssertThrow(false, "Plugin not implemented.");

      return NULL;
    }
  }
}

