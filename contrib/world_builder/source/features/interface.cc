/*
  Copyright (C) 2018 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/

#include <world_builder/features/interface.h>
#include <world_builder/features/continental_plate.h>

#include <deal.II/base/exceptions.h>

#include <boost/algorithm/string.hpp>

using dealii::StandardExceptions::ExcMessage;


  namespace WorldBuilder
  {
    namespace Features
    {
      Interface::Interface()
      {

      }

      Interface *
      create_feature(const std::string name, World *world)
      {
        std::string feature_name = boost::algorithm::to_lower_copy(name);
        boost::algorithm::trim(feature_name);
        if (feature_name == "continental plate")
          return new Features::ContinentalPlate(world);
        else
          AssertThrow(false, ExcMessage("Plugin not implemented."));

        return NULL;
      }
    }
  }

