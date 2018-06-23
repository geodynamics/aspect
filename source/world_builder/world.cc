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

#include <aspect/world_builder/world.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/utilities.h>
#include <sstream>
#include <boost/property_tree/json_parser.hpp>

using dealii::StandardExceptions::ExcMessage;

namespace aspect
{
  namespace WorldBuilder
  {
    World::World(std::string filename)
    {
      // Get world builder file and check wether it exists
      AssertThrow(access( filename.c_str(), F_OK ) != -1,
                  ExcMessage("Could not find the world builder file at the specified location: " + filename));


      // Now read in the world builder file into a file stream and
      // put it into a boost property tree.
      std::ifstream json_input_stream(filename.c_str());
      ptree property_tree;
      boost::property_tree::json_parser::read_json (json_input_stream, property_tree);
      this->read(property_tree);
    }

    World::~World()
    {}

    void
    World::read(ptree &tree)
    {

      //todo: wrap this into a convient function

      boost::optional<std::string> value  = tree.get_optional<std::string> ("Surface rotation angle");
      AssertThrow (value, ExcMessage("Entry undeclared:  Surface rotation angle"));
      surface_rotation_angle = dealii::Utilities::string_to_double(value.get());

      boost::optional<ptree &> child = tree.get_child("Surface rotation point");
      AssertThrow (child, ExcMessage("Entry undeclared:  Surface rotation point"));
      for (boost::property_tree::ptree::const_iterator it = child.get().begin(); it != child.get().end(); ++it)
        {
          surface_rotation_point.push_back(dealii::Utilities::string_to_double(it->second.get<std::string>("")));
        }
      AssertThrow (surface_rotation_point.size() == 2, ExcMessage("Only 2d coordinates allowed for Surface rotation point."));

      value  = tree.get_optional<std::string> ("Minimum parts per distance unit");
      AssertThrow (value, ExcMessage("Entry undeclared:  Minimum parts per distance unit"));
      minimum_parts_per_distance_unit = dealii::Utilities::string_to_int(value.get());

      value  = tree.get_optional<std::string> ("Minimum distance points");
      AssertThrow (value, ExcMessage("Entry undeclared:  Minimum distance points"));
      minimum_distance_points = dealii::Utilities::string_to_double(value.get());

      child = tree.get_child("Surface objects");
      AssertThrow (child, ExcMessage("Entry undeclared:  Surface rotation point"));
      for (boost::property_tree::ptree::iterator it = child.get().begin(); it != child.get().end(); ++it)
        {
          features.push_back(Features::create_feature(it->first));
          features.back()->read(it->second);
        }



    }

    double
    World::temperature(const std::array<double,2> point) const
    {
      // turn it into a 3d coordinate and call the 3d temperature function
      return 1;
    }

    double
    World::temperature(const std::array<double,3> point) const
    {
      double temperature = 0;
      for (std::vector<Features::Interface *>::const_iterator it = features.begin(); it != features.end(); ++it)
        {
          (*it)->temperature(point,temperature);
        }
      return 1;
    }
  }
}
