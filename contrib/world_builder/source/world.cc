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

#include <sstream>

#include <boost/property_tree/json_parser.hpp>

#include <world_builder/world.h>
#include <world_builder/utilities.h>
#include <world_builder/assert.h>
#include <world_builder/point.h>


namespace WorldBuilder
{
  World::World(std::string filename)
    :
    potential_mantle_temperature(1600),
    thermal_expansion_coefficient_alfa(3.5e-5),
    specific_heat_Cp(1250)

  {
    // Get world builder file and check wether it exists
    AssertThrow(access( filename.c_str(), F_OK ) != -1,
                "Could not find the world builder file at the specified location: " + filename);


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
    AssertThrow (value, "Entry undeclared:  Surface rotation angle");
    surface_rotation_angle = Utilities::string_to_double(value.get());

    boost::optional<ptree &> child = tree.get_child("Surface rotation point");
    AssertThrow (child, "Entry undeclared:  Surface rotation point");
    for (boost::property_tree::ptree::const_iterator it = child.get().begin(); it != child.get().end(); ++it)
      {
        surface_rotation_point.push_back(Utilities::string_to_double(it->second.get<std::string>("")));
      }
    AssertThrow (surface_rotation_point.size() == 2, "Only 2d coordinates allowed for Surface rotation point.");

    value  = tree.get_optional<std::string> ("Minimum parts per distance unit");
    AssertThrow (value, "Entry undeclared:  Minimum parts per distance unit");
    minimum_parts_per_distance_unit = Utilities::string_to_unsigned_int(value.get());

    value  = tree.get_optional<std::string> ("Minimum distance points");
    AssertThrow (value, "Entry undeclared:  Minimum distance points");
    minimum_distance_points = Utilities::string_to_double(value.get());

    child = tree.get_child("Surface objects");
    AssertThrow (child, "Entry undeclared:  Surface rotation point");
    for (boost::property_tree::ptree::iterator it = child.get().begin(); it != child.get().end(); ++it)
      {
        features.push_back(Features::create_feature(it->first,this));
        features.back()->read(it->second);
      }


    child = tree.get_child("Coordinate system");
    if (child)
      {
        value = child.get().get_optional<std::string>("name");

        AssertThrow (value, "Entry undeclared:  Coordinate system.name");
        coordinate_system = CoordinateSystems::create_coordinate_system(value.get());
      }
    else
      {
        coordinate_system = CoordinateSystems::create_coordinate_system("cartesian");
      }



  }

  double
  World::temperature(const std::array<double,2> /*point*/, const double /*depth*/, const double /*gravity*/) const
  {
    // turn it into a 3d coordinate and call the 3d temperature function
    return 1.0;
  }

  double
  World::temperature(const std::array<double,3> point_, const double depth, const double gravity_norm) const
  {
    Point<3> point(point_);
    double temperature = potential_mantle_temperature + (((potential_mantle_temperature * thermal_expansion_coefficient_alfa * gravity_norm) / specific_heat_Cp) * 1000.0) * ((depth) / 1000.0);;
    for (std::vector<Features::Interface *>::const_iterator it = features.begin(); it != features.end(); ++it)
      {
        temperature = (*it)->temperature(point,depth,gravity_norm,temperature);
      }

    return temperature;
  }

  double
  World::composition(const std::array<double,2> /*point*/, const double /*depth*/, const unsigned int /*composition_number*/) const
  {
    // turn it into a 3d coordinate and call the 3d temperature function
    return 1.0;
  }

  double
  World::composition(const std::array<double,3> point_, const double depth, const unsigned int composition_number) const
  {
    Point<3> point(point_);
    double composition = 0;
    for (std::vector<Features::Interface *>::const_iterator it = features.begin(); it != features.end(); ++it)
      {
        composition = (*it)->composition(point,depth,composition_number, composition);
      }

    return composition;
  }

  WorldBuilder::CoordinateSystems::Interface *
  World::get_coordinate_system() const
  {
    return coordinate_system;
  }
}

