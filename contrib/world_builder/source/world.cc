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

#include <sstream>

#include <boost/property_tree/json_parser.hpp>

#include <world_builder/world.h>
#include <world_builder/utilities.h>
#include <world_builder/assert.h>
#include <world_builder/point.h>
#include <world_builder/nan.h>


namespace WorldBuilder
{

  using namespace Utilities;

  World::World(std::string filename)
    :
    dim(3),
    surface_rotation_angle(0.0),
    minimum_parts_per_distance_unit(NaN::ISNAN),
    minimum_distance_points(NaN::ISNAN),
    potential_mantle_temperature(1600),
    thermal_expansion_coefficient_alpha(3.5e-5),
    specific_heat_Cp(1250),
    coordinate_system(NULL)
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

    //todo: wrap this into a convenient function
    boost::optional<ptree &> child;

    /**
     * get the cross section
     */
    child = tree.get_child_optional("Cross section");
    if (child)
      {
        dim = 2;
        for (boost::property_tree::ptree::const_iterator it = child.get().begin(); it != child.get().end(); ++it)
          {
            std::vector<double> tmp;
            boost::optional<const ptree &> child2 = it->second.get_child("");
            AssertThrow (child, "Cross section: This should be a 2d array, but only one dimension found.");
            for (boost::property_tree::ptree::const_iterator it2 = child2.get().begin(); it2 != child2.get().end(); ++it2)
              {
                tmp.push_back(stod(it2->second.get<std::string>("")));
              }
            AssertThrow (tmp.size() == 2, "Cross section: These represent 2d coordinates, but there are " <<
                         tmp.size() <<
                         " coordinates specified.");

            std::array<double,2> tmp_array;
            std::copy(tmp.begin(), tmp.end(), tmp_array.begin());
            cross_section.push_back(tmp_array);
          }
        AssertThrow (cross_section.size() == 2, "Cross section: This feature requires at exactly 2 coordinates, but " <<
                     cross_section.size() <<
                     " where provided.");

        /**
         * pre-compute stuff for the cross section
         */
        const Point<2> diff_points = Point<2>(cross_section[0])-Point<2>(cross_section[1]);
        const double one_over_cross_section_length = 1/(diff_points.norm());
        surface_coord_conversions = {diff_points[0] *one_over_cross_section_length,diff_points[1] *one_over_cross_section_length};
      }
    else
      {
        dim = 3;
      }


    // Get the rotation angle
    surface_rotation_angle = string_to_double(get_from_ptree(tree, "", "Surface rotation angle"));

    /**
     * Get the point to rotate about when the rotation angle is not zero.
     */
    if (std::abs(surface_rotation_angle) < std::numeric_limits<double>::epsilon())
      {
        child = tree.get_child("Surface rotation point");
        AssertThrow (child, "Entry undeclared:  Surface rotation point. "
                     "Need a surface rotation point when rotation the surface.");
        for (boost::property_tree::ptree::const_iterator it = child.get().begin(); it != child.get().end(); ++it)
          {
            surface_rotation_point.push_back(Utilities::string_to_double(it->second.get<std::string>("")));
          }
        AssertThrow (surface_rotation_point.size() == 2, "Only 2d coordinates allowed for Surface rotation point.");
      }

    // Get variables needed for the interpolation (todo: not yet used).
    minimum_parts_per_distance_unit = string_to_unsigned_int(get_from_ptree(tree, "", "Minimum parts per distance unit"));
    minimum_distance_points = string_to_double(get_from_ptree(tree, "", "Minimum distance points"));

    // Now reading the main paramters from the file
    child = tree.get_child("Surface objects");
    AssertThrow (child, "Entry undeclared:  Surface rotation point");
    for (boost::property_tree::ptree::const_iterator it = child.get().begin(); it != child.get().end(); ++it)
      {
        features.push_back(Features::create_feature(it->first,*this));
        std::string path = "Surface objects" + World::path_seperator + it->first;
        features.back()->read(it->second,path);
      }


    child = tree.get_child("Coordinate system");
    if (child)
      {
        const boost::optional<std::string> value = child.get().get_optional<std::string>("name");

        AssertThrow (value, "Entry undeclared:  Coordinate system.name");
        coordinate_system = CoordinateSystems::create_coordinate_system(value.get());
      }
    else
      {
        coordinate_system = CoordinateSystems::create_coordinate_system("cartesian");
      }



  }

  double
  World::temperature(const std::array<double,2> &point, const double depth, const double gravity_norm) const
  {
    // turn it into a 3d coordinate and call the 3d temperature function
    AssertThrow(dim == 2, "This function can only be called when the cross section "
                "variable in the world builder file has been set.");

    Point<3> coord_3d(cross_section[0][0] + point[0] * surface_coord_conversions[0],
                      cross_section[0][1] + point[0] * surface_coord_conversions[1],
                      point[1]);

    return temperature(coord_3d.get_array(), depth, gravity_norm);
  }

  double
  World::temperature(const std::array<double,3> &point_, const double depth, const double gravity_norm) const
  {
    Point<3> point(point_);
    double temperature = potential_mantle_temperature + (((potential_mantle_temperature * thermal_expansion_coefficient_alpha * gravity_norm) / specific_heat_Cp) * 1000.0) * ((depth) / 1000.0);;
    for (std::vector<Features::Interface *>::const_iterator it = features.begin(); it != features.end(); ++it)
      {
        temperature = (*it)->temperature(point,depth,gravity_norm,temperature);
      }

    return temperature;
  }

  bool
  World::composition(const std::array<double,2> &point, const double depth, const unsigned int composition_number) const
  {
    // turn it into a 3d coordinate and call the 3d temperature function
    AssertThrow(dim == 2, "This function can only be called when the cross section "
                "variable in the world builder file has been set.");

    Point<3> coord_3d(cross_section[0][0] + point[0] * surface_coord_conversions[0],
                      cross_section[0][1] + point[0] * surface_coord_conversions[1],
                      point[1]);

    return composition(coord_3d.get_array(), depth, composition_number);
  }

  bool
  World::composition(const std::array<double,3> &point_, const double depth, const unsigned int composition_number) const
  {
    Point<3> point(point_);
    double composition = 0;
    for (std::vector<Features::Interface *>::const_iterator it = features.begin(); it != features.end(); ++it)
      {
        composition = (*it)->composition(point,depth,composition_number, composition);
      }

    return composition;
  }

  WorldBuilder::CoordinateSystems::Interface &
  World::get_coordinate_system() const
  {
    return *coordinate_system;
  }
}

