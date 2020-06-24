/*
  Copyright (C) 2018 - 2020 by the authors of the World Builder code.

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

#include "rapidjson/pointer.h"

#include <world_builder/config.h>
#include <world_builder/world.h>
#include <world_builder/utilities.h>
#include <world_builder/assert.h>
#include <world_builder/point.h>
#include <world_builder/nan.h>
#include <world_builder/parameters.h>
#include <world_builder/coordinate_systems/interface.h>
#include <world_builder/types/interface.h>

#include <world_builder/types/array.h>
#include <world_builder/types/bool.h>
#include <world_builder/types/double.h>
#include <world_builder/types/string.h>
#include <world_builder/types/object.h>
#include <world_builder/types/point.h>
#include <world_builder/types/unsigned_int.h>
#include <world_builder/types/plugin_system.h>


namespace WorldBuilder
{


  using namespace Utilities;

  World::World(std::string filename, bool has_output_dir, std::string output_dir, unsigned long random_number_seed)
    :
    parameters(*this),
    surface_coord_conversions(invalid),
    dim(NaN::ISNAN),
    random_number_engine(random_number_seed)
  {
    this->declare_entries(parameters);

    parameters.initialize(filename, has_output_dir, output_dir);

    this->parse_entries(parameters);
  }

  World::~World()
  {}

  void World::declare_entries(Parameters &prm)
  {
    prm.enter_subsection("properties");
    {
      prm.declare_entry("", Types::Object({"version", "features"}), "Root object");

      prm.declare_entry("version", Types::String(""),"The major and minor version number for which the input file was written.");

      prm.declare_entry("cross section", Types::Array(Types::Point<2>(),2,2),"This is an array of two points along where the cross section is taken");

      prm.declare_entry("potential mantle temperature", Types::Double(1600),
                        "The potential temperature of the mantle at the surface in Kelvin.");
      prm.declare_entry("surface temperature", Types::Double(293.15),
                        "The temperature at the surface in Kelvin.");
      prm.declare_entry("force surface temperature", Types::Bool(false),
                        "Force the provided surface temperature to be set at the surface");
      prm.declare_entry("thermal expansion coefficient", Types::Double(3.5e-5),
                        "The thermal expansion coefficient in $K^{-1}$.");
      prm.declare_entry("specific heat", Types::Double(1250),
                        "The specific heat in $J kg^{-1} K^{-1}.$");
      prm.declare_entry("thermal diffusivity", Types::Double(0.804e-6),
                        "The thermal diffusivity in $m^{2} s^{-1}$.");

      prm.declare_entry("maximum distance between coordinates",Types::Double(0),
                        "This enforces a maximum distance (in degree for spherical coordinates "
                        "or meter in cartesian coordinates) between coordinates in the model. "
                        "If the distance is larger, extra points are added by interpolation. "
                        "Requires interpolation to be not 'none'.");

      prm.declare_entry("interpolation",Types::String("none"),
                        "What type of interpolation should be used to enforce the minimum points per "
                        "distance parameter. Options are none, linear and monotone spline.");


      prm.declare_entry("coordinate system", Types::PluginSystem("cartesian", CoordinateSystems::Interface::declare_entries, {"model"}, false),"A coordinate system. Cartesian or spherical.");
      prm.declare_entry("features", Types::PluginSystem("",Features::Interface::declare_entries, {"model", "coordinates"}),"A list of features.");


    }
    prm.leave_subsection();

  }


  void World::parse_entries(Parameters &prm)
  {
    using namespace rapidjson;

    /**
     * First load the major version number in the file and check the major
     * version number of the program.
     */

    WBAssertThrow((Version::MAJOR == "0"
                   && prm.get<std::string>("version") == Version::MAJOR + "." + Version::MINOR)
                  || (Version::MAJOR != "0"
                      && prm.get<std::string>("version") == Version::MAJOR),
                  "The major and minor version combination (for major version 0) or the major "
                  "version (for major versions after 0) for which is input file was written "
                  "is not the same as the version of the World Builder you are running. This means "
                  "That there may have been incompatible changes made between the versions. \n\n"
                  "Verify those changes and wheter they affect your model. If this is not "
                  "the case, adjust the version number in the input file. \n\nThe provided version "
                  "number is \"" << prm.get<std::string>("version") << "\", while the used world builder "
                  "has  (major.minor) version \"" << Version::MAJOR << "." << Version::MINOR << "\". "
                  "If you created this file from scratch, fill set the version number to \"" <<
                  Version::MAJOR << "." << Version::MINOR << "\" to continue. If you got the world builder "
                  "file from somewhere, make sure that the output is what you expect it to be, because "
                  "backwards incompatible changes may have been made to the code.");

    /**
     * Seconly load the coordinate system parameters.
     */
    prm.coordinate_system = prm.get_unique_pointer<CoordinateSystems::Interface>("coordinate system");
    prm.coordinate_system->parse_entries(prm);

    prm.get_unique_pointers<Features::Interface>("features",prm.features);

    const bool set_cross_section = prm.check_entry("cross section");

    const CoordinateSystem coordinate_system = prm.coordinate_system->natural_coordinate_system();

    if (set_cross_section == true)
      {
        dim = 2;
        std::vector<Point<2> > cross_section_natural = prm.get_vector<Point<2> >("cross section");

        WBAssertThrow(cross_section_natural.size() == 2, "The cross section should contain two points, but it contains "
                      << cross_section.size() << " points.");

        for (auto it : cross_section_natural)
          cross_section.push_back(it *  (coordinate_system == spherical ? const_pi / 180.0 : 1.0));


        /**
         * pre-compute stuff for the cross section
         */
        surface_coord_conversions = cross_section[0]-cross_section[1];
        surface_coord_conversions *= -1/(surface_coord_conversions.norm());


      }
    else
      {
        dim = 3;
      }

    /**
     * Temperature parameters.
     */
    potential_mantle_temperature = prm.get<double>("potential mantle temperature");
    surface_temperature = prm.get<double>("surface temperature");
    surface_temperature = prm.get<double>("surface temperature");
    force_surface_temperature = prm.get<bool>("force surface temperature");
    thermal_expansion_coefficient = prm.get<double>("thermal expansion coefficient");
    specific_heat = prm.get<double>("specific heat");
    thermal_diffusivity = prm.get<double>("thermal diffusivity");

    /**
     * Model discretiation paramters
     */
    maximum_distance_between_coordinates = prm.get<double>("maximum distance between coordinates");
    interpolation = prm.get<std::string>("interpolation");

    /**
     * Now load the features. Some features use for example temperature values,
     * so it is important that this is parsed the last.
     */
    prm.enter_subsection("features");
    {
      for (unsigned int i = 0; i < prm.features.size(); ++i)
        {
          prm.enter_subsection(std::to_string(i));
          {
            prm.features[i]->parse_entries(prm);
          }
          prm.leave_subsection();
        }
    }
    prm.leave_subsection();
  }

  double
  World::temperature(const std::array<double,2> &point,
                     const double depth,
                     const double gravity_norm) const
  {
    // turn it into a 3d coordinate and call the 3d temperature function
    WBAssertThrow(dim == 2, "This function can only be called when the cross section "
                  "variable in the world builder file has been set. Dim is "
                  << dim << ".");

    const CoordinateSystem coordinate_system = this->parameters.coordinate_system->natural_coordinate_system();

    Point<2> point_natural(point[0], point[1],coordinate_system);
    if (coordinate_system == spherical)
      {
        point_natural[1] = std::sqrt(point[0]*point[0]+point[1]*point[1]);
        point_natural[0] = std::atan2(point[1],point[0]);
      }

    Point<3> coord_3d(coordinate_system);
    if (coordinate_system == spherical)
      {
        coord_3d[0] = point_natural[1];
        coord_3d[1] = cross_section[0][0] + point_natural[0] * surface_coord_conversions[0];
        coord_3d[2] = cross_section[0][1] + point_natural[0] * surface_coord_conversions[1];
      }
    else
      {
        coord_3d[0] = cross_section[0][0] + point_natural[0] * surface_coord_conversions[0];
        coord_3d[1] = cross_section[0][1] + point_natural[0] * surface_coord_conversions[1];
        coord_3d[2] = point_natural[1];
      }


    std::array<double, 3> point_3d_cartesian = this->parameters.coordinate_system->natural_to_cartesian_coordinates(coord_3d.get_array());

    return temperature(point_3d_cartesian, depth, gravity_norm);
  }

  double
  World::temperature(const std::array<double,3> &point_,
                     const double depth,
                     const double gravity_norm) const
  {
    // We receive the cartesian points from the user.
    Point<3> point(point_,cartesian);

    if (std::fabs(depth) < 2.0 * std::numeric_limits<double>::epsilon() && force_surface_temperature == true)
      return this->surface_temperature;

    double temperature = potential_mantle_temperature *
                         std::exp(((thermal_expansion_coefficient * gravity_norm) /
                                   specific_heat) * depth);


    for (auto &&it : parameters.features)
      {
        temperature = it->temperature(point,depth,gravity_norm,temperature);

        WBAssert(!std::isnan(temperature), "Temparture is not a number: " << temperature
                 << ", based on a feature with the name " << it->get_name());
        WBAssert(std::isfinite(temperature), "Temparture is not a finite: " << temperature
                 << ", based on a feature with the name " << it->get_name());
      }

    WBAssert(!std::isnan(temperature), "Temparture is not a number: " << temperature);
    WBAssert(std::isfinite(temperature), "Temparture is not a finite: " << temperature);

    return temperature;
  }

  double
  World::composition(const std::array<double,2> &point,
                     const double depth,
                     const unsigned int composition_number) const
  {
    // turn it into a 3d coordinate and call the 3d temperature function
    WBAssertThrow(dim == 2, "This function can only be called when the cross section "
                  "variable in the world builder file has been set. Dim is "
                  << dim << ".");

    const CoordinateSystem coordinate_system = this->parameters.coordinate_system->natural_coordinate_system();

    Point<2> point_natural(point[0], point[1],coordinate_system);
    if (coordinate_system == spherical)
      {
        point_natural[1] = std::sqrt(point[0]*point[0]+point[1]*point[1]);
        point_natural[0] = std::atan2(point[1],point[0]);
      }

    Point<3> coord_3d(coordinate_system);
    if (coordinate_system == spherical)
      {
        coord_3d[0] = point_natural[1];
        coord_3d[1] = cross_section[0][0] + point_natural[0] * surface_coord_conversions[0];
        coord_3d[2] = cross_section[0][1] + point_natural[0] * surface_coord_conversions[1];
      }
    else
      {
        coord_3d[0] = cross_section[0][0] + point_natural[0] * surface_coord_conversions[0];
        coord_3d[1] = cross_section[0][1] + point_natural[0] * surface_coord_conversions[1];
        coord_3d[2] = point_natural[1];
      }

    std::array<double, 3> point_3d_cartesian = this->parameters.coordinate_system->natural_to_cartesian_coordinates(coord_3d.get_array());

    return composition(point_3d_cartesian, depth, composition_number);
  }

  double
  World::composition(const std::array<double,3> &point_,
                     const double depth,
                     const unsigned int composition_number) const
  {
    // We receive the cartesian points from the user.
    Point<3> point(point_,cartesian);
    double composition = 0;
    for (auto &&it : parameters.features)
      {
        composition = it->composition(point,depth,composition_number, composition);

        WBAssert(!std::isnan(composition), "Composition is not a number: " << composition
                 << ", based on a feature with the name " << it->get_name());
        WBAssert(std::isfinite(composition), "Composition is not a finite: " << composition
                 << ", based on a feature with the name " << it->get_name());
      }

    WBAssert(!std::isnan(composition), "Composition is not a number: " << composition);
    WBAssert(std::isfinite(composition), "Composition is not a finite: " << composition);


    return composition;
  }



  WorldBuilder::grains
  World::grains(const std::array<double,2> &point,
                const double depth,
                const unsigned int composition_number,
                size_t number_of_grains) const
  {
    // turn it into a 3d coordinate and call the 3d temperature function
    WBAssertThrow(dim == 2, "This function can only be called when the cross section "
                  "variable in the world builder file has been set. Dim is "
                  << dim << ".");

    const CoordinateSystem coordinate_system = this->parameters.coordinate_system->natural_coordinate_system();

    Point<2> point_natural(point[0], point[1],coordinate_system);
    if (coordinate_system == spherical)
      {
        point_natural[1] = std::sqrt(point[0]*point[0]+point[1]*point[1]);
        point_natural[0] = std::atan2(point[1],point[0]);
      }

    Point<3> coord_3d(coordinate_system);
    if (coordinate_system == spherical)
      {
        coord_3d[0] = point_natural[1];
        coord_3d[1] = cross_section[0][0] + point_natural[0] * surface_coord_conversions[0];
        coord_3d[2] = cross_section[0][1] + point_natural[0] * surface_coord_conversions[1];
      }
    else
      {
        coord_3d[0] = cross_section[0][0] + point_natural[0] * surface_coord_conversions[0];
        coord_3d[1] = cross_section[0][1] + point_natural[0] * surface_coord_conversions[1];
        coord_3d[2] = point_natural[1];
      }

    std::array<double, 3> point_3d_cartesian = this->parameters.coordinate_system->natural_to_cartesian_coordinates(coord_3d.get_array());

    return grains(point_3d_cartesian, depth, composition_number,number_of_grains);
  }

  WorldBuilder::grains
  World::grains(const std::array<double,3> &point_,
                const double depth,
                const unsigned int composition_number,
                size_t number_of_grains) const
  {
    // We receive the cartesian points from the user.
    Point<3> point(point_,cartesian);
    WorldBuilder::grains grains;
    grains.sizes.resize(number_of_grains,0);
    grains.rotation_matrices.resize(number_of_grains);
    for (std::vector<std::unique_ptr<Features::Interface> >::const_iterator it = parameters.features.begin(); it != parameters.features.end(); ++it)
      {
        grains = (*it)->grains(point,depth,composition_number, grains);

        /*WBAssert(!std::isnan(composition), "Composition is not a number: " << composition
                 << ", based on a feature with the name " << (*it)->get_name());
        WBAssert(std::isfinite(composition), "Composition is not a finite: " << composition
                 << ", based on a feature with the name " << (*it)->get_name());*/
      }

    /*WBAssert(!std::isnan(composition), "Composition is not a number: " << composition);
    WBAssert(std::isfinite(composition), "Composition is not a finite: " << composition);*/


    return grains;
  }

  std::mt19937 &
  World::get_random_number_engine()
  {
    return random_number_engine;
  }

}

