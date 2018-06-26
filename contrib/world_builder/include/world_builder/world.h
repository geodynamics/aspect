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


#ifndef _world_builder_world_h
#define _world_builder_world_h


#include <boost/property_tree/ptree.hpp>

#include <world_builder/features/interface.h>
#include <world_builder/coordinate_systems/interface.h>

using boost::property_tree::ptree;


  namespace WorldBuilder
  {
    namespace Features
    {
      class Interface;
    }

    class World
    {
      public:
        /**
         * constructor
         */
        World(std::string filename);

        /**
         * Destructor
         */
        ~World();

        /**
         * read in the world builder file
         */
        void read(ptree &property_tree);

        double temperature(const std::array<double, 2> point, const double depth, const double gravity_norm) const;

        double temperature(const std::array<double, 3> point, const double depth, const double gravity_norm) const;


        double composition(const std::array<double, 2> point, const double depth, const unsigned int composition_number) const;

        double composition(const std::array<double, 3> point, const double depth, const unsigned int composition_number) const;

        /**
         * returs a pointer to the coordinate system
         */
        WorldBuilder::CoordinateSystems::Interface *get_coordinate_system() const;

      private:
        const char path_seperator = '.';
        /**
         * These are the top level parameters
         */
        std::vector<double> surface_rotation_point;
        double surface_rotation_angle;
        unsigned int minimum_parts_per_distance_unit;
        double minimum_distance_points;

        /**
         * adiabatic parameters
         */
        double potential_mantle_temperature;
        double thermal_expansion_coefficient_alfa;
        double specific_heat_Cp;

        /**
         * contains all the plugins.
         * todo: make a unique or shared pointer?
         */
        std::vector<WorldBuilder::Features::Interface *> features;

        // coordinate system
        WorldBuilder::CoordinateSystems::Interface *coordinate_system;



    };
  }

#endif
