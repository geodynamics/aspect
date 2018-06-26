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


#ifndef _world_builder_features_interface_h
#define _world_builder_features_interface_h


#include <boost/property_tree/ptree.hpp>

#include <world_builder/world.h>
#include <world_builder/point.h>


using boost::property_tree::ptree;


namespace WorldBuilder
{
  class World;

  namespace Features
  {

    class Interface
    {
      public:
        /**
         * constructor
         */
        Interface();

        /**
         * Destructor
         */
        ~Interface();

        /**
         * read in the world builder file
         */
        virtual
        void read(ptree &property_tree) = 0;

        /**
         * takes temperature and position and returns a temperature.
         */
        virtual
        double temperature(const Point<3> position,
                           const double depth,
                           const double gravity,
                           double temperature) const = 0;
        /**
         * Returns a value for the reqeusted composition (0 is not present,
         * 1 is present) based on the given position and
         */
        virtual
        double composition(const Point<3> position,
                           const double depth,
                           const unsigned int composition_number,
                           double temperature) const = 0;


      protected:
        WorldBuilder::World *world;

        std::string name;
        std::vector<Point<2> > coordinates;
        std::string temperature_submodule_name;
        std::string composition_submodule_name;

    };


    /**
     * factory function
     */
    Interface *
    create_feature(const std::string name, WorldBuilder::World *world);

  }
}

#endif
