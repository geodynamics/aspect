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


#ifndef _aspect_world_builder_world_h
#define _aspect_world_builder_world_h


#include <boost/property_tree/ptree.hpp>

#include <aspect/world_builder/feature.h>

using boost::property_tree::ptree;

namespace aspect
{
  namespace WorldBuilder
  {

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

        double temperature(const std::array<double, 2> point) const;

        double temperature(const std::array<double, 3> point) const;

      private:
        const char path_seperator = '.';
        /**
         * These are the top level parameters
         */
        std::vector<double> surface_rotation_point;
        double surface_rotation_angle;
        unsigned int minimum_parts_per_distance;
        double minimum_distance_points;

        /**
         * contains all the plugins.
         * todo: make a unique or shared pointer?
         */
        std::vector<Feature *> features;



    };
  }
}

#endif
