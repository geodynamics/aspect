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


#ifndef _aspect_world_builder_features_interface_h
#define _aspect_world_builder_features_interface_h


#include <aspect/world_builder/world.h>

#include <boost/property_tree/ptree.hpp>

using boost::property_tree::ptree;

namespace aspect
{
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
          double temperature(const std::array<double,3> position,
                             const double depth,
                             const double gravity,
                             double temperature) const = 0;


        protected:
          WorldBuilder::World* world;

          std::string name;
          std::vector<std::array<double,2> > coordinates;
          std::string temperature_submodule_name;
          std::string composition_submodule_name;

      };


      /**
       * factory function
       */
      Interface *
      create_feature(const std::string name, WorldBuilder::World* world);

    }
  }
}

#endif
