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

#include <boost/property_tree/ptree.hpp>

using boost::property_tree::ptree;

#ifndef _aspect_world_feature_world_h
#define _aspect_world_feature_world_h

namespace aspect
{
  namespace WorldBuilder
  {
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
          double temperature(const std::array<double,3> position, double temperature) const = 0;


        protected:
          std::string name;
          std::vector<std::vector<double> > coordinates;
          std::string temperature_submodule_name;
          std::string composition_submodule_name;

      };


      /**
       * factory function
       */
      Interface *
      create_feature(const std::string name);

    }
  }
}

#endif
