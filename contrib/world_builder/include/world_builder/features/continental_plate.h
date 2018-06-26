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



#ifndef _world_feature_features_continental_plate_h
#define _world_feature_features_continental_plate_h

#include <world_builder/features/interface.h>
#include <world_builder/world.h>


namespace WorldBuilder
{
  namespace Features
  {

    class ContinentalPlate : public Interface
    {
      public:
        /**
         * constructor
         */
        ContinentalPlate(WorldBuilder::World *world);

        /**
         * Destructor
         */
        ~ContinentalPlate();

        /**
         * Read in the world builder file
         */
        virtual
        void read(ptree &property_tree);

        /**
         * Returns a temperature based on the given position
         */
        virtual
        double temperature(const Point<3> position,
                           const double depth,
                           const double gravity,
                           double temperature) const;

        /**
         * Returns a value for the reqeusted composition (0 is not present,
         * 1 is present) based on the given position and
         */
        virtual
        double composition(const Point<3> position,
                           const double depth,
                           const unsigned int composition_number,
                           double temperature) const;



      private:
        // local parameters
        double temperature_submodule_depth;
        double temperature_submodule_temperature;
        double composition_submodule_depth;
        unsigned int composition_submodule_composition;

    };
  }
}

#endif
