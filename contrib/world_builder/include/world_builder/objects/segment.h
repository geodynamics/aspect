/*
  Copyright (C) 2018-2024 by the authors of the World Builder code.

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

#ifndef WORLD_BUILDER_OBJECTS_SEGMENT_H
#define WORLD_BUILDER_OBJECTS_SEGMENT_H


#include "world_builder/types/plugin_system.h"


namespace WorldBuilder
{
  class Parameters;

  namespace Objects
  {

    /**
      * This class represents an actual segment
      */
    template <class A, class B, class C>
    class Segment
    {
      public:

        /**
         * A constructor for the clone and set_entry function
         */
        Segment(const double default_length,
                const WorldBuilder::Point<2> &default_thickness,
                const WorldBuilder::Point<2> &default_top_truncation,
                const WorldBuilder::Point<2> &default_angle,
                std::vector<std::shared_ptr<A> > temperature_systems,
                std::vector<std::shared_ptr<B> > composition_systems,
                std::vector<std::shared_ptr<C> > grains_systems);

        /**
         * Copy constructor
         */
        Segment(Segment const &other);

        /**
         * Destructor
         */
        ~Segment();

        double value_length;
        double default_length;
        WorldBuilder::Point<2> value_thickness;
        WorldBuilder::Point<2> value_top_truncation;
        WorldBuilder::Point<2> value_angle;
        std::vector<std::shared_ptr<A> > temperature_systems;
        std::vector<std::shared_ptr<B> > composition_systems;
        std::vector<std::shared_ptr<C> > grains_systems;

      protected:
      private:

    };
  } // namespace Objects
} // namespace WorldBuilder

#endif
