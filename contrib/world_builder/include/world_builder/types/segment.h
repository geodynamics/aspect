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

#ifndef WORLD_BUILDER_TYPES_SEGMENT_H
#define WORLD_BUILDER_TYPES_SEGMENT_H


#include "world_builder/types/plugin_system.h"


namespace WorldBuilder
{
  class Parameters;

  namespace Types
  {

    /**
     * This class represents a segment value with documentation
     */
    class Segment : public Interface
    {
      public:
        /**
         * A constructor
         */
        Segment(const double default_length,
                const WorldBuilder::Point<2> &default_thickness,
                const WorldBuilder::Point<2> &default_top_truncation,
                const WorldBuilder::Point<2> &default_angle,
                const Types::Interface &temperature_plugin_system,
                const Types::Interface &composition_plugin_system,
                const Types::Interface &grains_plugin_system_);

        /**
         * A constructor for the load_entry function
         */
        Segment(double default_length,
                WorldBuilder::Point<2> default_thickness,
                WorldBuilder::Point<2> default_angle,
                std::string description);

        /**
         * Copy constructor
         */
        Segment(Segment const &other);


        /**
         * Destructor
         */
        ~Segment() override;

        /**
         * Todo
         */
        void write_schema(Parameters &prm,
                          const std::string &name,
                          const std::string &documentation) const override final;


        double value_length;
        double default_length;
        WorldBuilder::Point<2> value_thickness;
        WorldBuilder::Point<2> default_thickness;
        WorldBuilder::Point<2> default_top_truncation;
        WorldBuilder::Point<2> value_angle;
        WorldBuilder::Point<2> default_angle;
        std::unique_ptr<Types::Interface> temperature_plugin_system;
        std::unique_ptr<Types::Interface> composition_plugin_system;
        std::unique_ptr<Types::Interface> grains_plugin_system;

      protected:
        Segment *clone_impl() const override final
        {
          return new Segment(*this);
        };
      private:

    };
  } // namespace Types
} // namespace WorldBuilder

#endif
