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

#ifndef WORLD_BUILDER_TYPES_VALUE_AT_POINTS_H
#define WORLD_BUILDER_TYPES_VALUE_AT_POINTS_H

#include "world_builder/types/point.h"
#include <vector>


namespace WorldBuilder
{
  class Parameters;

  namespace Types
  {
    /**
     * This class represents a depth surface value with documentation
     */
    class ValueAtPoints : public Interface
    {
      public:
        /**
         * A constructor
         */
        ValueAtPoints(const double default_value,
                      uint64_t max_values_in_array,
                      std::vector<Point<2>> default_points_ = std::vector<Point<2>>());

        /**
         * Copy constructor
         */
        ValueAtPoints(ValueAtPoints const &other);


        /**
         * Destructor
         */
        ~ValueAtPoints() override;

        /**
         * Write schema
         */
        void write_schema(Parameters &prm,
                          const std::string &name,
                          const std::string &documentation) const override final;

        double default_value;
        uint64_t max_values_in_array;
        std::vector<Point<2> > default_points;

      protected:
        ValueAtPoints *clone_impl() const override final
        {
          return new ValueAtPoints(*this);
        };
      private:

    };
  } // namespace Types

} // namespace WorldBuilder

#endif
