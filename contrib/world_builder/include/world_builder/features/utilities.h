/*
  Copyright (C) 2020 by the authors of the World Builder code.

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


#ifndef _world_builder_features_utilities_h
#define _world_builder_features_utilities_h

#include <string>
#include <limits>

#include <world_builder/assert.h>

namespace WorldBuilder
{
  namespace Features
  {
    namespace Utilities
    {
      enum class Operations
      {
        REPLACE,ADD,SUBTRACT
      };

      /**
       * transforms string operations into enums.
       */
      Operations
      string_operations_to_enum(const std::string &operation);


      /**
       * Applies different opertions such as replace, add and subtract to the original values
       */
      inline double
      apply_operation(const Utilities::Operations operation,
                      const double old_value,
                      const double new_value)
      {
        switch (operation)
          {
            case Utilities::Operations::REPLACE:
              return new_value;
              break;

            case Utilities::Operations::ADD:
              return old_value + new_value;
              break;

            case Utilities::Operations::SUBTRACT:
              return old_value - new_value;

            default:
              WBAssert(false,"Operation not found.");
          }

        return std::numeric_limits<double>::signaling_NaN();
      }
    }
  }
}

#endif