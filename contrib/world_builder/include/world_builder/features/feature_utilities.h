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


#ifndef WORLD_BUILDER_FEATURES_FEATURE_UTILITIES_H
#define WORLD_BUILDER_FEATURES_FEATURE_UTILITIES_H

#include <limits>
#include <vector>

#include "world_builder/assert.h"

namespace WorldBuilder
{
  namespace Features
  {
    namespace FeatureUtilities
    {
      enum class Operations
      {
        REPLACE,ADD,SUBTRACT,REPLACE_DEFINED_ONLY
      };

      /**
       * transforms string operations into enums.
       */
      Operations
      string_operations_to_enum(const std::string &operation);


      /**
       * Applies different operations such as replace, add and subtract to the original values
       */
      inline double
      apply_operation(const Operations operation,
                      const double old_value,
                      const double new_value)
      {
        switch (operation)
          {
            case Operations::REPLACE:
            case Operations::REPLACE_DEFINED_ONLY:
              return new_value;
              break;

            case Operations::ADD:
              return old_value + new_value;
              break;

            case Operations::SUBTRACT:
              return old_value - new_value;

            default:
              WBAssert(false,"Operation not found.");
          }

        return std::numeric_limits<double>::signaling_NaN();
      }

      /**
       * A struct that is used to hold additional values based on the output of
       *  the function distance_point_from_curved_planes().
       */
      struct AdditionalParameters
      {
        // The total length of all the segments at the location of the plane.
        double total_local_segment_length;

        // The local thickness of the segment at the location of the plane.
        double local_thickness;
      };


      /**
       * Add a string to a vector of strings, if the exact string isn't already
       * in the vector. Returns the location of the string in the final vector.
       */
      size_t
      add_vector_unique(std::vector<std::string> &vector,const std::string &add_string);

    } // namespace Utilities
  } // namespace Features
} // namespace WorldBuilder

#endif
