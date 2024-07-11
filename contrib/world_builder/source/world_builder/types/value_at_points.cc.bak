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
#include <utility>

#include "world_builder/types/value_at_points.h"

#include "world_builder/parameters.h"


namespace WorldBuilder
{
  namespace Types
  {

    ValueAtPoints::ValueAtPoints(const double default_value_,
                                 uint64_t max_values_in_array_,
                                 std::vector<Point<2>> default_points_)
      :
      default_value(default_value_),
      max_values_in_array(max_values_in_array_),
      default_points(std::move(std::move(default_points_)))
    {
      this->type_name = Types::type::ValueAtPoints;

      WBAssert(default_points_.size() == 0, "default points are not implemented.");
    }

    ValueAtPoints::ValueAtPoints(ValueAtPoints const &other)
      :
      default_value(other.default_value),
      max_values_in_array(other.max_values_in_array),
      default_points(other.default_points)
    {
      this->type_name = Types::type::ValueAtPoints;
    }


    ValueAtPoints::~ValueAtPoints ()
      = default;


    void
    ValueAtPoints::write_schema(Parameters &prm,
                                const std::string &name,
                                const std::string &documentation) const
    {
      using namespace rapidjson;
      prm.enter_subsection(name);
      {
        Document &declarations = prm.declarations;
        const std::string base = prm.get_full_json_path();

        Pointer((base + "/type").c_str()).Set(declarations,"array");
        Pointer((base + "/additionalProperties").c_str()).Set(declarations,false);
        Pointer((base + "/minItems").c_str()).Set(declarations,1);
        Pointer((base + "/maxItems").c_str()).Set(declarations, max_values_in_array);
        Pointer((base + "/description").c_str()).Set(declarations,documentation.c_str());

        {
          Pointer((base + "/items/anyOf/0/type").c_str()).Set(declarations,"number");
          Pointer((base + "/items/anyOf/0/default value").c_str()).Set(declarations,default_value);

          Pointer((base + "/items/anyOf/1/type").c_str()).Set(declarations,"array");
          Pointer((base + "/items/anyOf/1/minItems").c_str()).Set(declarations,1);
          Pointer((base + "/items/anyOf/1/maxItems").c_str()).Set(declarations,std::numeric_limits<unsigned int>::max());

          Pointer((base + "/items/anyOf/1/items/type").c_str()).Set(declarations,"array");
          Pointer((base + "/items/anyOf/1/items/minItems").c_str()).Set(declarations,1);
          Pointer((base + "/items/anyOf/1/items/maxItems").c_str()).Set(declarations, max_values_in_array);

          Pointer((base + "/items/anyOf/1/items/items/type").c_str()).Set(declarations,"number");
        }
      }
      prm.leave_subsection();
    }

  } // namespace Types
} // namespace WorldBuilder


