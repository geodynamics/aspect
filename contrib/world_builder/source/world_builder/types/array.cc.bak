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
#include "world_builder/types/array.h"

#include "world_builder/parameters.h"

namespace WorldBuilder
{
  namespace Types
  {
    Array::Array(const Interface &type,
                 const unsigned int min_items_,
                 const unsigned int max_items_,
                 const bool unique_items_)
      :
      inner_type(type.get_type()),
      inner_type_ptr(type.clone()),
      required(false),
      min_items(min_items_),
      max_items(max_items_),
      unique_items(unique_items_)
    {
      this->type_name = Types::type::Array;
    }

    Array::Array(Array const &other)
      :
      inner_type(other.inner_type),
      inner_type_ptr(other.inner_type_ptr->clone()),
      required(other.required),
      min_items(other.min_items),
      max_items(other.max_items),
      unique_items(other.unique_items)
    {
      this->type_name = Types::type::Array;
    }

    Array::~Array ()
      = default;


    void
    Array::write_schema(Parameters &prm,
                        const std::string &name,
                        const std::string &documentation) const
    {
      using namespace rapidjson;
      Document &declarations = prm.declarations;
      const std::string &base = prm.get_full_json_path() + "/" + name;

      Pointer((base + "/type").c_str()).Set(declarations,"array");
      Pointer((base + "/minItems").c_str()).Set(declarations,min_items);
      Pointer((base + "/maxItems").c_str()).Set(declarations,max_items);
      Pointer((base + "/uniqueItems").c_str()).Set(declarations,unique_items);
      Pointer((base + "/description").c_str()).Set(declarations,documentation.c_str());

      prm.enter_subsection(name);
      {
        WBAssertThrow(this->inner_type_ptr != nullptr, "Internal error, inner pointer is NULL.");
        this->inner_type_ptr->write_schema(prm, "items", "");
      }
      prm.leave_subsection();


    }
  } // namespace Types
} // namespace WorldBuilder

