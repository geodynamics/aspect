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
#include "world_builder/types/string.h"


#include "world_builder/parameters.h"

namespace WorldBuilder
{
  namespace Types
  {

    String::String(std::string default_value_)
      :
      default_value(std::move(default_value_)),

      restricted_values({})
    {
      this->type_name = Types::type::String;
    }

    String::String(std::string default_value_,
                   const std::string &restricted_value_)
      :
      default_value(std::move(default_value_)),
      restricted_values({restricted_value_})
    {
      this->type_name = Types::type::String;
    }

    String::String(String const &other)
      :
      value(other.value),
      default_value(other.default_value)
    {
      this->type_name = Types::type::String;
    }


    String::String(std::string default_value_,
                   std::vector<std::string> restricted_values_)
      :
      default_value(std::move(default_value_)),
      restricted_values(std::move(restricted_values_))
    {
      this->type_name = Types::type::String;
    }

    /*String::String(std::string default_value, std::string description)
      :
      value(default_value),
      default_value(default_value),
      description(description)
    {
      this->type_name = Types::type::String;
    }*/

    String::String(std::string value_,
                   std::string default_value_,
                   std::string description_)
      :
      value(std::move(value_)),
      default_value(std::move(default_value_)),
      description(std::move(description_))
    {
      this->type_name = Types::type::String;
    }

    String::~String ()
      = default;

    void
    String::write_schema(Parameters &prm,
                         const std::string &name,
                         const std::string &documentation) const
    {
      using namespace rapidjson;
      Document &declarations = prm.declarations;
      const std::string base = prm.get_full_json_path() + "/" + name;
      Pointer((base + "/default value").c_str()).Set(declarations,default_value.c_str());
      Pointer((base + "/type").c_str()).Set(declarations,"string");
      Pointer((base + "/description").c_str()).Set(declarations,documentation.c_str());
      for (unsigned int i = 0; i < restricted_values.size(); ++i)
        {
          if (!restricted_values[i].empty())
            {
              if (i == 0 && Pointer((base + "/enum").c_str()).Get(declarations) == nullptr)
                {
                  // The enum array doesn't exist yet, so we create it and fill it.
                  Pointer((base + "/enum/0").c_str()).Create(declarations);
                  Pointer((base + "/enum/0").c_str()).Set(declarations, restricted_values[i].c_str());
                }
              else
                {
                  // The enum array already exist yet, so we add an element to the end.
                  Pointer((base + "/enum/-").c_str()).Set(declarations, restricted_values[i].c_str());
                }
            }
        }
    }

  } // namespace Types
} // namespace WorldBuilder

