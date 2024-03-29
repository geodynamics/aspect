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

#include "world_builder/types/object.h"


#include "world_builder/parameters.h"

namespace WorldBuilder
{
  namespace Types
  {
    Object::Object(std::vector<std::string> required_,
                   const bool additional_properties_)
      :
      required(std::move(required_)),
      additional_properties(additional_properties_)
    {
      this->type_name = Types::type::Object;
    }

    Object::Object(Object const &other)
      :
      required(other.required),
      additional_properties(other.additional_properties)
    {
      this->type_name = Types::type::Object;
    }

    Object::~Object ()
      = default;

    void
    Object::write_schema(Parameters &prm,
                         const std::string & /*name*/,
                         const std::string &documentation) const
    {
      // The object has to be defined just above the properties section.
      // because it is very inconvenient to make this object there, we
      // take care here to do that.
      WBAssertThrow(!prm.path.empty() && prm.path.back() == "properties",
                    "The Object must be defined in a subsection called properties.");
      prm.leave_subsection();
      {
        using namespace rapidjson;
        Document &declarations = prm.declarations;
        const std::string path = prm.get_full_json_path();

        Pointer((path + "/type").c_str()).Set(declarations,"object");
        Pointer((path + "/description").c_str()).Set(declarations,documentation.c_str());
        Pointer((path + "/additionalProperties").c_str()).Set(declarations,additional_properties);

        if (!required.empty())
          {
            for (unsigned int i = 0; i < required.size(); ++i)
              {
                if (i == 0 && Pointer((path + "/required").c_str()).Get(declarations) == nullptr)
                  {
                    // The required array doesn't exist yet, so we create it and fill it.
                    Pointer((path + "/required/0").c_str()).Create(declarations);
                    Pointer((path + "/required/0").c_str()).Set(declarations, required[i].c_str());
                  }
                else
                  {
                    // The required array already exist yet, so we add an element to the end.
                    Pointer((path + "/required/-").c_str()).Set(declarations, required[i].c_str());
                  }
              }
          }
      }
      prm.enter_subsection("properties");
    }
  } // namespace Types
} // namespace WorldBuilder

