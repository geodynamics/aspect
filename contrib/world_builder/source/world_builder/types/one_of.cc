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
#include "world_builder/types/one_of.h"

#include "world_builder/parameters.h"

namespace WorldBuilder
{
  namespace Types
  {
    OneOf::OneOf(const Interface &type_1,
                 const Interface &type_2)
    {
      this->type_name = Types::type::OneOf;
      inner_types_ptr.emplace_back(type_1.clone());
      inner_types_ptr.emplace_back(type_2.clone());

    }

    OneOf::OneOf(OneOf const &other)
    {
      this->type_name = Types::type::OneOf;
      for (const auto &i : other.inner_types_ptr)
        {
          inner_types_ptr.emplace_back(i->clone());
        }
    }

    OneOf::~OneOf ()
      = default;


    void
    OneOf::write_schema(Parameters &prm,
                        const std::string &name,
                        const std::string &documentation) const
    {
      using namespace rapidjson;
      Document &declarations = prm.declarations;
      const std::string &base = prm.get_full_json_path() + "/" + name;

      prm.enter_subsection(name);
      {
        Pointer((base + "/description").c_str()).Set(declarations,documentation.c_str());
        prm.enter_subsection("oneOf");
        {
          for (size_t i = 0; i < inner_types_ptr.size(); i++)
            {
              WBAssertThrow(inner_types_ptr[i] != nullptr, "Internal error, inner pointer is NULL.");
              inner_types_ptr[i]->write_schema(prm, std::to_string(i), "");
            }
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();


    }
  } // namespace Types
} // namespace WorldBuilder

