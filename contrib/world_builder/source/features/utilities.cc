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

#include <world_builder/features/utilities.h>

namespace WorldBuilder
{
  namespace Features
  {
    namespace Utilities
    {
      Operations
      string_operations_to_enum(const std::string &operation)
      {
        if (operation == "add") return Operations::ADD;
        if (operation == "subtract") return Operations::SUBTRACT;

        WBAssert(operation == "replace", "Could not find operation: " << operation << ".");
        return Operations::REPLACE;
      }
    }
  }
}
