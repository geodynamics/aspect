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
#ifndef WORLD_BUILDER_ASSERT_H_
#define WORLD_BUILDER_ASSERT_H_

#include <sstream>

namespace WorldBuilder
{
#ifndef NDEBUG
#   define WBAssert(condition, message) \
  do { \
      if (! (condition)) { \
          std::stringstream smessage; \
          smessage << "Assert `" #condition "` failed in " << __FILE__ \
                   << " at line " << __LINE__ << ": " << message << std::endl << std::endl << "Error not recoverable, aborting program."; \
          throw std::runtime_error(smessage.str()); \
        } \
    } while (false)
#else
#   define WBAssert(condition, message) do { } while (false)
#endif

#   define WBAssertThrow(condition, message) \
  do { \
      if (! (condition)) { \
          std::stringstream smessage; \
          smessage << "AssertThrow `" #condition "` failed in " << __FILE__ \
                   << " at line " << __LINE__ << ": " << message << std::endl << std::endl << "Error not recoverable, aborting program."; \
          throw std::runtime_error(smessage.str()); \
        } \
    } while (false)


#   define WBAssertThrowExc(condition, exc, message) \
  do { \
      if (! (condition)) { \
          exc \
          std::stringstream smessage; \
          smessage << "AssertThrow `" #condition "` failed in " << __FILE__ \
                   << " at line " << __LINE__ << ": " << message << std::endl; \
          throw std::runtime_error(smessage.str()); \
        } \
    } while (false)
} // namespace WorldBuilder

#endif
