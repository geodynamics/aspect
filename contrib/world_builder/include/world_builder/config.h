/*
  Copyright (C) 2018-2021 by the authors of the World Builder code.

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

#ifndef WORLD_BUILDER_CONFIG_H_
#define WORLD_BUILDER_CONFIG_H_


#include <string>

namespace WorldBuilder
{
  struct Version
  {
    static const std::string MAJOR;
    static const std::string MINOR;
    static const std::string PATCH;
    static const std::string LABEL;

    static const std::string GIT_SHA1;
    static const std::string GIT_BRANCH;
    static const std::string GIT_DATE;
    static const std::string GIT_COMMIT_SUBJECT;
  };

  struct Data
  {
    static const std::string WORLD_BUILDER_SOURCE_DIR;
  };
}


#endif /* WORLD_BUILDER_CONFIG_H_ */
