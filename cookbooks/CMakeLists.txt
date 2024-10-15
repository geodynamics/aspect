# Copyright (C) 2023 - 2024 by the authors of the ASPECT code.
#
# This file is part of ASPECT.
#
# ASPECT is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# ASPECT is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ASPECT; see the file doc/COPYING.  If not see
# <http://www.gnu.org/licenses/>.

CMAKE_MINIMUM_REQUIRED(VERSION 3.13.4)

MESSAGE(STATUS "===== Configuring ASPECT cookbooks =================")

SET(Aspect_DIR ${CMAKE_BINARY_DIR})

# Step 1: Find all separate plugins.
FILE(GLOB_RECURSE _cmakelists
     ${CMAKE_CURRENT_SOURCE_DIR}/*/CMakeLists.txt)

# Step 2: Add all subdirectories of plugins.
FOREACH(_file ${_cmakelists})
  GET_FILENAME_COMPONENT(_dir ${_file} DIRECTORY)
  ADD_SUBDIRECTORY(${_dir})
ENDFOREACH()
