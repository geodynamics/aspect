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

cmake_minimum_required(VERSION 3.13.4)

message(STATUS "===== Configuring ASPECT cookbooks =================")

set(Aspect_DIR ${CMAKE_BINARY_DIR})

# Step 1: Find all separate plugins.
file(GLOB_RECURSE _cmakelists
     ${CMAKE_CURRENT_SOURCE_DIR}/*/CMakeLists.txt)

# Step 2: Add all subdirectories of plugins.
foreach(_file ${_cmakelists})
  get_filename_component(_dir ${_file} DIRECTORY)
  add_subdirectory(${_dir})
endforeach()
