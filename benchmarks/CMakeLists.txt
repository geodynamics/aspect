# ------------------------------------------------------
#
# SPDX-License-Identifier: GPL-2.0-or-later
#
# SPDX-FileCopyrightText: Copyright (C) 2023-2025 by the ASPECT authors.
#
# This file is part of ASPECT.
#
# Detailed license information governing the source code and
# contributions can be found in the folder LICENSES and in
# LICENSE and CONTRIBUTING.md at the top level directory.
#
# ------------------------------------------------------

cmake_minimum_required(VERSION 3.13.4)

message(STATUS "===== Configuring ASPECT benchmarks ================")

set(Aspect_DIR ${CMAKE_BINARY_DIR})

# Step 1: Find all separate plugins.
file(GLOB_RECURSE _cmakelists
     ${CMAKE_CURRENT_SOURCE_DIR}/*/CMakeLists.txt)

# Step 2: Add all subdirectories of plugins.
foreach(_file ${_cmakelists})
  get_filename_component(_dir ${_file} DIRECTORY)
  add_subdirectory(${_dir})
endforeach()
