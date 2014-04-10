# Copyright (C) 2013, 2014 by the authors of the ASPECT code.
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

#  $Id$


MESSAGE("====================================================")
MESSAGE("============ Configuring ASPECT ====================")
MESSAGE("====================================================")

# Set the name of the project and target:
SET(TARGET "aspect")

FILE(GLOB_RECURSE TARGET_SRC  "source/*.cc")
INCLUDE_DIRECTORIES(include)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

FIND_PACKAGE(deal.II 8.0 QUIET
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()

IF(NOT DEAL_II_WITH_TRILINOS OR NOT DEAL_II_WITH_P4EST)
  MESSAGE(FATAL_ERROR
    "\nAspect requires a deal.II installation built with support for Trilinos and p4est but one or both of these appears to be missing!\n"
    )
  IF(NOT DEAL_II_WITH_TRILINOS)
    MESSAGE(FATAL_ERROR
      "\n-- deal.II was built without support for Trilinos!\n"
      )
  ENDIF()

  IF(NOT DEAL_II_WITH_P4EST)
    MESSAGE(FATAL_ERROR
      "\n-- deal.II was built without support for p4est!\n"
      )
  ENDIF()
ENDIF()

DEAL_II_INITIALIZE_CACHED_VARIABLES()
PROJECT(${TARGET})


# Pass down the source directory to the sources. This can be used
# to hard-code the location of data files, such as in
# $ASPECT_SOURCE_DIR/data/velocity-boundary-conditions/gplates/*
FOREACH(_source_file ${TARGET_SRC})
  SET_PROPERTY(SOURCE ${_source_file}
    APPEND PROPERTY COMPILE_DEFINITIONS ASPECT_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
ENDFOREACH()


# Configure a cmake fragment that plugins can use to
# set up compiler flags, include paths, etc to compile an
# Aspect plugin
CONFIGURE_FILE(
  ${CMAKE_SOURCE_DIR}/AspectConfig.cmake.in
  ${CMAKE_BINARY_DIR}/AspectConfig.cmake
  @ONLY
)


# Next, set up the testsuite
IF(EXISTS ${CMAKE_SOURCE_DIR}/tests/CMakeLists.txt)
  ENABLE_TESTING()
  ADD_SUBDIRECTORY(tests)
ENDIF()


MESSAGE("====================================================")

# Depending on whether we link statically or allow for shared libs,
# we can or can not load plugins via external shared libs. Pass this
# down during compilation so we can disable it in the code
IF (DEAL_II_STATIC_EXECUTABLE STREQUAL "ON")
  MESSAGE(STATUS "Creating a statically linked executable")
  MESSAGE(STATUS "Disabling dynamic loading of plugins from the input file")
  FOREACH(_source_file ${TARGET_SRC})
    SET_PROPERTY(SOURCE ${_source_file}
      APPEND PROPERTY COMPILE_DEFINITIONS ASPECT_USE_SHARED_LIBS=0)
  ENDFOREACH()
ELSE()
  MESSAGE(STATUS "Enabling dynamic loading of plugins from the input file")
  FOREACH(_source_file ${TARGET_SRC})
    SET_PROPERTY(SOURCE ${_source_file}
      APPEND PROPERTY COMPILE_DEFINITIONS ASPECT_USE_SHARED_LIBS=1)
  ENDFOREACH()
ENDIF()

DEAL_II_INVOKE_AUTOPILOT()

MESSAGE(STATUS "writing config into detailed.log...")
LIST(APPEND CMAKE_MODULE_PATH
  ${CMAKE_SOURCE_DIR}
  )

INCLUDE(write_config)

