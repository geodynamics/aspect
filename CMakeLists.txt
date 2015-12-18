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


MESSAGE(STATUS "====================================================")
MESSAGE(STATUS "============ Configuring ASPECT ====================")
MESSAGE(STATUS "====================================================")

# Set the name of the project and target:
SET(TARGET "aspect")

FILE(GLOB_RECURSE TARGET_SRC  "source/*.cc" "include/*.h")
INCLUDE_DIRECTORIES(include)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

FIND_PACKAGE(deal.II 8.2 QUIET
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n*** Could not find a suitably recent version of deal.II. ***\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake "
    "or set an environment variable \"DEAL_II_DIR\" that contains a path to a "
    "sufficiently recent version of deal.II."
    )
ENDIF()

MESSAGE(STATUS "Found deal.II version ${DEAL_II_PACKAGE_VERSION} at '${deal.II_DIR}'")

SET(ASPECT_USE_PETSC OFF CACHE BOOL "Use PETSc instead of Trilinos if set to 'on'.")

MESSAGE(STATUS "Using PETSc = '${ASPECT_USE_PETSC}'")

SET(_DEALII_GOOD ON)

IF(ASPECT_USE_PETSC AND NOT DEAL_II_WITH_PETSC)
    MESSAGE(SEND_ERROR
      "\n-- deal.II was built without support for PETSc!\n"
      )
    SET(_DEALII_GOOD OFF)
ENDIF()

IF(NOT DEAL_II_WITH_P4EST)
    MESSAGE(SEND_ERROR
      "\n-- deal.II was built without support for p4est!\n"
      )
    SET(_DEALII_GOOD OFF)
ENDIF()

IF(NOT ASPECT_USE_PETSC AND NOT DEAL_II_WITH_TRILINOS)
    MESSAGE(SEND_ERROR
      "\n-- deal.II was built without support for Trilinos!\n"
      )
    SET(_DEALII_GOOD OFF)
ENDIF()

IF (NOT _DEALII_GOOD)
  MESSAGE(FATAL_ERROR
    "\nAspect requires a deal.II installation built with support for Trilinos/PETSc and p4est but one or both of these appears to be missing!\n"
    )
ENDIF()

DEAL_II_INITIALIZE_CACHED_VARIABLES()
PROJECT(${TARGET})

IF (ASPECT_USE_PETSC)
FOREACH(_source_file ${TARGET_SRC})
  SET_PROPERTY(SOURCE ${_source_file}
    APPEND PROPERTY COMPILE_DEFINITIONS ASPECT_USE_PETSC="1")
ENDFOREACH()
ENDIF()


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
  ${CMAKE_SOURCE_DIR}/cmake/AspectConfig.cmake.in
  ${CMAKE_BINARY_DIR}/AspectConfig.cmake
  @ONLY
)


# Next, set up the testsuite

SET(ASPECT_RUN_ALL_TESTS OFF CACHE BOOL "Set up complete test suite to run.")

CONFIGURE_FILE(
  ${CMAKE_SOURCE_DIR}/tests/cmake/CTestCustom.ctest.in
  ${CMAKE_BINARY_DIR}/CTestCustom.ctest
  @ONLY
)
CONFIGURE_FILE(
  ${CMAKE_SOURCE_DIR}/tests/cmake/print_test_info.sh
  ${CMAKE_BINARY_DIR}/print_test_info.sh
  @ONLY
)

ADD_CUSTOM_TARGET(setup_tests
  COMMAND ${CMAKE_COMMAND} -D ASPECT_RUN_ALL_TESTS=ON -D PRINT_TEST_SUMMARY_AS_IMPORTANT=ON . >/dev/null
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Enabling all tests ...")

IF(EXISTS ${CMAKE_SOURCE_DIR}/tests/CMakeLists.txt)
  ENABLE_TESTING()
  ADD_SUBDIRECTORY(tests)
ENDIF()


ADD_CUSTOM_TARGET(generate_reference_output
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/cmake/generate_reference_output.sh 
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

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

# See whether we can verify that every plugin we load is compiled against
# the same deal.II library
INCLUDE (CheckIncludeFiles)
CHECK_INCLUDE_FILES ("link.h" _HAVE_LINK_H)
IF (_HAVE_LINK_H)
  MESSAGE(STATUS "Enabling checking of compatible deal.II library when loading plugins")
  FOREACH(_source_file ${TARGET_SRC})
    SET_PROPERTY(SOURCE ${_source_file}
      APPEND PROPERTY COMPILE_DEFINITIONS ASPECT_HAVE_LINK_H=1)
  ENDFOREACH()
ENDIF()

DEAL_II_INVOKE_AUTOPILOT()


# Check if we can raise floating point exceptions.
#
# Note that some library we link with in ASPECT on some platforms will trigger
# floating point exceptions when converting -numeric_limits<double>::max to a
# string. The only thing we can do is a configure check and disable the
# exceptions. This is done here:
SET(ASPECT_USE_FP_EXCEPTIONS ON CACHE BOOL "If ON, floating point exception are raised in debug mode.")

IF (ASPECT_USE_FP_EXCEPTIONS)
  INCLUDE(${CMAKE_SOURCE_DIR}/cmake/fpe_check.cmake)

  IF (HAVE_FP_EXCEPTIONS)
    MESSAGE(STATUS "Runtime floating point checks enabled.")
    FOREACH(_source_file ${TARGET_SRC})
      SET_PROPERTY(SOURCE ${_source_file}
        APPEND PROPERTY COMPILE_DEFINITIONS ASPECT_USE_FP_EXCEPTIONS=1)
    ENDFOREACH()
  ELSE()
    SET(ASPECT_USE_FP_EXCEPTIONS OFF CACHE BOOL "" FORCE)
    MESSAGE(STATUS "No support for feenableexcept(), disabling runtime floating point exception checks.")
  ENDIF()
ENDIF()



MESSAGE(STATUS "Writing config into detailed.log...")
LIST(APPEND CMAKE_MODULE_PATH
  ${CMAKE_SOURCE_DIR}
  )

INCLUDE(cmake/write_config)

