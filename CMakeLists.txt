# Copyright (C) 2013 - 2024 by the authors of the ASPECT code.
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
project(aspect CXX)

# load in version info and export it
file(STRINGS "${CMAKE_SOURCE_DIR}/VERSION" ASPECT_PACKAGE_VERSION LIMIT_COUNT 1)


message(STATUS "")
message(STATUS "====================================================")
message(STATUS "===== Configuring ASPECT ${ASPECT_PACKAGE_VERSION}")
message(STATUS "====================================================")
message(STATUS "")



# ####################################################################
# First figure out some global stuff
# ####################################################################


if(EXISTS ${CMAKE_SOURCE_DIR}/CMakeCache.txt)
  message(FATAL_ERROR  "Detected the file\n"
  "${CMAKE_SOURCE_DIR}/CMakeCache.txt\n"
  "in your source directory, which is a left-over from an in-source build. "
  "Please delete the file before running cmake from a separate build directory.")
endif()

if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  message(FATAL_ERROR  "ASPECT does not support in-source builds in ${CMAKE_BINARY_DIR}. "
  "Please run cmake from a separate build directory."
  "\n"
  "Note that CMake created a file called CMakeCache.txt and a folder called "
  "CMakeFiles in the source directory that you have to remove before you can "
  "begin a build in a different directory.")
endif()

# Set up CMAKE_BUILD_TYPE. Newer cmake versions (such as 3.29) set it
# to DebugRelease in the call to project() above, whereas older ones
# set it to the empty string or not at all. So if it's empty, set it
# to DebugRelease.
if ("${CMAKE_BUILD_TYPE}" STREQUAL "")
  set(CMAKE_BUILD_TYPE
      "DebugRelease"
      CACHE STRING
      "Choose the type of build, options are: Debug, Release and DebugRelease."
      FORCE
  )
endif()

# Then convert the awkward tri-state CMAKE_BUILD_TYPE into two bi-state
# variables ASPECT_BUILD_DEBUG and ASPECT_BUILD_RELEASE that are either
# ON or OFF.
if("${CMAKE_BUILD_TYPE}" STREQUAL "Release" OR
   "${CMAKE_BUILD_TYPE}" STREQUAL "Debug" OR
   "${CMAKE_BUILD_TYPE}" STREQUAL "DebugRelease" )
  message(STATUS "Setting up ASPECT for ${CMAKE_BUILD_TYPE} mode.")
else()
  message(FATAL_ERROR
    "CMAKE_BUILD_TYPE must either be 'Release', 'Debug', or 'DebugRelease', but is set to '${CMAKE_BUILD_TYPE}'.")
endif()

if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug" OR
   "${CMAKE_BUILD_TYPE}" STREQUAL "DebugRelease" )
  set(ASPECT_BUILD_DEBUG "ON")
else()
  set(ASPECT_BUILD_DEBUG "OFF")
endif()

if("${CMAKE_BUILD_TYPE}" STREQUAL "Release" OR
   "${CMAKE_BUILD_TYPE}" STREQUAL "DebugRelease" )
  set(ASPECT_BUILD_RELEASE "ON")
else()
  set(ASPECT_BUILD_RELEASE "OFF")
endif()

# Some other ASPECT configuration options:
set(ASPECT_MAX_NUM_PARTICLE_SYSTEMS 2 CACHE STRING "The maximum number of particle systems supported at runtime.")



# ##############################################################################
# Find external libraries and deal with external tools
# ##############################################################################

message(STATUS "")
message(STATUS "===== Configuring external libraries ===============")

list(APPEND CMAKE_MODULE_PATH
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/cmake/modules
  )


# deal.II first
find_package(deal.II 9.5.0 QUIET
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} $ENV{DEAL_II_DIR}
  )
if(NOT ${deal.II_FOUND})
  message(FATAL_ERROR "\n*** Could not find a suitably recent version of deal.II. ***\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake "
    "or set an environment variable \"DEAL_II_DIR\" that contains a path to a "
    "sufficiently recent version of deal.II."
    )
endif()

message(STATUS "Found deal.II version ${DEAL_II_PACKAGE_VERSION} at '${deal.II_DIR}'")

set(_DEALII_GOOD ON)

if(NOT DEAL_II_WITH_P4EST)
    message(SEND_ERROR
      "\n-- deal.II was built without support for p4est!\n"
      )
    set(_DEALII_GOOD OFF)
endif()

if(NOT DEAL_II_WITH_TRILINOS)
    message(SEND_ERROR
      "\n-- deal.II was built without support for Trilinos!\n"
      )
    set(_DEALII_GOOD OFF)
endif()

if(NOT DEAL_II_WITH_SUNDIALS)
    message(SEND_ERROR
      "\n-- deal.II was built without support for SUNDIALS!\n"
      )
    set(_DEALII_GOOD OFF)
endif()

if(NOT _DEALII_GOOD)
  message(FATAL_ERROR
    "\nASPECT requires a deal.II installation built with certain features enabled that seem to be missing (see above)!\n"
    )
endif()

deal_ii_initialize_cached_variables()


# zlib
set(ZLIB_DIR "" CACHE PATH "An optional hint to a ZLIB installation")
if("${ZLIB_DIR}" STREQUAL "")
  set(ZLIB_DIR "$ENV{ZLIB_DIR}" CACHE PATH "An optional hint to a ZLIB installation" FORCE)
endif()
find_package(ZLIB)
if(${ZLIB_FOUND})
  message(STATUS "Using ASPECT_WITH_ZLIB = 'ON'")
  set(ASPECT_WITH_ZLIB ON)
else()
  message(STATUS "Using ASPECT_WITH_ZLIB = 'OFF'")
  set(ASPECT_WITH_ZLIB OFF)
endif()


# PerpleX
find_package(PerpleX QUIET
  HINTS ./contrib/perplex/install/ ../ ../../ ${PERPLEX_DIR} $ENV{PERPLEX_DIR})
if(${PerpleX_FOUND})
  message(STATUS "PerpleX found at ${PerpleX_INCLUDE_DIR}")
  set(ASPECT_WITH_PERPLEX ON)
else()
  set(ASPECT_WITH_PERPLEX OFF)
endif()


# libdap so that ASPECT can connect to the OPeNDAP servers
set(ASPECT_WITH_LIBDAP OFF CACHE BOOL "Check if the user wants to compile ASPECT with the libdap libraries.")
message(STATUS "Using ASPECT_WITH_LIBDAP = '${ASPECT_WITH_LIBDAP}'")
if(ASPECT_WITH_LIBDAP)
  find_package(LIBDAP)
  if(${LIBDAP_FOUND})
    include_directories(${LIBDAP_INCLUDE_DIRS})
    message(STATUS "LIBDAP found at ${LIBDAP_LIBRARY}\n")
  else()
    message(FATAL_ERROR "LIBDAP not found. Disable ASPECT_WITH_LIBDAP or specify a hint to your installation directory with LIBDAP_DIR.")
  endif()
endif()


# FastScape library and link it to ASPECT if requested
set(ASPECT_WITH_FASTSCAPE OFF CACHE BOOL "Whether the user wants to compile ASPECT with the landscape evolution code FastScape, or not.")
message(STATUS "Using ASPECT_WITH_FASTSCAPE = '${ASPECT_WITH_FASTSCAPE}'")
if(ASPECT_WITH_FASTSCAPE)
  find_library(FASTSCAPE NAMES fastscapelib_fortran PATHS $ENV{FASTSCAPE_DIR} ${FASTSCAPE_DIR} PATH_SUFFIXES lib NO_DEFAULT_PATH)
  if(FASTSCAPE)
     message(STATUS "FastScape library found at ${FASTSCAPE_DIR}")

    # Get the fastscape source path so we can check the version.
    file(STRINGS "${FASTSCAPE_DIR}/Makefile" _fastscape_makefile)
    foreach(_line ${_fastscape_makefile})
        if("${_line}" MATCHES "^CMAKE_SOURCE_DIR")
          string(REPLACE "CMAKE_SOURCE_DIR = " "" FASTSCAPE_SOURCE_DIR ${_line})
        endif()
    endforeach()

    # Now get the version from setup.py
    message(STATUS "Parsing '${FASTSCAPE_SOURCE_DIR}/setup.py' for version information")
    file(STRINGS "${FASTSCAPE_SOURCE_DIR}/setup.py" _fastscape_info)
    foreach(_line ${_fastscape_info})
      string(STRIP ${_line} _line)
      if("${_line}" MATCHES "^version")
        string(REGEX REPLACE "[^0-9.]" "" FASTSCAPE_VERSION ${_line})
      endif()
    endforeach()

    # Throw an error if the version is too old.
    if(${FASTSCAPE_VERSION} VERSION_LESS 2.8)
      message(FATAL_ERROR "The linked FastScape version is ${FASTSCAPE_VERSION}, however at least 2.8.0 is required.")
    endif()
  else()
     message(FATAL_ERROR "Trying to link with FastScape but libfastscapelib_fortran.so was not found in ${FASTSCAPE_DIR}")
  endif()
endif()


# NetCDF (c including parallel)
set(ASPECT_WITH_NETCDF ON CACHE BOOL "Check if the user wants to compile ASPECT with the NetCDF libraries.")

if(ASPECT_WITH_NETCDF)
  find_package(NETCDF)
  if(${NETCDF_FOUND})
    message(STATUS "Using ASPECT_WITH_NETCDF = '${ASPECT_WITH_NETCDF}'")
    message(STATUS "  NETCDF_INCLUDE_DIR: ${NETCDF_INCLUDE_DIR}")
    message(STATUS "  NETCDF_LIBRARY: ${NETCDF_LIBRARY}")
    message(STATUS "  NETCDF_VERSION: ${NETCDF_VERSION}")
  else()
    message(STATUS "NetCDF not found. Disabling ASPECT_WITH_NETCDF. You can specify a hint to your installation directory with NETCDF_DIR.")
    set(ASPECT_WITH_NETCDF OFF CACHE BOOL "" FORCE)
  endif()
else()
  message(STATUS "Using ASPECT_WITH_NETCDF = 'OFF'")
endif()


# WorldBuilder
set(ASPECT_WITH_WORLD_BUILDER ON CACHE BOOL "Whether to enable compiling aspect with the Geodynamic World Builder.")
message(STATUS "Using ASPECT_WITH_WORLD_BUILDER = '${ASPECT_WITH_WORLD_BUILDER}'")
if(ASPECT_WITH_WORLD_BUILDER)
  # Check whether we can find the WorldBuilder.
  set(WORLD_BUILDER_SOURCE_DIR "" CACHE PATH "Provide an external World Builder directory to be compiled with ASPECT. If the path is not provided or the World Builder is not found in the provided location, the version in the contrib folder is used.")
  if(NOT EXISTS ${WORLD_BUILDER_SOURCE_DIR}/VERSION)
    message(STATUS "World Builder not found. Using internal version.")
    set(WORLD_BUILDER_SOURCE_DIR "${CMAKE_SOURCE_DIR}/contrib/world_builder/" CACHE PATH "" FORCE)
  endif()

  # Add the WorldBuilder include dirs. This includes both the WorldBuilder
  # header files, as well as where the cmake run we do below puts the
  # WorldBuilder config.h file.
  include_directories("${WORLD_BUILDER_SOURCE_DIR}/include/")
  include_directories("${CMAKE_BINARY_DIR}/world_builder/include")


  include("${WORLD_BUILDER_SOURCE_DIR}/cmake/version.cmake")

  message(STATUS "Using World Builder version ${WORLD_BUILDER_VERSION} found at ${WORLD_BUILDER_SOURCE_DIR}.")

  if(WORLD_BUILDER_VERSION VERSION_GREATER_EQUAL 0.6.0)
    # if we configured with 0.5.0 before, we have a stray include file here. delete it.
    file(REMOVE "${CMAKE_BINARY_DIR}/include/world_builder/config.h")

    set(WB_ENABLE_TESTS OFF)
    set(WB_ENABLE_APPS OFF)
    set(WB_ENABLE_HELPER_TARGETS OFF)
    set(WB_ENABLE_PYTHON OFF)
    set(WB_MAKE_FORTRAN_WRAPPER OFF)

    # First for debug mode...
    if(${ASPECT_BUILD_DEBUG} STREQUAL "ON")
      set(_build_type ${CMAKE_BUILD_TYPE})
      set(CMAKE_BUILD_TYPE Debug)

      message(STATUS "")
      message(STATUS "===== Setting up Geodynamic World Builder version ${WORLD_BUILDER_VERSION} in ${CMAKE_BUILD_TYPE} mode")

      set(WB_TARGET "WorldBuilderDebug")
      add_subdirectory("${WORLD_BUILDER_SOURCE_DIR}" ${CMAKE_BINARY_DIR}/world_builder/)
      separate_arguments(SEPARATED_FLAGS NATIVE_COMMAND ${ASPECT_ADDITIONAL_CXX_FLAGS})
      target_compile_options(WorldBuilderDebug PRIVATE "-g" "${SEPARATED_FLAGS}")

      set(CMAKE_BUILD_TYPE ${_build_type})
    endif()

    # ...then for release mode
    if(${ASPECT_BUILD_RELEASE} STREQUAL "ON")
      set(_build_type ${CMAKE_BUILD_TYPE})
      set(CMAKE_BUILD_TYPE Release)

      message(STATUS "")
      message(STATUS "===== Setting up Geodynamic World Builder version ${WORLD_BUILDER_VERSION} in ${CMAKE_BUILD_TYPE} mode")

      set(WB_TARGET "WorldBuilderRelease")
      add_subdirectory("${WORLD_BUILDER_SOURCE_DIR}" ${CMAKE_BINARY_DIR}/world_builder_release/)
      target_compile_definitions(WorldBuilderRelease PUBLIC "NDEBUG")

      separate_arguments(SEPARATED_FLAGS NATIVE_COMMAND ${ASPECT_ADDITIONAL_CXX_FLAGS})
      target_compile_options(WorldBuilderRelease PRIVATE "-O3" "${SEPARATED_FLAGS}")

      set(CMAKE_BUILD_TYPE ${_build_type})
    endif()

  else()

    file(GLOB_RECURSE wb_files "${WORLD_BUILDER_SOURCE_DIR}/source/world_builder/*.cc")
    list(APPEND ASPECT_SOURCE_FILES ${wb_files})
    add_definitions(-DWB_WITH_MPI)

    # generate config.cc and include it:
    configure_file("${WORLD_BUILDER_SOURCE_DIR}/include/world_builder/config.h.in"
                    "${CMAKE_BINARY_DIR}/include/world_builder/config.h" @ONLY)
    include_directories("${CMAKE_BINARY_DIR}/include/")
    file(REMOVE "${CMAKE_BINARY_DIR}/source/world_builder/config.cc")

    # Move some file to the end for unity builds to make sure other file come
    # "before". Note: The current design keeps all ASPECT files before all GWB
    # files. Mixing them will causes many issues with non-unique namespace
    # names like Utilities.
    set(UNITY_WB_LAST_FILES
        "${WORLD_BUILDER_SOURCE_DIR}/source/world_builder/parameters.cc")

    foreach(_source_file ${UNITY_WB_LAST_FILES})
      list(FIND ASPECT_SOURCE_FILES ${_source_file} _index)
      if(_index EQUAL -1)
        message(FATAL_ERROR "could not find ${_source_file}.")
      endif()

      list(REMOVE_ITEM ASPECT_SOURCE_FILES ${_source_file})
      list(APPEND ASPECT_SOURCE_FILES ${_source_file})
    endforeach()

    foreach(_source_file ${wb_files})
      # exclude the world builder files from including precompiled headers, they
      # do not include ASPECT's dependencies at all.
      set_property(SOURCE ${_source_file} PROPERTY SKIP_PRECOMPILE_HEADERS TRUE )
      # Temporarily disable world builder unity builds:
      set_property(SOURCE ${_source_file} PROPERTY SKIP_UNITY_BUILD_INCLUSION TRUE )
    endforeach()

    set_property(SOURCE "${CMAKE_BINARY_DIR}/world_builder_config.cc" PROPERTY SKIP_PRECOMPILE_HEADERS TRUE )
  endif()
endif()
message(STATUS "")


#
# Other stuff about external tools and how we interact with the system:
#

# Depending on whether we link statically or allow for shared libs,
# we can or can not load plugins via external shared libs. Pass this
# down during compilation so we can disable it in the code.
set(ASPECT_USE_SHARED_LIBS ON CACHE BOOL "If ON, we support loading shared plugin files.")
if(DEAL_II_STATIC_EXECUTABLE STREQUAL "ON")
  message(STATUS "Creating a statically linked executable")
  set(ASPECT_USE_SHARED_LIBS OFF CACHE BOOL "" FORCE)
endif()


# Generate compile_commands.json for tooling (VS Code, etc.)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(FORCE_COLORED_OUTPUT ON
    CACHE BOOL "Forces colored output when compiling with gcc and clang.")

# Query the current git commit so we can put it into the start-up
# message of ASPECT.
include(${CMAKE_SOURCE_DIR}/cmake/macro_aspect_query_git_information.cmake)
ASPECT_QUERY_GIT_INFORMATION("ASPECT")
configure_file(${CMAKE_SOURCE_DIR}/include/aspect/revision.h.in
               ${CMAKE_BINARY_DIR}/include/aspect/revision.h @ONLY)


# Check if we can raise floating point exceptions.
#
# Note that some library we link with in ASPECT on some platforms will trigger
# floating point exceptions when converting -numeric_limits<double>::max to a
# string. The only thing we can do is a configure check and disable the
# exceptions. This is done here:
set(ASPECT_USE_FP_EXCEPTIONS ON CACHE BOOL "If ON, floating point exception are raised in debug mode.")

# Clang 6.0 throws random floating point exceptions, which we could not
# track down. Disable the exceptions for now.
if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" AND NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 6.0.0)
  set(ASPECT_USE_FP_EXCEPTIONS OFF CACHE BOOL "" FORCE)
endif()

if(ASPECT_USE_FP_EXCEPTIONS)
  message(STATUS "Determining whether we can use floating point exceptions...")
  include(${CMAKE_SOURCE_DIR}/cmake/fpe_check.cmake)

  if(HAVE_FP_EXCEPTIONS)
    message(STATUS "Runtime floating point checks enabled.")
  else()
    set(ASPECT_USE_FP_EXCEPTIONS OFF CACHE BOOL "" FORCE)
    message(STATUS "No support for feenableexcept(), disabling runtime floating point exception checks.")
  endif()
endif()

# Check if we want to precompile header files. This speeds up compile time,
# but can fail on some machines with old CMake. Starting with CMake 3.16
# there is native support inside CMake and we can precompile headers.
#
# While there, also set up the variable to use "unity" builds.
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16)
  set(ASPECT_PRECOMPILE_HEADERS ON CACHE BOOL "Precompile external header files to speedup compile time. Currently only supported for CMake 3.16 and newer versions.")
  set(ASPECT_UNITY_BUILD ON CACHE BOOL "Combine source files into less compile targets to speedup compile time. Currently only supported for CMake 3.16 and newer versions.")
else()
  set(ASPECT_PRECOMPILE_HEADERS OFF CACHE BOOL "Precompile external header files to speedup compile time. Currently only supported for CMake 3.16 and newer versions." FORCE)
  set(ASPECT_UNITY_BUILD OFF CACHE BOOL "Combine source files into less compile targets to speedup compile time. Currently only supported for CMake 3.16 and newer versions." FORCE)
endif()

if(ASPECT_PRECOMPILE_HEADERS AND CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  # Intel 19.1 produces internal compiler errors inside bundled boost with
  # precompiled headers, so we deactivate it:
  set(ASPECT_PRECOMPILE_HEADERS OFF CACHE BOOL "" FORCE)
endif()




# ##############################################################################
# Configuring the ASPECT targets themselves
# ##############################################################################

message(STATUS "")
message(STATUS "===== Configuring ASPECT build targets =============")


# ##############################################################################
# Generate config.h

# First, define macro with the source directory that will be exported in
# config.h. This can be used to hard-code the location of data files, such as
# in $ASPECT_SOURCE_DIR/data/velocity-boundary-conditions/gplates/*
set(ASPECT_SOURCE_DIR ${CMAKE_SOURCE_DIR})

# And finally generate the file
configure_file(
  ${CMAKE_SOURCE_DIR}/include/aspect/config.h.in
  ${CMAKE_BINARY_DIR}/include/aspect/config.h
  )


# ##############################################################################
# Collect source files

file(GLOB_RECURSE ASPECT_SOURCE_FILES    "source/*.cc" "include/*.h")
list(REMOVE_ITEM  ASPECT_SOURCE_FILES    ${CMAKE_CURRENT_SOURCE_DIR}/source/main.cc)
file(GLOB_RECURSE ASPECT_UNIT_TEST_FILES "unit_tests/*.cc" "contrib/catch/catch.hpp")
set(ASPECT_MAIN_FILE "source/main.cc")

# Special treatment of some files for unity builds. We need to exclude them
# because otherwise template class instantiations split across several
# .cc files (for example Simulator<dim> in core.cc, helper_functions.cc,
# solver_schemes.cc, etc. might end up in the same unity file. This would
# lead to errors with duplicate instantiations. The only reliable option
# is to exclude the files that instantiate the whole class.
set(UNITY_SEPARATE_FILES
  "source/simulator/core.cc;source/volume_of_fluid/handler.cc")

foreach(_source_file ${UNITY_SEPARATE_FILES})
  set(_full_name "${CMAKE_SOURCE_DIR}/${_source_file}")
  set_property(SOURCE ${_full_name} PROPERTY SKIP_UNITY_BUILD_INCLUSION TRUE)
endforeach()



# Set up include directories. Put the ASPECT header files
# to the front of the list, whereas the 'catch' headers can
# be preempted by the system
include_directories(BEFORE ${CMAKE_BINARY_DIR}/include include)
include_directories(AFTER  contrib/catch)

# Set the name of the main targets in the form of
# TARGET_EXE/LIB_DEBUG/RELEASE. Set TARGET_EXE as the debug build,
# unless we have only release mode enabled.
if(${ASPECT_BUILD_DEBUG} STREQUAL "ON")
  set(TARGET_EXE_DEBUG "aspect.exe.debug")
  list(APPEND TARGET_EXECUTABLES "${TARGET_EXE_DEBUG}")

  set(TARGET_EXE "${TARGET_EXE_DEBUG}")
endif()

if(${ASPECT_BUILD_RELEASE} STREQUAL "ON")
  set(TARGET_EXE_RELEASE "aspect.exe.release")
  list(APPEND TARGET_EXECUTABLES "${TARGET_EXE_RELEASE}")

  if ("${TARGET_EXE}" STREQUAL "")
    set(TARGET_EXE "${TARGET_EXE_RELEASE}")
  endif()
endif()


message(STATUS "Writing externally readable configuration information")
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  "${CMAKE_BINARY_DIR}/AspectConfigVersion.cmake"
  VERSION ${ASPECT_PACKAGE_VERSION}
  COMPATIBILITY AnyNewerVersion
)

# Configure a cmake fragment that plugins can use to
# set up compiler flags, include paths, etc to compile an
# ASPECT plugin.
# Config for the build dir:
set(CONFIG_INCLUDE_DIRS "${CMAKE_BINARY_DIR}/include" "${CMAKE_SOURCE_DIR}/include")
set(CONFIG_DIR "${CMAKE_BINARY_DIR}")
configure_file(
  ${CMAKE_SOURCE_DIR}/cmake/AspectConfig.cmake.in
  ${CMAKE_BINARY_DIR}/AspectConfig.cmake
  @ONLY
)
# Config for the install dir:
set(CONFIG_INCLUDE_DIRS "${CMAKE_INSTALL_PREFIX}/include")
set(CONFIG_DIR "${CMAKE_INSTALL_PREFIX}/bin")
configure_file(
  ${CMAKE_SOURCE_DIR}/cmake/AspectConfig.cmake.in
  ${CMAKE_BINARY_DIR}/forinstall/AspectConfig.cmake
  @ONLY
)

# Next, set up the testsuite

set(ASPECT_RUN_ALL_TESTS OFF CACHE BOOL "Set up complete test suite to run.")
set(ASPECT_NEED_TEST_CONFIGURE ON CACHE BOOL "If true, reconfigure test project.")
set(ASPECT_COMPARE_TEST_RESULTS ON CACHE BOOL "Compare test results with high accuracy.")

configure_file(
  ${CMAKE_SOURCE_DIR}/tests/cmake/CTestCustom.ctest.in
  ${CMAKE_BINARY_DIR}/CTestCustom.ctest
  @ONLY
)
configure_file(
  ${CMAKE_SOURCE_DIR}/tests/cmake/print_test_info.sh
  ${CMAKE_BINARY_DIR}/print_test_info.sh
  @ONLY
)

# enable all tests and force a re-run of the cmake configuration in the test/ folder:
add_custom_target(setup_tests
  COMMAND ${CMAKE_COMMAND} -D ASPECT_RUN_ALL_TESTS=ON -D ASPECT_NEED_TEST_CONFIGURE=ON . >/dev/null
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Enabling all tests ...")

set(ASPECT_TEST_GENERATOR "Unix Makefiles" CACHE STRING
  "Generator to use for the test cmake project. Using ninja instead of make is not recommended.")


# This is an empty target but we will make it depend on tests/ and unit_tests/ next:
add_custom_target(test)

# Generate CTestTestfile.cmake in the main build folder that lists all subfolders
# that contain tests. This way you can call "ctest" in the build directory.
file(WRITE ${CMAKE_BINARY_DIR}/CTestTestfile.cmake "# auto-generated ctest file\n")

if(EXISTS ${CMAKE_SOURCE_DIR}/unit_tests/CMakeLists.txt)
  # If we have the unit_tests directory, add it so ctest picks it up and configure
  # the test project in the subfolder:
  file(APPEND ${CMAKE_BINARY_DIR}/CTestTestfile.cmake "SUBDIRS(unit_tests)\n")
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/unit_tests)
  execute_process(
    COMMAND ${CMAKE_COMMAND}
        -D CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -D ASPECT_BINARY=${CMAKE_BINARY_DIR}/aspect
        ${CMAKE_CURRENT_SOURCE_DIR}/unit_tests
    OUTPUT_FILE setup_unit_tests.log
    RESULT_VARIABLE test_cmake_result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/unit_tests
   )
  if(NOT test_cmake_result EQUAL 0)
    message(FATAL_ERROR "ERROR: unittest/ project could not be configured.")
  endif()

  # Finally create a custom target to run the tests:
  add_custom_target(run_unit_tests
    COMMAND ${CMAKE_BINARY_DIR}/aspect --test
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/unit_tests
    DEPENDS ${TARGET_EXE}
    COMMENT "Running unit_tests ...")
endif()

if(EXISTS ${CMAKE_SOURCE_DIR}/tests/CMakeLists.txt)
  # Hook up the tests:
  file(APPEND ${CMAKE_BINARY_DIR}/CTestTestfile.cmake "SUBDIRS(tests)\n")
endif()

if(EXISTS ${CMAKE_SOURCE_DIR}/tests/CMakeLists.txt
  AND ${ASPECT_NEED_TEST_CONFIGURE})

  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/tests)

  # set a flag so we don't need to rerun this configuration step every time:
  set(ASPECT_NEED_TEST_CONFIGURE OFF CACHE BOOL "" FORCE)

  message(STATUS "Setting up test project, see tests/setup_tests.log for details.")
  execute_process(
    COMMAND ${CMAKE_COMMAND} -G ${ASPECT_TEST_GENERATOR}
        -D ASPECT_RUN_ALL_TESTS=${ASPECT_RUN_ALL_TESTS}
        -D ASPECT_COMPARE_TEST_RESULTS=${ASPECT_COMPARE_TEST_RESULTS}
        -D Aspect_DIR=${CMAKE_BINARY_DIR}
        -D CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -D ASPECT_BINARY=${CMAKE_BINARY_DIR}/aspect
        ${CMAKE_CURRENT_SOURCE_DIR}/tests
    OUTPUT_FILE setup_tests.log
    RESULT_VARIABLE test_cmake_result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/tests
   )
  if(NOT test_cmake_result EQUAL 0)
    message(FATAL_ERROR "ERROR: tests/ project could not be configured.")
  endif()

  # Finally create a custom target:
  add_custom_target(run_tests
    COMMAND ${CMAKE_COMMAND} --build . --target test
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/tests
    DEPENDS ${TARGET_EXE}
    COMMENT "Running tests ...")
endif()

# Provide the "generate_reference_output" target:
add_custom_target(generate_reference_output
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/cmake/generate_reference_output.sh
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})


# Provide "indent" target for indenting all headers and source files
add_custom_target(indent
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  COMMAND ./contrib/utilities/indent
  COMMENT "Indenting all ASPECT header and source files..."
  )

if(CMAKE_GENERATOR MATCHES "Ninja")
  set(_make_command "$ ninja")
else()
  set(_make_command " $ make")
endif()

# Provide "release" and "debug" targets to switch compile mode
if(${DEAL_II_BUILD_TYPE} MATCHES "DebugRelease")
add_custom_target(release
  COMMAND ${CMAKE_COMMAND} -D CMAKE_BUILD_TYPE=Release .
  COMMAND ${CMAKE_COMMAND} -E echo "***"
  COMMAND ${CMAKE_COMMAND} -E echo "*** Switched to Release mode. Now recompile with: ${_make_command}"
  COMMAND ${CMAKE_COMMAND} -E echo "***"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  VERBATIM
  COMMENT "switching to RELEASE mode..."
  )
add_custom_target(debug
  COMMAND ${CMAKE_COMMAND} -D CMAKE_BUILD_TYPE=Debug .
  COMMAND ${CMAKE_COMMAND} -E echo "***"
  COMMAND ${CMAKE_COMMAND} -E echo "*** Switched to Debug mode. Now recompile with: ${_make_command}"
  COMMAND ${CMAKE_COMMAND} -E echo "***"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  VERBATIM
  COMMENT "switching to DEBUG mode..."
  )
add_custom_target(debugrelease
  COMMAND ${CMAKE_COMMAND} -D CMAKE_BUILD_TYPE=DebugRelease .
  COMMAND ${CMAKE_COMMAND} -E echo "***"
  COMMAND ${CMAKE_COMMAND} -E echo "*** Switched to Debug and Release mode. Now recompile with: ${_make_command}"
  COMMAND ${CMAKE_COMMAND} -E echo "***"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  VERBATIM
  COMMENT "switching to DEBUG/RELEASE mode..."
  )
endif()

# Provide a "distclean" target (like it is done in deal.II):
add_custom_target(distclean
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target clean
  COMMAND ${CMAKE_COMMAND} -E remove_directory CMakeFiles
  COMMAND ${CMAKE_COMMAND} -E remove
    CMakeCache.txt cmake_install.cmake Makefile
    build.ninja rules.ninja .ninja_deps .ninja_log
  COMMENT "distclean invoked"
  )



file(WRITE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/print_usage.cmake
"message(
\"###
#
#  Project ${TARGET_EXE} set up with  ${DEAL_II_PACKAGE_NAME}-${DEAL_II_PACKAGE_VERSION}  found at
#      ${DEAL_II_PATH}
#
#  CMAKE_BUILD_TYPE:          ${CMAKE_BUILD_TYPE}
#
#  You can now run
#      ${_make_command}                - to compile and link ${TARGET_EXE}
#      ${_make_command} debug          - to switch the build type to 'Debug'
#      ${_make_command} release        - to switch the build type to 'Release'
#      ${_make_command} debugrelease   - to switch the build type to compile both
#      ${_make_command} clean          - to remove the generated executable as well as
#                               all intermediate compilation files
#      ${_make_command} distclean      - to clean the directory from all generated
#                               files (includes clean, runclean and the removal
#                               of the generated build system)
#      ${_make_command} setup_tests    - enable all tests and re-run test detection
#      ${_make_command} indent         - fix indentation of all source files
#      ${_make_command} info           - to view this message again
\")")

# Provide "info" target
add_custom_target(info
  COMMAND ${CMAKE_COMMAND} -P ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/print_usage.cmake
  )


INCLUDE (CheckCXXSourceCompiles)

set(_backup_libs ${CMAKE_REQUIRED_LIBRARIES})
list(APPEND CMAKE_REQUIRED_LIBRARIES ${CMAKE_DL_LIBS})
check_cxx_source_compiles("
#include <cstddef>
#include <dlfcn.h>

int main()
{
  void *handle = dlopen (\"somelib.so\", RTLD_LAZY);
  return handle == NULL || dlerror();
}
" HAVE_DLOPEN)
set(CMAKE_REQUIRED_LIBRARIES ${_backup_libs})

if(NOT HAVE_DLOPEN)
  message(STATUS "dlopen() test failed, disabling dynamic plugin loading")
  set(ASPECT_USE_SHARED_LIBS OFF CACHE BOOL "" FORCE)
endif()

if(ASPECT_USE_SHARED_LIBS)
  message(STATUS "Enabling dynamic loading of plugins from the input file")
else()
  message(STATUS "Disabling dynamic loading of plugins from the input file")
endif()

# See whether we can verify that every plugin we load is compiled against
# the same deal.II library
set(ASPECT_HAVE_LINK_H ON CACHE BOOL "If ON, link.h exists and is usable.")
include(CheckIncludeFileCXX)
check_include_file_cxx("link.h" _HAVE_LINK_H)
if(NOT _HAVE_LINK_H)
  set(ASPECT_HAVE_LINK_H OFF CACHE BOOL "" FORCE)
endif()
if(ASPECT_HAVE_LINK_H)
  message(STATUS "Enabling checking of compatible deal.II library when loading plugins")
endif()

if (${FORCE_COLORED_OUTPUT})
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER MATCHES "AppleClang")
    string(APPEND DEAL_II_CXX_FLAGS_DEBUG " -fcolor-diagnostics")
    string(APPEND DEAL_II_CXX_FLAGS_RELEASE " -fcolor-diagnostics")
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    string(APPEND DEAL_II_CXX_FLAGS_DEBUG " -fdiagnostics-color=always")
    string(APPEND DEAL_II_CXX_FLAGS_RELEASE " -fdiagnostics-color=always")
  endif()
endif()

# deal.II versions >=9.5 disable deprecation warnings in user code. Enable
# the warnings again by removing the flag that disables them.
string(REPLACE "-Wno-deprecated-declarations" "" DEAL_II_WARNING_FLAGS "${DEAL_II_WARNING_FLAGS}")

set(ASPECT_ADDITIONAL_CXX_FLAGS "" CACHE STRING "Additional CMAKE_CXX_FLAGS applied after the deal.II options.")

if(NOT ASPECT_ADDITIONAL_CXX_FLAGS STREQUAL "")
  message(STATUS "Appending ASPECT_ADDITIONAL_CXX_FLAGS: '${ASPECT_ADDITIONAL_CXX_FLAGS}':")
  string(APPEND DEAL_II_CXX_FLAGS_DEBUG " ${ASPECT_ADDITIONAL_CXX_FLAGS}")
  string(APPEND DEAL_II_CXX_FLAGS_RELEASE " ${ASPECT_ADDITIONAL_CXX_FLAGS}")
  message(STATUS "  DEAL_II_WARNING_FLAGS: ${DEAL_II_WARNING_FLAGS}")
  message(STATUS "  DEAL_II_CXX_FLAGS_DEBUG: ${DEAL_II_CXX_FLAGS_DEBUG}")
  message(STATUS "  DEAL_II_CXX_FLAGS_RELEASE: ${DEAL_II_CXX_FLAGS_RELEASE}")
endif()


# Create the combined set of source files that make up the project
set(TARGET_SRC
    ${ASPECT_MAIN_FILE}
    ${ASPECT_UNIT_TEST_FILES}
    ${ASPECT_SOURCE_FILES})



# Set up targets for ASPECT executables:
if(${ASPECT_BUILD_DEBUG} STREQUAL "ON")
  add_executable(${TARGET_EXE_DEBUG} ${TARGET_SRC})
  set_property(TARGET ${TARGET_EXE_DEBUG}
               PROPERTY OUTPUT_NAME aspect-debug)
  deal_ii_setup_target(${TARGET_EXE_DEBUG} DEBUG)
endif()

if(${ASPECT_BUILD_RELEASE} STREQUAL "ON")
  add_executable(${TARGET_EXE_RELEASE} ${TARGET_SRC})
  set_property(TARGET ${TARGET_EXE_RELEASE}
               PROPERTY OUTPUT_NAME aspect-release)
  deal_ii_setup_target(${TARGET_EXE_RELEASE} RELEASE)
endif()

# Make sure that we have an executable 'aspect' that
# links to the debug version, or the release version if
# we did not compile in debug mode. TARGET_EXE describes
# which target this references.
#
# We do this in the build directory with the first of the following
# two commands, and in the install directory with the second.
add_custom_target(aspect ALL
                  COMMAND ${CMAKE_COMMAND} -E create_symlink
                           $<TARGET_PROPERTY:${TARGET_EXE},OUTPUT_NAME>  # executable file name
                           ${CMAKE_BINARY_DIR}/aspect                    # link name
                  DEPENDS ${TARGET_EXE}
                  COMMENT "Creating symlink from ${TARGET_EXE} to './aspect'")
install(PROGRAMS ${CMAKE_BINARY_DIR}/aspect
	TYPE BIN)

# This should be done in DEAL_II_SETUP_TARGET, but deal.II 9.5
# does not do so, so do it here. Remove this when we require deal.II 9.6.
if(${ASPECT_BUILD_RELEASE} STREQUAL "ON")
  if(DEAL_II_PACKAGE_VERSION VERSION_LESS 9.6.0)
    target_compile_definitions(${TARGET_EXE_RELEASE} PUBLIC "NDEBUG")
  endif()
endif()

# turn on debug checks like invalid element access in c++ standard library:
if(${ASPECT_BUILD_DEBUG} STREQUAL "ON")
  target_compile_definitions(${TARGET_EXE_DEBUG} PUBLIC "_GLIBCXX_ASSERTIONS")
endif()


# ##############################################################################
# Make sure we correctly link with external libraries
# ##############################################################################

# zlib
if(${ZLIB_FOUND})
  message(STATUS "Linking ASPECT against zlib")
  include_directories(${ZLIB_INCLUDE_DIR})
  foreach(_T ${TARGET_EXECUTABLES})
    target_link_libraries(${_T} ${ZLIB_LIBRARY})
  endforeach()
endif()

# PerpleX
if(${PerpleX_FOUND})
  message(STATUS "Linking ASPECT against PerpleX")
  include_directories(${PerpleX_INCLUDE_DIR})
  foreach(_T ${TARGET_EXECUTABLES})
    target_link_libraries(${_T} ${PerpleX_LIBRARIES})
  endforeach()
endif()

# libdap
if(${LIBDAP_FOUND})
  message(STATUS "Linking ASPECT against libdap")
  include_directories(${LIBDAP_INCLUDE_DIRS})
  foreach(_T ${TARGET_EXECUTABLES})
    target_link_libraries(${_T} ${LIBDAP_LIBRARIES})
  endforeach()
endif()

# Fastscape
if (FASTSCAPE)
  message(STATUS "Linking ASPECT against Fastscape")
  foreach(_T ${TARGET_EXECUTABLES})
    target_link_libraries(${_T} ${FASTSCAPE})
  endforeach()
endif()

# NetCDF
if(${NETCDF_FOUND})
  message(STATUS "Linking ASPECT against NetCDF")
  include_directories(${NETCDF_INCLUDE_DIRS})
  foreach(_T ${TARGET_EXECUTABLES})
    target_link_libraries(${_T} ${NETCDF_LIBRARIES})
  endforeach()
endif()

# WorldBuilder
if(WORLD_BUILDER_VERSION VERSION_GREATER_EQUAL 0.6.0 AND ASPECT_WITH_WORLD_BUILDER)
  message(STATUS "Linking ASPECT against WorldBuilder")
  if(${ASPECT_BUILD_DEBUG} STREQUAL "ON")
    target_include_directories(${TARGET_EXE_DEBUG} PUBLIC "${CMAKE_BINARY_DIR}/world_builder/include/")
    target_link_libraries(${TARGET_EXE_DEBUG} WorldBuilderDebug)
  endif()

  if(${ASPECT_BUILD_RELEASE} STREQUAL "ON")
    target_include_directories(${TARGET_EXE_RELEASE} PUBLIC "${CMAKE_BINARY_DIR}/world_builder_release/include/")
    target_link_libraries(${TARGET_EXE_RELEASE} WorldBuilderRelease)
  endif()
endif()

# Some systems need to explicitly link to some libraries to use dlopen
if(ASPECT_USE_SHARED_LIBS)
  message(STATUS "Linking ASPECT against dlopen")
  foreach(_T ${TARGET_EXECUTABLES})
    target_link_libraries(${_T} ${CMAKE_DL_LIBS})
  endforeach()
endif()



if(ASPECT_PRECOMPILE_HEADERS)
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16)
    message(STATUS "Precompiling common header files.")
    # Use the native cmake support to precompile some common headers
    # from ASPECT and deal.II that are frequently included, but rarely changed.
    foreach(_T ${TARGET_EXECUTABLES})
      TARGET_PRECOMPILE_HEADERS(${_T} PRIVATE
        <aspect/global.h> <aspect/plugins.h> <aspect/introspection.h> <aspect/parameters.h>)
      TARGET_PRECOMPILE_HEADERS(${_T} PRIVATE
        <deal.II/base/table_handler.h> <deal.II/base/timer.h>
        <deal.II/base/conditional_ostream.h> <deal.II/distributed/tria.h>
        <deal.II/dofs/dof_handler.h> <deal.II/fe/fe.h> <deal.II/fe/mapping_q.h>
        <deal.II/particles/particle_handler.h>)
    endforeach()
  else()
    message(FATAL_ERROR "ASPECT_PRECOMPILE_HEADERS is currently only supported for CMake 3.16 and newer versions.")
  endif()
endif()

if(ASPECT_UNITY_BUILD)
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16)
    foreach(_T ${TARGET_EXECUTABLES})
      set_property(TARGET ${_T} PROPERTY UNITY_BUILD TRUE)
    endforeach()
    message(STATUS "Combining source files into unity build.")
  elseif(CMAKE_VERSION VERSION_LESS 3.16)
    message(FATAL_ERROR "ASPECT_UNITY_BUILD is currently only supported for CMake 3.16 and newer versions.")
  endif()
else()
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16)
    foreach(_T ${TARGET_EXECUTABLES})
      set_property(TARGET ${_T} PROPERTY UNITY_BUILD FALSE)
    endforeach()
  endif()
endif()


#
## installation
#
# binaries:
install(TARGETS ${TARGET_EXECUTABLES}
  RUNTIME DESTINATION bin
  COMPONENT runtime)


# make sure we have the rpath to our dependencies set:
foreach(_T ${TARGET_EXECUTABLES})
  set_property(TARGET ${_T} PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)
endforeach()

# headers:
install(DIRECTORY include/ ${CMAKE_BINARY_DIR}/include/
  DESTINATION include
  COMPONENT includes
  FILES_MATCHING PATTERN "*.h")

# examples:
set(ASPECT_INSTALL_EXAMPLES OFF CACHE BOOL "If ON all cookbooks and benchmarks will be built and installed.")

if(ASPECT_INSTALL_EXAMPLES)
  add_subdirectory(benchmarks)
  add_subdirectory(cookbooks)

  install(DIRECTORY cookbooks/
    DESTINATION cookbooks
    COMPONENT examples
    FILES_MATCHING PATTERN "*")
  install(DIRECTORY benchmarks/
    DESTINATION benchmarks
    COMPONENT examples
    FILES_MATCHING PATTERN "*")
endif()

# data files:
install(DIRECTORY ${CMAKE_SOURCE_DIR}/data/
  DESTINATION data
  COMPONENT data)

# cmake stuff:
install(FILES ${CMAKE_BINARY_DIR}/forinstall/AspectConfig.cmake ${CMAKE_BINARY_DIR}/AspectConfigVersion.cmake
        DESTINATION "lib/cmake/Aspect/")
install(FILES ${CMAKE_BINARY_DIR}/include/aspect/revision.h DESTINATION "include/aspect/")

message(STATUS "Writing configuration details into detailed.log...")


message(STATUS "")
MESSAGE(STATUS "===== Configuring ASPECT documentation =============")

###########################################################
# Having configured most of the compilation phase, now also
# deal with other parts of the system.
###########################################################

# Start with the documentation

add_subdirectory(doc)

include(cmake/write_config)

# print "info" if run for the first time:
if(NOT USAGE_PRINTED)
  include(${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/print_usage.cmake)
  set(USAGE_PRINTED TRUE CACHE INTERNAL "")
else()
  message(STATUS "Run  ${_make_command} info  to print a detailed help message")
endif()
