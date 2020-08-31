# Copyright (C) 2013 - 2020 by the authors of the ASPECT code.
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

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)

MESSAGE(STATUS "====================================================")
MESSAGE(STATUS "============ Configuring ASPECT ====================")
MESSAGE(STATUS "====================================================")

# Set the name of the project and target:
SET(TARGET "aspect")

FILE(GLOB_RECURSE TARGET_SRC
     "source/*.cc" "unit_tests/*.cc" "include/*.h" "contrib/catch/catch.hpp")

# Special treatment of some files for unity builds. We move the following
# files to the end of the list of all .cc files (TARGET_SRC). This ordering is
# required to not have specializations (for example helper_functions.cc)
# appear after explicit class instantiation (for example core.cc) of the same
# class that happen to be included in the same unity cxx file.
SET(UNITY_LAST_FILES
  "source/simulator/core.cc;source/volume_of_fluid/handler.cc")

FOREACH(_source_file ${UNITY_LAST_FILES})
  SET(_full_name "${CMAKE_SOURCE_DIR}/${_source_file}")
  LIST(FIND TARGET_SRC ${_full_name} _index)
  IF(_index EQUAL -1)
    MESSAGE(FATAL_ERROR "could not find ${_full_name}.")
  ENDIF()

  LIST(REMOVE_ITEM TARGET_SRC ${_full_name})
  LIST(APPEND TARGET_SRC ${_full_name})
ENDFOREACH()

# Set up include directories. Put the ASPECT header files
# to the front of the list, whereas the 'catch' headers can
# be preempted by the system
INCLUDE_DIRECTORIES(BEFORE ${CMAKE_BINARY_DIR}/include include)
INCLUDE_DIRECTORIES(AFTER  contrib/catch)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

LIST(APPEND CMAKE_MODULE_PATH
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/cmake/modules
  )

FIND_PACKAGE(deal.II 9.2.0 QUIET
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

IF(NOT DEAL_II_WITH_CXX11)
  MESSAGE(FATAL_ERROR "\n*** ASPECT requires C++11 but your version of deal.II is not "
    "configured with it. This likely means that your compiler is too old. Check "
    "DEAL_II_WITH_CXX11 in deal.II.")
ENDIF()


SET(ASPECT_USE_PETSC OFF CACHE BOOL "Use PETSc instead of Trilinos if set to 'on'.")
IF(ASPECT_USE_PETSC)
    MESSAGE(FATAL_ERROR "Using PETSc by setting ASPECT_USE_PETSC to ON is currently not supported, see https://github.com/geodynamics/aspect/issues/3338 for details.")
ENDIF()

MESSAGE(STATUS "Using ASPECT_USE_PETSC = '${ASPECT_USE_PETSC}'")

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
    "\nASPECT requires a deal.II installation built with support for Trilinos/PETSc and p4est but one or both of these appears to be missing!\n"
    )
ENDIF()

DEAL_II_INITIALIZE_CACHED_VARIABLES()
PROJECT(${TARGET} CXX)

SET(ASPECT_WITH_WORLD_BUILDER ON CACHE BOOL "Whether to enable compiling aspect with the Geodynamic World Builder.")
MESSAGE(STATUS "ASPECT_WITH_WORLD_BUILDER = '${ASPECT_WITH_WORLD_BUILDER}'")
if(ASPECT_WITH_WORLD_BUILDER)
  # Check whether we can find the WorldBuilder.
  SET(WORLD_BUILDER_SOURCE_DIR "" CACHE PATH "Provide an external World Builder directory to be compiled with ASPECT. \
    If the path is not not provided or the World Builder is not found in the provided location, the version in the contrib folder is is used.")
  IF (NOT EXISTS ${WORLD_BUILDER_SOURCE_DIR}/VERSION)
    MESSAGE(STATUS "World Builder not found. Using internal version.")
    SET(WORLD_BUILDER_SOURCE_DIR "${CMAKE_SOURCE_DIR}/contrib/world_builder/" CACHE PATH "" FORCE)
  ENDIF()

  # Always include the header files.
  INCLUDE_DIRECTORIES("${WORLD_BUILDER_SOURCE_DIR}/include/")

  INCLUDE("${WORLD_BUILDER_SOURCE_DIR}/cmake/version.cmake")
  MESSAGE(STATUS "Using World Builder version ${WORLD_BUILDER_VERSION} found at ${WORLD_BUILDER_SOURCE_DIR}.")

  # add source and include dirs:
  FILE(GLOB_RECURSE wb_files "${WORLD_BUILDER_SOURCE_DIR}/source/*.cc")
  LIST(APPEND TARGET_SRC ${wb_files})

  # generate config.cc and include it:
  CONFIGURE_FILE("${WORLD_BUILDER_SOURCE_DIR}/source/config.cc.in" "${CMAKE_BINARY_DIR}/world_builder_config.cc" @ONLY)
  LIST(INSERT TARGET_SRC 0 "${CMAKE_BINARY_DIR}/world_builder_config.cc")

  # Move some file to the end for unity builds to make sure other file come
  # "before". Note: The current design keeps all ASPECT files (including
  # ASPECT_UNITY_LAST files) before all GWB files. Mixing them will causes
  # many issues with non-unique namespace names like Utilities.
  SET(UNITY_WB_LAST_FILES
  "${WORLD_BUILDER_SOURCE_DIR}/source/parameters.cc")

  FOREACH(_source_file ${UNITY_WB_LAST_FILES})
    LIST(FIND TARGET_SRC ${_source_file} _index)
    IF(_index EQUAL -1)
      MESSAGE(FATAL_ERROR "could not find ${_source_file}.")
    ENDIF()

    LIST(REMOVE_ITEM TARGET_SRC ${_source_file})
    LIST(APPEND TARGET_SRC ${_source_file})
  ENDFOREACH()

  FOREACH(_source_file ${wb_files})
    # exclude the world builder files from including precompiled headers, they
    # do not include ASPECT's dependencies at all.
    SET_PROPERTY(SOURCE ${_source_file} PROPERTY COTIRE_EXCLUDED TRUE )
    SET_PROPERTY(SOURCE ${_source_file} PROPERTY SKIP_PRECOMPILE_HEADERS TRUE )
    # Temporarily disable world builder unity builds:
    SET_PROPERTY(SOURCE ${_source_file} PROPERTY SKIP_UNITY_BUILD_INCLUSION TRUE )
  ENDFOREACH()

  SET_PROPERTY(SOURCE "${CMAKE_BINARY_DIR}/world_builder_config.cc" PROPERTY COTIRE_EXCLUDED TRUE )
  SET_PROPERTY(SOURCE "${CMAKE_BINARY_DIR}/world_builder_config.cc" PROPERTY SKIP_PRECOMPILE_HEADERS TRUE )
ENDIF()



# load in version info and export it
FILE(STRINGS "${CMAKE_SOURCE_DIR}/VERSION" ASPECT_PACKAGE_VERSION LIMIT_COUNT 1)

INCLUDE(${CMAKE_SOURCE_DIR}/cmake/macro_aspect_query_git_information.cmake)
ASPECT_QUERY_GIT_INFORMATION("ASPECT")
CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/include/aspect/revision.h.in ${CMAKE_BINARY_DIR}/include/aspect/revision.h @ONLY)

INCLUDE(CMakePackageConfigHelpers)
WRITE_BASIC_PACKAGE_VERSION_FILE(
  "${CMAKE_BINARY_DIR}/AspectConfigVersion.cmake"
  VERSION ${ASPECT_PACKAGE_VERSION}
  COMPATIBILITY AnyNewerVersion
)

# Configure a cmake fragment that plugins can use to
# set up compiler flags, include paths, etc to compile an
# ASPECT plugin.
# Config for the build dir:
SET(CONFIG_INCLUDE_DIRS "${CMAKE_BINARY_DIR}/include" "${CMAKE_SOURCE_DIR}/include")
SET(CONFIG_DIR "${CMAKE_BINARY_DIR}")
CONFIGURE_FILE(
  ${CMAKE_SOURCE_DIR}/cmake/AspectConfig.cmake.in
  ${CMAKE_BINARY_DIR}/AspectConfig.cmake
  @ONLY
)
# Config for the install dir:
SET(CONFIG_INCLUDE_DIRS "${CMAKE_INSTALL_PREFIX}/include")
SET(CONFIG_DIR "${CMAKE_INSTALL_PREFIX}/bin")
CONFIGURE_FILE(
  ${CMAKE_SOURCE_DIR}/cmake/AspectConfig.cmake.in
  ${CMAKE_BINARY_DIR}/forinstall/AspectConfig.cmake
  @ONLY
)

# Next, set up the testsuite

SET(ASPECT_RUN_ALL_TESTS OFF CACHE BOOL "Set up complete test suite to run.")
SET(ASPECT_NEED_TEST_CONFIGURE ON CACHE BOOL "If true, reconfigure test project.")
SET(ASPECT_COMPARE_TEST_RESULTS ON CACHE BOOL "Compare test results with high accuracy.")

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

# enable all tests and force a re-run of the cmake configuration in the test/ folder:
ADD_CUSTOM_TARGET(setup_tests
  COMMAND ${CMAKE_COMMAND} -D ASPECT_RUN_ALL_TESTS=ON -D ASPECT_NEED_TEST_CONFIGURE=ON . >/dev/null
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Enabling all tests ...")

# Disable the ability to run all tests in an in-source build:
IF (("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}") AND ASPECT_RUN_ALL_TESTS)
  SET(ASPECT_RUN_ALL_TESTS OFF CACHE BOOL "" FORCE)
  MESSAGE(FATAL_ERROR "\nEnabling all tests is not supported in in-source builds. Please create a separate build directory!\n")
ENDIF()


SET(ASPECT_TEST_GENERATOR "Unix Makefiles" CACHE STRING
  "Generator to use for the test cmake project. Using ninja instead of make is not recommended.")


# Allow us to make a global "test" target:
IF(POLICY CMP0037)
  # allow to override "test" target
  CMAKE_POLICY(SET CMP0037 OLD)
ENDIF()

# This is an empty target but we will make it depend on tests/ and unit_tests/ next:
ADD_CUSTOM_TARGET(test)

# Generate CTestTestfile.cmake in the main build folder that lists all subfolders
# that contain tests. This way you can call "ctest" in the build directory.
FILE(WRITE ${CMAKE_BINARY_DIR}/CTestTestfile.cmake "# auto-generated ctest file\n")

IF(EXISTS ${CMAKE_SOURCE_DIR}/unit_tests/CMakeLists.txt)
  # If we have the unit_tests directory, add it so ctest picks it up and configure
  # the test project in the subfolder:
  FILE(APPEND ${CMAKE_BINARY_DIR}/CTestTestfile.cmake "SUBDIRS(unit_tests)\n")
  FILE(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/unit_tests)
  EXECUTE_PROCESS(
    COMMAND ${CMAKE_COMMAND}
        -D CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -D ASPECT_BINARY=${CMAKE_BINARY_DIR}/aspect
        ${CMAKE_CURRENT_SOURCE_DIR}/unit_tests
    OUTPUT_FILE setup_unit_tests.log
    RESULT_VARIABLE test_cmake_result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/unit_tests
   )
  IF(NOT test_cmake_result EQUAL 0)
    MESSAGE(FATAL_ERROR "ERROR: unittest/ project could not be configured.")
  ENDIF()

  # Finally hook up that "make test" will execute the tests:
  ADD_CUSTOM_TARGET(run_unit_tests
    COMMAND ${CMAKE_BINARY_DIR}/aspect --test
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/unit_tests
    DEPENDS aspect
    COMMENT "Running unit_tests ...")

  ADD_DEPENDENCIES(test run_unit_tests)
ENDIF()

IF(EXISTS ${CMAKE_SOURCE_DIR}/tests/CMakeLists.txt
  AND ${ASPECT_NEED_TEST_CONFIGURE})

  # Hook up the tests:
  FILE(APPEND ${CMAKE_BINARY_DIR}/CTestTestfile.cmake "SUBDIRS(tests)\n")

  FILE(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/tests)

  SET(aspect_binary ${CMAKE_BINARY_DIR}/aspect)

  # set a flag so we don't need to rerun this configuration step every time:
  SET(ASPECT_NEED_TEST_CONFIGURE OFF CACHE BOOL "" FORCE)

  MESSAGE(STATUS "Setting up test project, see tests/setup_tests.log for details.")
  EXECUTE_PROCESS(
    COMMAND ${CMAKE_COMMAND} -G ${ASPECT_TEST_GENERATOR}
        -D ASPECT_RUN_ALL_TESTS=${ASPECT_RUN_ALL_TESTS}
        -D ASPECT_COMPARE_TEST_RESULTS=${ASPECT_COMPARE_TEST_RESULTS}
        -D Aspect_DIR=${CMAKE_BINARY_DIR}
        -D CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -D ASPECT_BINARY=${aspect_binary}
        ${CMAKE_CURRENT_SOURCE_DIR}/tests
    OUTPUT_FILE setup_tests.log
    RESULT_VARIABLE test_cmake_result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/tests
   )
  IF(NOT test_cmake_result EQUAL 0)
    MESSAGE(FATAL_ERROR "ERROR: tests/ project could not be configured.")
  ENDIF()

  # Finally hook up that "make test" will execute the tests:
  ADD_CUSTOM_TARGET(run_tests
    COMMAND ${CMAKE_COMMAND} --build . --target test
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/tests
    DEPENDS aspect
    COMMENT "Running tests ...")

  ADD_DEPENDENCIES(test run_tests)
ENDIF()

# Provide the "generate_reference_output" target:
ADD_CUSTOM_TARGET(generate_reference_output
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/cmake/generate_reference_output.sh
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})


# Provide "indent" target for indenting all headers and source files
ADD_CUSTOM_TARGET(indent
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  COMMAND ./contrib/utilities/indent
  COMMENT "Indenting all ASPECT header and source files..."
  )

IF(CMAKE_GENERATOR MATCHES "Ninja")
  SET(_make_command "$ ninja")
ELSE()
  SET(_make_command " $ make")
ENDIF()

# Provide "release" and "debug" targets to switch compile mode
IF(${DEAL_II_BUILD_TYPE} MATCHES "DebugRelease")
ADD_CUSTOM_TARGET(release
  COMMAND ${CMAKE_COMMAND} -D CMAKE_BUILD_TYPE=Release .
  COMMAND ${CMAKE_COMMAND} -E echo "***"
  COMMAND ${CMAKE_COMMAND} -E echo "*** Switched to Release mode. Now recompile with: ${_make_command}"
  COMMAND ${CMAKE_COMMAND} -E echo "***"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  VERBATIM
  COMMENT "switching to RELEASE mode..."
  )
ADD_CUSTOM_TARGET(debug
  COMMAND ${CMAKE_COMMAND} -D CMAKE_BUILD_TYPE=Debug .
  COMMAND ${CMAKE_COMMAND} -E echo "***"
  COMMAND ${CMAKE_COMMAND} -E echo "*** Switched to Debug mode. Now recompile with: ${_make_command}"
  COMMAND ${CMAKE_COMMAND} -E echo "***"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  VERBATIM
  COMMENT "switching to DEBUG mode..."
  )
ENDIF()

# Provide a "distclean" target (like it is done in deal.II):
ADD_CUSTOM_TARGET(distclean
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target clean
  COMMAND ${CMAKE_COMMAND} -E remove_directory CMakeFiles
  COMMAND ${CMAKE_COMMAND} -E remove
    CMakeCache.txt cmake_install.cmake Makefile
    build.ninja rules.ninja .ninja_deps .ninja_log
  COMMENT "distclean invoked"
  )



FILE(WRITE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/print_usage.cmake
"MESSAGE(
\"###
#
#  Project ${TARGET} set up with  ${DEAL_II_PACKAGE_NAME}-${DEAL_II_PACKAGE_VERSION}  found at
#      ${DEAL_II_PATH}
#
#  CMAKE_BUILD_TYPE:          ${CMAKE_BUILD_TYPE}
#
#  You can now run
#      ${_make_command}                - to compile and link ${TARGET}
#      ${_make_command} debug          - to switch the build type to 'Debug'
#      ${_make_command} release        - to switch the build type to 'Release'
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
ADD_CUSTOM_TARGET(info
  COMMAND ${CMAKE_COMMAND} -P ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/print_usage.cmake
  )


# Depending on whether we link statically or allow for shared libs,
# we can or can not load plugins via external shared libs. Pass this
# down during compilation so we can disable it in the code
SET(ASPECT_USE_SHARED_LIBS ON CACHE BOOL "If ON, we support loading shared plugin files.")
IF (DEAL_II_STATIC_EXECUTABLE STREQUAL "ON")
  MESSAGE(STATUS "Creating a statically linked executable")
  SET(ASPECT_USE_SHARED_LIBS OFF CACHE BOOL "" FORCE)
ENDIF()

INCLUDE (CheckCXXSourceCompiles)

SET(_backup_libs ${CMAKE_REQUIRED_LIBRARIES})
LIST(APPEND CMAKE_REQUIRED_LIBRARIES ${CMAKE_DL_LIBS})
CHECK_CXX_SOURCE_COMPILES("
#include <cstddef>
#include <dlfcn.h>

int main()
{
  void *handle = dlopen (\"somelib.so\", RTLD_LAZY);
  return handle == NULL || dlerror();
}
" HAVE_DLOPEN)
SET(CMAKE_REQUIRED_LIBRARIES ${_backup_libs})

IF (NOT HAVE_DLOPEN)
  MESSAGE(STATUS "dlopen() test failed, disabling dynamic plugin loading")
  SET(ASPECT_USE_SHARED_LIBS OFF CACHE BOOL "" FORCE)
ENDIF()

IF (ASPECT_USE_SHARED_LIBS)
  MESSAGE(STATUS "Enabling dynamic loading of plugins from the input file")
ELSE()
  MESSAGE(STATUS "Disabling dynamic loading of plugins from the input file")
ENDIF()

# See whether we can verify that every plugin we load is compiled against
# the same deal.II library
SET(ASPECT_HAVE_LINK_H ON CACHE BOOL "If ON, link.h exists and is usable.")
INCLUDE (CheckIncludeFileCXX)
CHECK_INCLUDE_FILE_CXX ("link.h" _HAVE_LINK_H)
IF (NOT _HAVE_LINK_H)
  SET(ASPECT_HAVE_LINK_H OFF CACHE BOOL "" FORCE)
ENDIF()
IF (ASPECT_HAVE_LINK_H)
  MESSAGE(STATUS "Enabling checking of compatible deal.II library when loading plugins")
ENDIF()

SET(ASPECT_ADDITIONAL_CXX_FLAGS "" CACHE STRING "Additional CMAKE_CXX_FLAGS applied after the deal.II options.")

IF(NOT ASPECT_ADDITIONAL_CXX_FLAGS STREQUAL "")
  MESSAGE(STATUS "Appending ASPECT_ADDITIONAL_CXX_FLAGS: '${ASPECT_ADDITIONAL_CXX_FLAGS}':")
  STRING(APPEND DEAL_II_CXX_FLAGS_DEBUG " ${ASPECT_ADDITIONAL_CXX_FLAGS}")
  STRING(APPEND DEAL_II_CXX_FLAGS_RELEASE " ${ASPECT_ADDITIONAL_CXX_FLAGS}")
  MESSAGE(STATUS "  DEAL_II_CXX_FLAGS_DEBUG: ${DEAL_II_CXX_FLAGS_DEBUG}")
  MESSAGE(STATUS "  DEAL_II_CXX_FLAGS_RELEASE: ${DEAL_II_CXX_FLAGS_RELEASE}")
ENDIF()

ADD_EXECUTABLE(${TARGET} ${TARGET_SRC})
DEAL_II_SETUP_TARGET(${TARGET})

FIND_PACKAGE(PerpleX QUIET HINTS ./contrib/perplex/install/ ../ ../../ $ENV{PERPLEX_DIR})
IF(${PerpleX_FOUND})
  MESSAGE(STATUS "PerpleX found at ${PerpleX_INCLUDE_DIR}")
  INCLUDE_DIRECTORIES(${PerpleX_INCLUDE_DIR})
  TARGET_LINK_LIBRARIES(${TARGET} ${PerpleX_LIBRARIES})
  SET(ASPECT_WITH_PERPLEX ON)
ELSE()
  SET(ASPECT_WITH_PERPLEX OFF)
ENDIF()


# Find the libdap package so that ASPECT can connect to the OPeNDAP servers
# Author: Kodi Neumiller
SET(ASPECT_WITH_LIBDAP OFF CACHE BOOL "Check if the user wants to compile ASPECT with the libdap libraries.")
MESSAGE(STATUS "Using ASPECT_WITH_LIBDAP = '${ASPECT_WITH_LIBDAP}'")
IF(ASPECT_WITH_LIBDAP)
  FIND_PACKAGE(LIBDAP)
  IF(${LIBDAP_FOUND})
    INCLUDE_DIRECTORIES(${LIBDAP_INCLUDE_DIRS})
    TARGET_LINK_LIBRARIES(${TARGET} ${LIBDAP_LIBRARIES})
    MESSAGE(STATUS "LIBDAP found at ${LIBDAP_LIBRARY}\n")
  ELSE()
    MESSAGE(FATAL_ERROR "LIBDAP not found. Disable ASPECT_WITH_LIBDAP or specify a hint to your installation directory with LIBDAP_DIR.")
  ENDIF()
ENDIF()

#
# NETCDF (c including parallel)
#
SET(ASPECT_WITH_NETCDF OFF CACHE BOOL "Check if the user wants to compile ASPECT with the NETCDF libraries.")

MESSAGE(STATUS "Using ASPECT_WITH_NETCDF = '${ASPECT_WITH_NETCDF}'")

IF(ASPECT_WITH_NETCDF)
  FIND_PACKAGE(NETCDF)
  IF(${NETCDF_FOUND})
    MESSAGE(STATUS "  NETCDF_INCLUDE_DIR: ${NETCDF_INCLUDE_DIR}")
    MESSAGE(STATUS "  NETCDF_LIBRARY: ${NETCDF_LIBRARY}")
    MESSAGE(STATUS "  NETCDF_VERSION: ${NETCDF_VERSION}")

    INCLUDE_DIRECTORIES(${NETCDF_INCLUDE_DIRS})
    TARGET_LINK_LIBRARIES(${TARGET} ${NETCDF_LIBRARIES})
  ELSE()
    MESSAGE(FATAL_ERROR "NETCDF not found. Disable ASPECT_WITH_NETCDF or specify a hint to your installation directory with NETCDF_DIR.")
  ENDIF()
ENDIF()


IF (ASPECT_USE_SHARED_LIBS)
  # some systems need to explicitly link to some libraries to use dlopen
  TARGET_LINK_LIBRARIES(aspect ${CMAKE_DL_LIBS})
ENDIF()



# Check if we can raise floating point exceptions.
#
# Note that some library we link with in ASPECT on some platforms will trigger
# floating point exceptions when converting -numeric_limits<double>::max to a
# string. The only thing we can do is a configure check and disable the
# exceptions. This is done here:
SET(ASPECT_USE_FP_EXCEPTIONS ON CACHE BOOL "If ON, floating point exception are raised in debug mode.")

# Clang 6.0 throws random floating point exceptions, which we could not
# track down. Disable the exceptions for now.
IF("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" AND NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 6.0.0)
  SET(ASPECT_USE_FP_EXCEPTIONS OFF CACHE BOOL "" FORCE)
ENDIF()

IF (ASPECT_USE_FP_EXCEPTIONS)
  INCLUDE(${CMAKE_SOURCE_DIR}/cmake/fpe_check.cmake)

  IF (HAVE_FP_EXCEPTIONS)
    MESSAGE(STATUS "Runtime floating point checks enabled.")
  ELSE()
    SET(ASPECT_USE_FP_EXCEPTIONS OFF CACHE BOOL "" FORCE)
    MESSAGE(STATUS "No support for feenableexcept(), disabling runtime floating point exception checks.")
  ENDIF()
ENDIF()


#
# generate config.h
#

# Define macro with the source directory that will be exported in
# config.h. This can be used to hard-code the location of data files, such as
# in $ASPECT_SOURCE_DIR/data/velocity-boundary-conditions/gplates/*
SET(ASPECT_SOURCE_DIR ${CMAKE_SOURCE_DIR})

# And finally generate the file
CONFIGURE_FILE(
  ${CMAKE_SOURCE_DIR}/include/aspect/config.h.in
  ${CMAKE_BINARY_DIR}/include/aspect/config.h
  )


# Check if we want to precompile header files. This speeds up compile time,
# but can fail on some machines with old cmake, therefore it is disabled by
# default. Also cotire currently has issues with cmake generator expressions
# used in deal.II for cmake version 3.9 and above. Disable cotire for cmake >=
# 3.9. Starting with cmake 3.16 there is native support inside cmake and we can
# use that.
IF(CMAKE_VERSION VERSION_LESS 3.9)
  SET(ASPECT_PRECOMPILE_HEADERS OFF CACHE BOOL "Precompile external header files using cotire to speedup compile time. Currently only supported for cmake 3.8 and older versions, or 3.16 and newer versions.")
  SET(ASPECT_UNITY_BUILD OFF CACHE BOOL "Combine source files into less compile targets to speedup compile time. Currently only supported for cmake 3.16 and newer versions." FORCE)
ELSEIF(NOT CMAKE_VERSION VERSION_LESS 3.16)
  SET(ASPECT_PRECOMPILE_HEADERS ON CACHE BOOL "Precompile external header files to speedup compile time. Currently only supported for cmake 3.8 and older versions or 3.16 and newer versions.")
  SET(ASPECT_UNITY_BUILD ON CACHE BOOL "Combine source files into less compile targets to speedup compile time. Currently only supported for cmake 3.16 and newer versions.")
ELSE()
  SET(ASPECT_PRECOMPILE_HEADERS OFF CACHE BOOL "Precompile external header files to speedup compile time. Currently only supported for cmake 3.8 and older versions or 3.16 and newer versions." FORCE)
  SET(ASPECT_UNITY_BUILD OFF CACHE BOOL "Combine source files into less compile targets to speedup compile time. Currently only supported for cmake 3.16 and newer versions." FORCE)
ENDIF()

IF (ASPECT_PRECOMPILE_HEADERS)
  IF (CMAKE_VERSION VERSION_LESS 3.9)
    MESSAGE(STATUS "Enabling cotire to precompile external header files.")
    INCLUDE(cmake/cotire)

    # Set some properties that are necessary for the header precompilations
    # and unity build to pass. We need to exclude some headers that
    # frequently cause problems, and undefine some macros between source
    # files that would otherwise interfere with each other.
    SET_TARGET_PROPERTIES (aspect PROPERTIES
        COTIRE_PREFIX_HEADER_IGNORE_PATH
        "/usr/include/;${CMAKE_SOURCE_DIR};${CMAKE_BINARY_DIR}")

    # This line activates the cotire module for the aspect target,
    # and therefore the whole precompilation
    cotire(aspect)

  ELSEIF (NOT CMAKE_VERSION VERSION_LESS 3.16)
    MESSAGE(STATUS "Precompiling common header files.")
    # Use the native cmake support to precompile some common headers
    # from ASPECT and deal.II that are frequently included, but rarely changed.
    TARGET_PRECOMPILE_HEADERS(aspect PRIVATE
      <aspect/global.h> <aspect/plugins.h> <aspect/introspection.h> <aspect/parameters.h>)
    TARGET_PRECOMPILE_HEADERS(aspect PRIVATE
      <deal.II/base/table_handler.h> <deal.II/base/timer.h>
      <deal.II/base/conditional_ostream.h> <deal.II/distributed/tria.h>
      <deal.II/dofs/dof_handler.h> <deal.II/fe/fe.h> <deal.II/fe/mapping_q.h>
      <deal.II/particles/particle_handler.h>)

  ELSE()
    MESSAGE(FATAL_ERROR "ASPECT_PRECOMPILE_HEADERS is currently only supported for CMake 3.8 and older versions or 3.16 and newer versions.")
  ENDIF()

ELSE()
  MESSAGE(STATUS "Disabling precompiling headers.")
ENDIF()


IF (ASPECT_UNITY_BUILD)
  IF (NOT CMAKE_VERSION VERSION_LESS 3.16)
    SET_PROPERTY(TARGET aspect PROPERTY UNITY_BUILD TRUE)
    MESSAGE(STATUS "Combining source files into unity build.")
  ELSEIF(CMAKE_VERSION VERSION_LESS 3.16)
    MESSAGE(FATAL_ERROR "ASPECT_UNITY_BUILD is currently only supported for CMake 3.16 and newer versions.")
  ENDIF()
ELSE()
  IF (NOT CMAKE_VERSION VERSION_LESS 3.16)
    SET_PROPERTY(TARGET aspect PROPERTY UNITY_BUILD FALSE)
  ENDIF()

  MESSAGE(STATUS "Disabling unity build.")
ENDIF()


# Find the deal.II parameter GUI and install helper script
FIND_PROGRAM(PARAMETER_GUI_EXECUTABLE
        parameter_gui
        HINTS $ENV{PARAMETER_GUI_DIR} $ENV{PARAMETER_GUI_DIR}/bin ${PARAMETER_GUI_DIR} ${PARAMETER_GUI_DIR}/bin
        PATH bin)
MARK_AS_ADVANCED(CLEAR PARAMETER_GUI_EXECUTABLE)

# did the user specify something that doesn't exist?
IF (PARAMETER_GUI_EXECUTABLE
        AND
    (NOT EXISTS ${PARAMETER_GUI_EXECUTABLE} OR IS_DIRECTORY ${PARAMETER_GUI_EXECUTABLE}))
  MESSAGE(STATUS "Warning: PARAMETER_GUI_EXECUTABLE '${PARAMETER_GUI_EXECUTABLE}' does not exist")
  SET(PARAMETER_GUI_EXECUTABLE "PARAMETER_GUI_EXECUTABLE-NOTFOUND" CACHE FILEPATH "" FORCE)
ENDIF()

IF (NOT PARAMETER_GUI_EXECUTABLE)
  MESSAGE(STATUS "Parameter GUI not found: install and provide a hint using -D PARAMETER_GUI_DIR or set -D PARAMETER_GUI_EXECUTABLE directly.")
ELSE()
  MESSAGE(STATUS "Found parameter GUI at: ${PARAMETER_GUI_EXECUTABLE}")
  CONFIGURE_FILE(
  ${CMAKE_SOURCE_DIR}/doc/aspect-gui
  ${CMAKE_BINARY_DIR}/aspect-gui
  @ONLY
  )
INSTALL(FILES ${CMAKE_BINARY_DIR}/aspect-gui
        DESTINATION bin
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
                GROUP_READ GROUP_EXECUTE
                WORLD_READ WORLD_EXECUTE)
ENDIF()

#
## installation
#
# binary:
INSTALL(TARGETS aspect
  RUNTIME DESTINATION bin
  COMPONENT runtime)

# make sure we have the rpath to our dependencies set:
SET_PROPERTY(TARGET aspect PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)

# headers:
INSTALL(DIRECTORY include/
  DESTINATION include
  COMPONENT includes
  FILES_MATCHING PATTERN "*.h")

# cmake stuff:
INSTALL(FILES ${CMAKE_BINARY_DIR}/forinstall/AspectConfig.cmake ${CMAKE_BINARY_DIR}/AspectConfigVersion.cmake
        DESTINATION "lib/cmake/Aspect/")

INSTALL(FILES ${CMAKE_BINARY_DIR}/include/aspect/revision.h DESTINATION "include/aspect/")

MESSAGE(STATUS "Writing config into detailed.log...")
INCLUDE(cmake/write_config)

# print "info" if run for the first time:
IF(NOT USAGE_PRINTED)
  INCLUDE(${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/print_usage.cmake)
  SET(USAGE_PRINTED TRUE CACHE INTERNAL "")
ELSE()
  MESSAGE(STATUS "Run  ${_make_command} info  to print a detailed help message")
ENDIF()
