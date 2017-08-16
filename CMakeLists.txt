# Copyright (C) 2013 - 2016 by the authors of the ASPECT code.
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

LIST(APPEND CMAKE_MODULE_PATH
  ${CMAKE_SOURCE_DIR}
  )

FIND_PACKAGE(deal.II 8.5.0 QUIET
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


# load in version info and export it
FILE(STRINGS "${CMAKE_SOURCE_DIR}/VERSION" ASPECT_PACKAGE_VERSION LIMIT_COUNT 1)

INCLUDE(CMakePackageConfigHelpers)
WRITE_BASIC_PACKAGE_VERSION_FILE(
  "${CMAKE_BINARY_DIR}/AspectConfigVersion.cmake"
  VERSION ${ASPECT_PACKAGE_VERSION}
  COMPATIBILITY AnyNewerVersion
)

# Configure a cmake fragment that plugins can use to
# set up compiler flags, include paths, etc to compile an
# Aspect plugin.
# Config for the build dir:
SET(CONFIG_INCLUDE_DIRS "${CMAKE_SOURCE_DIR}/include")
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

SET(ASPECT_TEST_GENERATOR "Unix Makefiles" CACHE STRING 
  "Generator to use for the test cmake project. Using ninja instead of make is not recommended.")

IF(EXISTS ${CMAKE_SOURCE_DIR}/tests/CMakeLists.txt)
  FILE(WRITE ${CMAKE_BINARY_DIR}/CTestTestfile.cmake "SUBDIRS(tests)")

  FILE(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/tests)

  SET(aspect_binary ${CMAKE_BINARY_DIR}/aspect)

  MESSAGE(STATUS "Setting up test project, see tests/setup_tests.log for details.")
  EXECUTE_PROCESS(
    COMMAND ${CMAKE_COMMAND} -G ${ASPECT_TEST_GENERATOR}
	-D ASPECT_RUN_ALL_TESTS=${ASPECT_RUN_ALL_TESTS}
	-D ASPECT_DIR=${CMAKE_BINARY_DIR}
	-D CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
	-D ASPECT_BINARY=${aspect_binary}
	${CMAKE_CURRENT_SOURCE_DIR}/tests
    OUTPUT_FILE setup_tests.log
    RESULT_VARIABLE test_cmake_result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/tests
   )
  IF(NOT test_cmake_result EQUAL 0)
    MESSAGE(FATAL_ERROR "ERROR: Test project could not be configured.")
  ENDIF()


  IF(POLICY CMP0037)
    # allow to override "test" target
    CMAKE_POLICY(SET CMP0037 OLD)
  ENDIF()

  ADD_CUSTOM_TARGET(test
    COMMAND ${CMAKE_COMMAND} --build . --target test
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/tests
    DEPENDS aspect
    COMMENT "Running tests ...")

ENDIF()


ADD_CUSTOM_TARGET(generate_reference_output
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/cmake/generate_reference_output.sh
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})


# Provide "indent" target for indenting all headers and source files
ADD_CUSTOM_TARGET(indent
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  COMMAND ./doc/indent
  COMMENT "Indenting all ASPECT header and source files..."
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
  FOREACH(_source_file ${TARGET_SRC})
    SET_PROPERTY(SOURCE ${_source_file}
      APPEND PROPERTY COMPILE_DEFINITIONS ASPECT_USE_SHARED_LIBS=1)
  ENDFOREACH()
ELSE()
  MESSAGE(STATUS "Disabling dynamic loading of plugins from the input file")
  FOREACH(_source_file ${TARGET_SRC})
    SET_PROPERTY(SOURCE ${_source_file}
      APPEND PROPERTY COMPILE_DEFINITIONS ASPECT_USE_SHARED_LIBS=0)
  ENDFOREACH()
ENDIF()

# See whether we can verify that every plugin we load is compiled against
# the same deal.II library
SET(ASPECT_HAVE_LINK_H ON CACHE BOOL "If ON, link.h exists and is usable.")
INCLUDE (CheckIncludeFiles)
CHECK_INCLUDE_FILES ("link.h" _HAVE_LINK_H)
IF (NOT _HAVE_LINK_H)
  SET(ASPECT_HAVE_LINK_H OFF CACHE BOOL "" FORCE)
ENDIF()
IF (ASPECT_HAVE_LINK_H)
  MESSAGE(STATUS "Enabling checking of compatible deal.II library when loading plugins")
  FOREACH(_source_file ${TARGET_SRC})
    SET_PROPERTY(SOURCE ${_source_file}
      APPEND PROPERTY COMPILE_DEFINITIONS ASPECT_HAVE_LINK_H=1)
  ENDFOREACH()
ENDIF()

# Check if we want to use cotire to precompile header files and prepare a
# unity build target. This speeds up compile time, but can fail on some
# machines, therefore it is disabled by default.
SET(ASPECT_PRECOMPILE_HEADERS OFF CACHE BOOL "Precompile external header files using cotire to speedup compile time.")

IF (ASPECT_PRECOMPILE_HEADERS)
# For the unity build we need to compile helper_functions.cc before core.cc
# ensure this by reordering the files
  FIND_FILE(HELPER_PATH helper_functions.cc HINTS "${CMAKE_SOURCE_DIR}/source/simulator")
  LIST(REMOVE_ITEM TARGET_SRC ${HELPER_PATH})
  LIST(INSERT TARGET_SRC 0 ${HELPER_PATH})
ENDIF()


DEAL_II_INVOKE_AUTOPILOT()


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

IF (ASPECT_PRECOMPILE_HEADERS)
  MESSAGE(STATUS "Enabling cotire to precompile external header files.")
  INCLUDE(cmake/cotire)

  # Set some properties that are necessary for the header precompilations 
  # and unity build to pass. We need to exclude some headers that
  # frequently cause problems, and undefine some macros between source 
  # files that would otherwise interfere with each other.
  set_target_properties (aspect PROPERTIES
      COTIRE_PREFIX_HEADER_IGNORE_PATH
      "/usr/include/;${CMAKE_SOURCE_DIR};${CMAKE_BINARY_DIR}")
  set_source_files_properties (${TARGET_SRC} PROPERTIES
      COTIRE_UNITY_SOURCE_POST_UNDEFS "INSTANTIATE")

  # This line activates the cotire module for the aspect target, 
  # and therefore the whole precompilation
  cotire(aspect)
ELSE()
  MESSAGE(STATUS "Disabling cotire.")
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
INSTALL(DIRECTORY include/ DESTINATION include
  COMPONENT includes)

# cmake stuff:
INSTALL(FILES ${CMAKE_BINARY_DIR}/forinstall/AspectConfig.cmake ${CMAKE_BINARY_DIR}/AspectConfigVersion.cmake
	DESTINATION "lib/cmake/Aspect/")

MESSAGE(STATUS "Writing config into detailed.log...")
INCLUDE(cmake/write_config)

