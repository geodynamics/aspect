#
# Copyright (C) 2013, 2014 by Matthias Maier
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
#
# $Id$


########################################################################
#                                                                      #
#                             Test setup:                              #
#                                                                      #
########################################################################

#
# This is the ctest script for running and submitting build and regression
# tests.
#
# Invoke it in a _build directory_ (or designated build directory) via:
#
#   ctest -S <...>/run_testsuite.cmake
#
# The following configuration variables can be overwritten with
#
#   ctest -D<variable>=<value> [...]
#
#
#   CTEST_SOURCE_DIRECTORY
#     - The source directory of Aspect
#       Note: This is _not_ the test directory ending in "[...]/tests"
#     - If unspecified, "../aspect" relative to the location of this
#       script is used. If this is not a source directory, an error is
#       thrown.
#
#   CTEST_BINARY_DIRECTORY
#     - The designated build directory (already configured, empty, or non
#       existent - see the information about TRACKs what will happen)
#     - If unspecified the current directory is used. If the current
#       directory is equal to CTEST_SOURCE_DIRECTORY or the "tests"
#       directory, an error is thrown.
#
#   CTEST_CMAKE_GENERATOR
#     - The CMake Generator to use (e.g. "Unix Makefiles", or "Ninja", see
#       $ man cmake)
#     - If unspecified the generator of a configured build directory will
#       be used, otherwise "Unix Makefiles".
#
#   TRACK
#     - The track the test should be submitted to. Defaults to
#       "Experimental". Possible values are:
#
#       "Experimental"     - all tests that are not specifically "build" or
#                            "regression" tests should go into this track
#
#       "Build Tests"      - Build tests that configure and build in a
#                            clean directory and run the build tests
#                            "build_tests/*"
#
#       "Nightly"          - Reserved for nightly regression tests for
#                            build bots on various architectures
#
#       "Regression Tests" - Reserved for the regression tester
#
#   CONFIG_FILE
#     - A configuration file (see ../deal.II/docs/development/Config.sample)
#       that will be used during the configuration stage (invokes
#       # cmake -C ${CONFIG_FILE}). This only has an effect if
#       CTEST_BINARY_DIRECTORY is empty.
#
#   DESCRIPTION
#     - A string that is appended to CTEST_BUILD_NAME
#
#   MAKEOPTS
#     - Additional options that will be passed directly to make (or ninja).
#
#   SUBMIT
#     - default OFF, set to ON to submit to our cdash instance
#
# Furthermore, the following variables controlling the testsuite can be set
# and will be automatically handed down to cmake:
#
#   TEST_DIFF
#   TEST_TIME_LIMIT
#   TEST_PICKUP_REGEX
#   TEST_OVERRIDE_LOCATION
#
# For details, consult the ./README file.
#

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)
MESSAGE("-- This is CTest ${CMAKE_VERSION}")

#
# TRACK: Default to Experimental:
#

IF("${TRACK}" STREQUAL "")
  SET(TRACK "Experimental")
ENDIF()

IF( NOT "${TRACK}" STREQUAL "Experimental"
    AND NOT "${TRACK}" STREQUAL "Build Tests"
    AND NOT "${TRACK}" STREQUAL "Nightly"
    AND NOT "${TRACK}" STREQUAL "Regression Tests" )
  MESSAGE(FATAL_ERROR "
Unknown TRACK \"${TRACK}\" - see the manual for valid values.
"
    )
ENDIF()

MESSAGE("-- TRACK:                  ${TRACK}")

#
# CTEST_SOURCE_DIRECTORY:
#

IF("${CTEST_SOURCE_DIRECTORY}" STREQUAL "")
  #
  # If CTEST_SOURCE_DIRECTORY is not set we just assume that this script
  # resides in the top level source directory
  #
  SET(CTEST_SOURCE_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})

  # The code above set the source directory of that of the run_testsuite.cmake
  # script, but we need the directory of aspect, which is simply one level
  # higher
  IF ("${CTEST_SOURCE_DIRECTORY}" MATCHES "/tests")
    SET(CTEST_SOURCE_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/../)
  ENDIF()

  IF(NOT EXISTS ${CTEST_SOURCE_DIRECTORY}/CMakeLists.txt)
    MESSAGE(FATAL_ERROR "
Could not find a suitable source directory. Please, set
CTEST_SOURCE_DIRECTORY manually to the appropriate source directory.
"
      )
  ENDIF()
ENDIF()

MESSAGE("-- CTEST_SOURCE_DIRECTORY: ${CTEST_SOURCE_DIRECTORY}")

#
# Read in custom config files:
#

CTEST_READ_CUSTOM_FILES(${CTEST_SOURCE_DIRECTORY})

#
# CTEST_BINARY_DIRECTORY:
#

IF("${CTEST_BINARY_DIRECTORY}" STREQUAL "")
  #
  # If CTEST_BINARY_DIRECTORY is not set we just use the current directory
  # except if it is equal to CTEST_SOURCE_DIRECTORY in which case we fail.
  #
  SET(CTEST_BINARY_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
ENDIF()

#
# Read in custom config files:
#

CTEST_READ_CUSTOM_FILES(${CTEST_BINARY_DIRECTORY})

# Make sure that for a build test the directory is empty:
FILE(GLOB _test ${CTEST_BINARY_DIRECTORY}/*)
IF( "${TRACK}" STREQUAL "Build Tests"
    AND NOT "${_test}" STREQUAL "" )
      MESSAGE(FATAL_ERROR "
TRACK was set to \"Build Tests\" which require an empty build directory.
But files were found in \"${CTEST_BINARY_DIRECTORY}\"
"
        )
ENDIF()

MESSAGE("-- CTEST_BINARY_DIRECTORY: ${CTEST_BINARY_DIRECTORY}")

#
# CTEST_CMAKE_GENERATOR:
#

# Query Generator from build directory (if possible):
IF(EXISTS ${CTEST_BINARY_DIRECTORY}/CMakeCache.txt)
  FILE(STRINGS ${CTEST_BINARY_DIRECTORY}/CMakeCache.txt _generator
    REGEX "^CMAKE_GENERATOR:"
    )
  STRING(REGEX REPLACE "^.*=" "" _generator ${_generator})
ENDIF()

IF("${CTEST_CMAKE_GENERATOR}" STREQUAL "")
  IF(NOT "${_generator}" STREQUAL "")
    SET(CTEST_CMAKE_GENERATOR ${_generator})
  ELSE()
    # default to "Unix Makefiles"
    SET(CTEST_CMAKE_GENERATOR "Unix Makefiles")
  ENDIF()
ELSE()
  # ensure that CTEST_CMAKE_GENERATOR (that was apparantly set) is
  # compatible with the build directory:
  IF( NOT "${CTEST_CMAKE_GENERATOR}" STREQUAL "${_generator}"
      AND NOT "${_generator}" STREQUAL "" )
    MESSAGE(FATAL_ERROR "
The build directory is already set up with Generator \"${_generator}\", but
CTEST_CMAKE_GENERATOR was set to a different Generator \"${CTEST_CMAKE_GENERATOR}\".
"
     )
  ENDIF()
ENDIF()

MESSAGE("-- CTEST_CMAKE_GENERATOR:  ${CTEST_CMAKE_GENERATOR}")

#
# CTEST_SITE:
#

FIND_PROGRAM(HOSTNAME_COMMAND NAMES hostname)
IF(NOT "${HOSTNAME_COMMAND}" MATCHES "-NOTFOUND")
  EXEC_PROGRAM(${HOSTNAME_COMMAND} OUTPUT_VARIABLE _hostname)
  SET(CTEST_SITE "${_hostname}")
ELSE()
  # Well, no hostname available. What about:
  SET(CTEST_SITE "BobMorane")
ENDIF()

MESSAGE("-- CTEST_SITE:             ${CTEST_SITE}")

IF( "${TRACK}" STREQUAL "Regression Tests"
    AND NOT CTEST_SITE MATCHES "tester" )
  MESSAGE(FATAL_ERROR "
I'm sorry ${CTEST_SITE}, I'm afraid I can't do that.
The TRACK \"Regression Tests\" is not for you.
"
    )
ENDIF()

#
# Assemble configuration options, we need it now:
#

IF("${MAKEOPTS}" STREQUAL "")
  SET(MAKEOPTS $ENV{MAKEOPTS})
ENDIF()

IF(NOT "${CONFIG_FILE}" STREQUAL "")
  SET(_options "-C${CONFIG_FILE}")
ENDIF()

IF("${TRACK}" STREQUAL "Build Tests")
  SET(TEST_PICKUP_REGEX "^build_tests")
ENDIF()

# Pass all relevant variables down to configure:
GET_CMAKE_PROPERTY(_variables VARIABLES)

#
# CTEST_BUILD_NAME:
#

# Append compiler information to CTEST_BUILD_NAME:
IF(NOT EXISTS ${CTEST_BINARY_DIRECTORY}/detailed.log)
  MESSAGE(FATAL_ERROR "could not find detailed.log")
ENDIF()

IF(EXISTS ${CTEST_BINARY_DIRECTORY}/detailed.log)
  FILE(STRINGS ${CTEST_BINARY_DIRECTORY}/detailed.log _compiler_id
    REGEX "CMAKE_CXX_COMPILER:"
    )
  STRING(REGEX REPLACE
    "^.*CMAKE_CXX_COMPILER:     \(.*\) on platform.*$" "\\1"
    _compiler_id ${_compiler_id}
    )
  STRING(REGEX REPLACE "^\(.*\) .*$" "\\1" _compiler_name ${_compiler_id})
  STRING(REGEX REPLACE "^.* " "" _compiler_version ${_compiler_id})
  STRING(REGEX REPLACE " " "-" _compiler_id ${_compiler_id})
  IF( NOT "${_compiler_id}" STREQUAL "" OR
      _compiler_id MATCHES "CMAKE_CXX_COMPILER" )
    SET(CTEST_BUILD_NAME "${_compiler_id}")
  ENDIF()


  FILE(STRINGS ${CTEST_BINARY_DIRECTORY}/detailed.log _use_petsc_line
    REGEX "ASPECT_USE_PETSC:"
    )
  STRING(REGEX REPLACE "^.*#.*ASPECT_USE_PETSC: *(.*)$" "\\1" USE_PETSC ${_use_petsc_line})

ENDIF()

#
# Query version information:
#

IF (NOT EXISTS  ${CTEST_SOURCE_DIRECTORY}/.git)
    MESSAGE(FATAL_ERROR "Could not find .git directory")
ENDIF()

FIND_PACKAGE(Git)

IF(NOT GIT_FOUND)
    MESSAGE(FATAL_ERROR "Could not find git.")
ENDIF()




#Git_WC_INFO(${CTEST_SOURCE_DIRECTORY} bla)

EXECUTE_PROCESS(
   COMMAND ${GIT_EXECUTABLE} log -n 1 --pretty=format:"%H %h %ae"
   WORKING_DIRECTORY ${CTEST_SOURCE_DIRECTORY}
   OUTPUT_VARIABLE _git_WC_INFO
   RESULT_VARIABLE _result
   OUTPUT_STRIP_TRAILING_WHITESPACE
   )

IF(NOT ${_result} EQUAL 0)
    MESSAGE(FATAL_ERROR "Could not get git revision.")
ENDIF()

STRING(REGEX REPLACE "^\"([^ ]+) ([^ ]+) ([^\"]+)\""
         "\\1" _git_WC_REV "${_git_WC_INFO}")

STRING(REGEX REPLACE "^\"([^ ]+) ([^ ]+) ([^\"]+)\""
         "\\2" _git_WC_SHORTREV "${_git_WC_INFO}")

STRING(REGEX REPLACE "^\"([^ ]+) ([^ ]+) ([^\"]+)\""
         "\\3" _git_WC_AUTHOR "${_git_WC_INFO}")

EXECUTE_PROCESS(
   COMMAND ${GIT_EXECUTABLE} symbolic-ref HEAD
   WORKING_DIRECTORY ${CTEST_SOURCE_DIRECTORY}
   OUTPUT_VARIABLE _git_WC_BRANCH
   RESULT_VARIABLE _result
   OUTPUT_STRIP_TRAILING_WHITESPACE
   )

STRING(REGEX REPLACE "refs/heads/" ""
  _git_WC_BRANCH "${_git_WC_BRANCH}")

IF(USE_PETSC)
  SET(CTEST_BUILD_NAME "${CTEST_BUILD_NAME}-PETSc")
ENDIF()

SET(CTEST_BUILD_NAME "${CTEST_BUILD_NAME}-${_git_WC_BRANCH}")



#
# Write revision log:
#

FILE(WRITE ${CTEST_BINARY_DIRECTORY}/revision.log
"###
#
#  Git information:
#    Branch: ${_git_WC_BRANCH}
#    Commit: ${_git_WC_REV}
#    Author: ${_git_WC_AUTHOR}
#
###"
  )

#
# Append config file name to CTEST_BUILD_NAME:
#

IF(NOT "${CONFIG_FILE}" STREQUAL "")
  GET_FILENAME_COMPONENT(_conf ${CONFIG_FILE} NAME_WE)
  STRING(REGEX REPLACE "#.*$" "" _conf ${_conf})
  SET(CTEST_BUILD_NAME "${CTEST_BUILD_NAME}-${_conf}")
ENDIF()

#
# Append DESCRIPTION string to CTEST_BUILD_NAME:
#

IF(NOT "${DESCRIPTION}" STREQUAL "")
  SET(CTEST_BUILD_NAME "${CTEST_BUILD_NAME}-${DESCRIPTION}")
ENDIF()

MESSAGE("-- CTEST_BUILD_NAME:       ${CTEST_BUILD_NAME}")

#
# Declare files that should be submitted as notes:
#

SET(CTEST_NOTES_FILES
  ${CTEST_BINARY_DIRECTORY}/revision.log
  ${CTEST_BINARY_DIRECTORY}/detailed.log
  )

MESSAGE("-- CMake Options:          ${_options}")

IF(NOT "${MAKEOPTS}" STREQUAL "")
  MESSAGE("-- MAKEOPTS:               ${MAKEOPTS}")
ENDIF()


########################################################################
#                                                                      #
#                          Run the testsuite:                          #
#                                                                      #
########################################################################

CTEST_START(Experimental TRACK ${TRACK})

MESSAGE("-- Running CTEST_CONFIGURE()")
CTEST_CONFIGURE(OPTIONS "${_options}" RETURN_VALUE _res)

IF("${_res}" STREQUAL "0")
  # Only run the build stage if configure was successful:

  MESSAGE("-- Running CTEST_BUILD()")
  CTEST_BUILD(TARGET ${MAKEOPTS} NUMBER_ERRORS _res)

  IF("${_res}" STREQUAL "0")
    # Only run tests if the build was successful:

    MESSAGE("-- Running CTEST_TESTS()")
    CTEST_TEST()
  ENDIF()
ENDIF()

#
# Inject compiler information and svn revision into xml files:
#

FILE(STRINGS ${CTEST_BINARY_DIRECTORY}/Testing/TAG _tag LIMIT_COUNT 1)
SET(_path "${CTEST_BINARY_DIRECTORY}/Testing/${_tag}")
IF(NOT EXISTS ${_path})
  MESSAGE(FATAL_ERROR "
Unable to determine test submission files from TAG. Bailing out.
"
    )
ENDIF()

IF(CMAKE_SYSTEM_NAME MATCHES "Linux")
  #
  # Only use the following sed command on GNU userlands:
  #
  # TODO: Come up with a more robust way to inject this that also works on
  # BSD and Mac
  #
  FILE(GLOB _xml_files ${_path}/*.xml)
  EXECUTE_PROCESS(COMMAND sed -i -e
    s/CompilerName=\\"\\"/CompilerName=\\"${_compiler_name}\\"\\n\\tCompilerVersion=\\"${_compiler_version}\\"/g
    ${_xml_files}
    OUTPUT_QUIET RESULT_VARIABLE  _res
    )
  IF(NOT "${_res}" STREQUAL "0")
    MESSAGE(FATAL_ERROR "
  \"sed\" failed. Bailing out.
  "
      )
  ENDIF()
ENDIF()

FILE(WRITE ${_path}/Update.xml
"<?xml version=\"1.0\" encoding=\"UTF-8\"?>
<Update mode=\"Client\" Generator=\"ctest-${CTEST_VERSION}\">
<Site>${CTEST_SITE}</Site>
<BuildName>${CTEST_BUILD_NAME}</BuildName>
<BuildStamp>${_tag}-${TRACK}</BuildStamp>
<UpdateType>GIT</UpdateType>
<Revision>${_git_WC_SHORTREV}</Revision>
<Path>${_git_WC_BRANCH}</Path>
</Update>"
  )

IF("${submit}")
MESSAGE("-- Running CTEST_SUBMIT()")
CTEST_SUBMIT(RETURN_VALUE _res)

IF("${_res}" STREQUAL "0")
  MESSAGE("-- Submission successful. Goodbye!")
ENDIF()
ELSE()
MESSAGE("-- Submission skipped. Run with submit=ON")
ENDIF()
