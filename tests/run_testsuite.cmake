#
# Copyright (C) 2013 - 2024 by Matthias Maier
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
# along with ASPECT; see the file LICENSE.  If not see
# <http://www.gnu.org/licenses/>.
#


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

cmake_minimum_required(VERSION 3.13.4)
message("-- This is CTest ${CMAKE_VERSION}")

#
# TRACK: Default to Experimental:
#

if("${TRACK}" STREQUAL "")
  set(TRACK "Experimental")
endif()

if( NOT "${TRACK}" STREQUAL "Experimental"
    AND NOT "${TRACK}" STREQUAL "Build Tests"
    AND NOT "${TRACK}" STREQUAL "Nightly"
    AND NOT "${TRACK}" STREQUAL "Regression Tests" )
  message(FATAL_ERROR "
Unknown TRACK \"${TRACK}\" - see the manual for valid values.
"
    )
endif()

message("-- TRACK:                  ${TRACK}")

#
# CTEST_SOURCE_DIRECTORY:
#

if("${CTEST_SOURCE_DIRECTORY}" STREQUAL "")
  #
  # If CTEST_SOURCE_DIRECTORY is not set we just assume that this script
  # resides in the top level source directory
  #
  set(CTEST_SOURCE_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})

  # The code above set the source directory of that of the run_testsuite.cmake
  # script, but we need the directory of aspect, which is simply one level
  # higher
  IF ("${CTEST_SOURCE_DIRECTORY}" MATCHES "/tests")
    set(CTEST_SOURCE_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/../)
  endif()

  if(NOT EXISTS ${CTEST_SOURCE_DIRECTORY}/CMakeLists.txt)
    message(FATAL_ERROR "
Could not find a suitable source directory. Please, set
CTEST_SOURCE_DIRECTORY manually to the appropriate source directory.
"
      )
  endif()
endif()

message("-- CTEST_SOURCE_DIRECTORY: ${CTEST_SOURCE_DIRECTORY}")

#
# Read in custom config files:
#

ctest_read_custom_files(${CTEST_SOURCE_DIRECTORY})

#
# CTEST_BINARY_DIRECTORY:
#

if("${CTEST_BINARY_DIRECTORY}" STREQUAL "")
  #
  # If CTEST_BINARY_DIRECTORY is not set we just use the current directory
  # except if it is equal to CTEST_SOURCE_DIRECTORY in which case we fail.
  #
  set(CTEST_BINARY_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endif()

#
# Read in custom config files:
#

ctest_read_custom_files(${CTEST_BINARY_DIRECTORY})

# Make sure that for a build test the directory is empty:
file(GLOB _test ${CTEST_BINARY_DIRECTORY}/*)
if( "${TRACK}" STREQUAL "Build Tests"
    AND NOT "${_test}" STREQUAL "" )
      message(FATAL_ERROR "
TRACK was set to \"Build Tests\" which require an empty build directory.
But files were found in \"${CTEST_BINARY_DIRECTORY}\"
"
        )
endif()

message("-- CTEST_BINARY_DIRECTORY: ${CTEST_BINARY_DIRECTORY}")

#
# CTEST_CMAKE_GENERATOR:
#

# Query Generator from build directory (if possible):
if(EXISTS ${CTEST_BINARY_DIRECTORY}/CMakeCache.txt)
  file(STRINGS ${CTEST_BINARY_DIRECTORY}/CMakeCache.txt _generator
    REGEX "^CMAKE_GENERATOR:"
    )
  string(REGEX REPLACE "^.*=" "" _generator ${_generator})
endif()

if("${CTEST_CMAKE_GENERATOR}" STREQUAL "")
  if(NOT "${_generator}" STREQUAL "")
    set(CTEST_CMAKE_GENERATOR ${_generator})
  else()
    # default to "Unix Makefiles"
    set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
  endif()
else()
  # ensure that CTEST_CMAKE_GENERATOR (that was apparently set) is
  # compatible with the build directory:
  if( NOT "${CTEST_CMAKE_GENERATOR}" STREQUAL "${_generator}"
      AND NOT "${_generator}" STREQUAL "" )
    message(FATAL_ERROR "
The build directory is already set up with Generator \"${_generator}\", but
CTEST_CMAKE_GENERATOR was set to a different Generator \"${CTEST_CMAKE_GENERATOR}\".
"
     )
  endif()
endif()

message("-- CTEST_CMAKE_GENERATOR:  ${CTEST_CMAKE_GENERATOR}")

#
# CTEST_SITE:
#

find_program(HOSTNAME_COMMAND NAMES hostname)
if(NOT "${HOSTNAME_COMMAND}" MATCHES "-NOTFOUND")
  exec_program(${HOSTNAME_COMMAND} OUTPUT_VARIABLE _hostname)
  set(CTEST_SITE "${_hostname}")
else()
  # Well, no hostname available. What about:
  set(CTEST_SITE "BobMorane")
endif()

message("-- CTEST_SITE:             ${CTEST_SITE}")

if( "${TRACK}" STREQUAL "Regression Tests"
    AND NOT CTEST_SITE MATCHES "tester" )
  message(FATAL_ERROR "
I'm sorry ${CTEST_SITE}, I'm afraid I can't do that.
The TRACK \"Regression Tests\" is not for you.
"
    )
endif()

#
# Assemble configuration options, we need it now:
#

if("${MAKEOPTS}" STREQUAL "")
  set(MAKEOPTS $ENV{MAKEOPTS})
endif()

if(NOT "${CONFIG_FILE}" STREQUAL "")
  set(_options "-C${CONFIG_FILE}")
endif()

if("${TRACK}" STREQUAL "Build Tests")
  set(TEST_PICKUP_REGEX "^build_tests")
endif()

# Pass all relevant variables down to configure:
get_cmake_property(_variables VARIABLES)

#
# CTEST_BUILD_NAME:
#

# Append compiler information to CTEST_BUILD_NAME:
if(NOT EXISTS ${CTEST_BINARY_DIRECTORY}/detailed.log)
  message(FATAL_ERROR "could not find detailed.log")
endif()

if(EXISTS ${CTEST_BINARY_DIRECTORY}/detailed.log)
  file(STRINGS ${CTEST_BINARY_DIRECTORY}/detailed.log _compiler_id
    REGEX "CMAKE_CXX_COMPILER:"
    )
  string(REGEX REPLACE
    "^.*CMAKE_CXX_COMPILER:     \(.*\) on platform.*$" "\\1"
    _compiler_id ${_compiler_id}
    )
  string(REGEX REPLACE "^\(.*\) .*$" "\\1" _compiler_name ${_compiler_id})
  string(REGEX REPLACE "^.* " "" _compiler_version ${_compiler_id})
  string(REGEX REPLACE " " "-" _compiler_id ${_compiler_id})
  if( NOT "${_compiler_id}" STREQUAL "" OR
      _compiler_id MATCHES "CMAKE_CXX_COMPILER" )
    set(CTEST_BUILD_NAME "${_compiler_id}")
  endif()

endif()

#
# Query version information:
#

IF (NOT EXISTS  ${CTEST_SOURCE_DIRECTORY}/.git)
    message(FATAL_ERROR "Could not find .git directory")
endif()

find_package(Git)

if(NOT GIT_FOUND)
    message(FATAL_ERROR "Could not find git.")
endif()




#Git_wc_info(${CTEST_SOURCE_DIRECTORY} bla)

execute_process(
   COMMAND ${GIT_EXECUTABLE} log -n 1 --pretty=format:"%H %h %ae"
   WORKING_DIRECTORY ${CTEST_SOURCE_DIRECTORY}
   OUTPUT_VARIABLE _git_WC_INFO
   RESULT_VARIABLE _result
   OUTPUT_STRIP_TRAILING_WHITESPACE
   )

if(NOT ${_result} EQUAL 0)
    message(FATAL_ERROR "Could not get git revision.")
endif()

string(REGEX REPLACE "^\"([^ ]+) ([^ ]+) ([^\"]+)\""
         "\\1" _git_WC_REV "${_git_WC_INFO}")

string(REGEX REPLACE "^\"([^ ]+) ([^ ]+) ([^\"]+)\""
         "\\2" _git_WC_SHORTREV "${_git_WC_INFO}")

string(REGEX REPLACE "^\"([^ ]+) ([^ ]+) ([^\"]+)\""
         "\\3" _git_WC_AUTHOR "${_git_WC_INFO}")

execute_process(
   COMMAND ${GIT_EXECUTABLE} symbolic-ref HEAD
   WORKING_DIRECTORY ${CTEST_SOURCE_DIRECTORY}
   OUTPUT_VARIABLE _git_WC_BRANCH
   RESULT_VARIABLE _result
   OUTPUT_STRIP_TRAILING_WHITESPACE
   )

string(REGEX REPLACE "refs/heads/" ""
  _git_WC_BRANCH "${_git_WC_BRANCH}")

set(CTEST_BUILD_NAME "${CTEST_BUILD_NAME}-${_git_WC_BRANCH}")



#
# Write revision log:
#

file(WRITE ${CTEST_BINARY_DIRECTORY}/revision.log
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

if(NOT "${CONFIG_FILE}" STREQUAL "")
  get_filename_component(_conf ${CONFIG_FILE} NAME_WE)
  string(REGEX REPLACE "#.*$" "" _conf ${_conf})
  set(CTEST_BUILD_NAME "${CTEST_BUILD_NAME}-${_conf}")
endif()

#
# Append DESCRIPTION string to CTEST_BUILD_NAME:
#

if(NOT "${DESCRIPTION}" STREQUAL "")
  set(CTEST_BUILD_NAME "${CTEST_BUILD_NAME}-${DESCRIPTION}")
endif()

message("-- CTEST_BUILD_NAME:       ${CTEST_BUILD_NAME}")

#
# Declare files that should be submitted as notes:
#

set(CTEST_NOTES_FILES
  ${CTEST_BINARY_DIRECTORY}/revision.log
  ${CTEST_BINARY_DIRECTORY}/detailed.log
  )

message("-- CMake Options:          ${_options}")

if(NOT "${MAKEOPTS}" STREQUAL "")
  message("-- MAKEOPTS:               ${MAKEOPTS}")
endif()


########################################################################
#                                                                      #
#                          Run the testsuite:                          #
#                                                                      #
########################################################################

ctest_start(Experimental TRACK ${TRACK})

message("-- Running ctest_configure()")
ctest_configure(OPTIONS "${_options}" RETURN_VALUE _res)

if("${_res}" STREQUAL "0")
  # Only run the build stage if configure was successful:

  message("-- Running ctest_build()")
  ctest_build(TARGET ${MAKEOPTS} NUMBER_ERRORS _res)

  if("${_res}" STREQUAL "0")
    # Only run tests if the build was successful:

    message("-- Running ctest_tests()")
    ctest_test()
  endif()
endif()

#
# Inject compiler information and svn revision into xml files:
#

file(STRINGS ${CTEST_BINARY_DIRECTORY}/Testing/TAG _tag LIMIT_COUNT 1)
set(_path "${CTEST_BINARY_DIRECTORY}/Testing/${_tag}")
if(NOT EXISTS ${_path})
  message(FATAL_ERROR "
Unable to determine test submission files from TAG. Bailing out.
"
    )
endif()

if(CMAKE_SYSTEM_NAME MATCHES "Linux")
  #
  # Only use the following sed command on GNU userlands:
  #
  # TODO: Come up with a more robust way to inject this that also works on
  # BSD and Mac
  #
  file(GLOB _xml_files ${_path}/*.xml)
  execute_process(COMMAND sed -i -e
    s/CompilerName=\\"\\"/CompilerName=\\"${_compiler_name}\\"\\n\\tCompilerVersion=\\"${_compiler_version}\\"/g
    ${_xml_files}
    OUTPUT_QUIET RESULT_VARIABLE  _res
    )
  if(NOT "${_res}" STREQUAL "0")
    message(FATAL_ERROR "
  \"sed\" failed. Bailing out.
  "
      )
  endif()
endif()

file(WRITE ${_path}/Update.xml
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

if("${submit}")
message("-- Running ctest_submit()")
ctest_submit(RETURN_VALUE _res)

if("${_res}" STREQUAL "0")
  message("-- Submission successful. Goodbye!")
endif()
else()
message("-- Submission skipped. Run with submit=ON")
endif()
