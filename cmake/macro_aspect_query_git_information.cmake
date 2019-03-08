# Copyright (C) 2018 by the authors of the ASPECT code.
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
#
# This file implements the ASPECT_QUERY_GIT_INFORMATION macro, which is
# part of ASPECT. It is a modified copy of the corresponding deal.II
# macro DEAL_II_QUERY_GIT_INFORMATION.
#
# Usage:
#       ASPECT_QUERY_GIT_INFORMATION()
#       ASPECT_QUERY_GIT_INFORMATION("CUSTOM_PREFIX")
#
# This will try to gather information about current branch, as well as
# short and long revision. If ${CMAKE_SOURCE_DIR} is the root of a git
# repository the following variables will be populated:
#
#       GIT_BRANCH
#       GIT_REVISION
#       GIT_SHORTREV
#
# The macro can be called with an optional PREFIX argument to prefix the
# variables:
#
#       PREFIX_GIT_BRANCH
#       PREFIX_GIT_REVISION
#       PREFIX_GIT_SHORTREV
#

MACRO(ASPECT_QUERY_GIT_INFORMATION)

  MESSAGE(STATUS "Query git repository information.")

  # Set prefix.
  SET(_prefix "")
  IF(NOT "${ARGN}" STREQUAL "")
    SET(_prefix "${ARGN}_")
  ENDIF()

  FIND_PACKAGE(Git)

  #
  # Only run the following if we have git and the source directory seems to
  # be under version control.
  #
  IF(GIT_FOUND AND EXISTS ${CMAKE_SOURCE_DIR}/.git/HEAD)
    #
    # Bogus configure_file calls to trigger a reconfigure, and thus an
    # update of branch and commit information every time HEAD has changed.
    #
    CONFIGURE_FILE(
      ${CMAKE_SOURCE_DIR}/.git/HEAD
      ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/HEAD
      )
    FILE(STRINGS ${CMAKE_SOURCE_DIR}/.git/HEAD _head_ref LIMIT_COUNT 1)
    STRING(REPLACE "ref: " "" _head_ref ${_head_ref})
    IF(EXISTS ${CMAKE_SOURCE_DIR}/.git/${_head_ref})
      CONFIGURE_FILE(
        ${CMAKE_SOURCE_DIR}/.git/${_head_ref}
        ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/HEAD_REF
        )
    ENDIF()

    #
    # Query for revision:
    #

    EXECUTE_PROCESS(
       COMMAND ${GIT_EXECUTABLE} log -n 1 --pretty=format:"%H %h"
       WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
       OUTPUT_VARIABLE _info
       RESULT_VARIABLE _result
       OUTPUT_STRIP_TRAILING_WHITESPACE
       )
    IF(${_result} EQUAL 0)
      STRING(REGEX REPLACE "^\"([^ ]+) ([^ ]+)\"$"
        "\\1" ${_prefix}GIT_REVISION "${_info}")
      STRING(REGEX REPLACE "^\"([^ ]+) ([^ ]+)\"$"
        "\\2" ${_prefix}GIT_SHORTREV "${_info}")
    ENDIF()

    #
    # Query for branch:
    #

    EXECUTE_PROCESS(
       COMMAND ${GIT_EXECUTABLE} symbolic-ref HEAD
       WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
       OUTPUT_VARIABLE _branch
       RESULT_VARIABLE _result
       OUTPUT_STRIP_TRAILING_WHITESPACE
       )
    IF(${_result} EQUAL 0)
      STRING(REGEX REPLACE "refs/heads/" "" ${_prefix}GIT_BRANCH "${_branch}")
    ENDIF()
  ENDIF()

ENDMACRO()
