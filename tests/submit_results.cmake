# ------------------------------------------------------
#
# SPDX-License-Identifier: GPL-2.0-or-later
#
# SPDX-FileCopyrightText: Copyright (C) 2013-2024 by the ASPECT authors.
#
# This file is part of ASPECT.
#
# Detailed license information governing the source code and
# contributions can be found in the folder LICENSES and in
# LICENSE and CONTRIBUTING.md at the top level directory.
#
# ------------------------------------------------------

#
# Submit existing test results to CDash.
#
# Usage:
#   Invoke this script in a directory with test results already present
#   under ./Testing, i.e. valid ./Testing/TAG pointing to test results:
#
#   ctest -S ../tests/submit_results.cmake
#
# You may specify CTEST_SOURCE_DIRECTORY to point to a directory containing
# CTestConfig.cmake
#

if("${CTEST_SOURCE_DIRECTORY}" STREQUAL "")
  # Get parent directory:
  get_filename_component(CTEST_SOURCE_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}" PATH)

  if(NOT EXISTS ${CTEST_SOURCE_DIRECTORY}/CTestConfig.cmake)
    set(CTEST_SOURCE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  endif()
endif()

message("-- CTEST_SOURCE_DIRECTORY: ${CTEST_SOURCE_DIRECTORY}")

set(CTEST_BINARY_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

message("-- CTEST_BINARY_DIRECTORY: ${CTEST_BINARY_DIRECTORY}")

file(STRINGS ${CTEST_BINARY_DIRECTORY}/Testing/TAG _tag)
list(GET _tag 1 _track)

if("${_track}" STREQUAL "")
  message(FATAL_ERROR "
No test results found. Bailing out.
"
    )
endif()

ctest_start(Experimental TRACK ${_track} APPEND)
ctest_submit()
