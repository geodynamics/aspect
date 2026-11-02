# ------------------------------------------------------
#
# SPDX-License-Identifier: GPL-2.0-or-later
#
# SPDX-FileCopyrightText: Copyright (C) 2018-2025 by the ASPECT authors.
#
# This file is part of ASPECT.
#
# Detailed license information governing the source code and
# contributions can be found in the folder LICENSES and in
# LICENSE and CONTRIBUTING.md at the top level directory.
#
# ------------------------------------------------------

################### top matter ################
cmake_minimum_required(VERSION 3.13.4)

project(unittests CXX)

set(ASPECT_BINARY "${ASPECT_BINARY}" CACHE STRING "" FORCE)

enable_testing()

message(STATUS "creating unit_tests ctest in ${CMAKE_CURRENT_BINARY_DIR} binary=${ASPECT_BINARY}")

add_test(NAME unit_tests
         COMMAND ${ASPECT_BINARY} --test
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      )
