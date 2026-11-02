# ------------------------------------------------------
#
# SPDX-License-Identifier: GPL-2.0-or-later
#
# SPDX-FileCopyrightText: Copyright (C) 2017-2025 by the ASPECT authors.
#
# This file is part of ASPECT.
#
# Detailed license information governing the source code and
# contributions can be found in the folder LICENSES and in
# LICENSE and CONTRIBUTING.md at the top level directory.
#
# ------------------------------------------------------

cmake_minimum_required(VERSION 3.13.4)

find_package(Aspect 2.4.0 QUIET HINTS ${Aspect_DIR} ../ ../../ $ENV{ASPECT_DIR})

if (NOT Aspect_FOUND)
  message(FATAL_ERROR "\n"
    "Could not find a valid ASPECT build/installation directory. "
    "Please specify the directory where you are building ASPECT by passing\n"
    "   -D Aspect_DIR=<path to ASPECT>\n"
    "to cmake or by setting the environment variable ASPECT_DIR in your shell "
    "before calling cmake. See the section 'How to write a plugin' in the "
    "manual for more information.")
endif ()

DEAL_II_INITIALIZE_CACHED_VARIABLES()

set(TARGET "inclusion_compositional_fields")
project(${TARGET})

add_library(${TARGET} SHARED inclusion_compositional_fields.cc)
ASPECT_SETUP_PLUGIN(${TARGET})
