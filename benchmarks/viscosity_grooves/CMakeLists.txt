# Copyright (C) 2015 - 2024 by the authors of the ASPECT code.
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

set(TARGET "viscosity_grooves")
project(${TARGET})

add_library(${TARGET} SHARED viscosity_grooves.cc)
ASPECT_SETUP_PLUGIN(${TARGET})
