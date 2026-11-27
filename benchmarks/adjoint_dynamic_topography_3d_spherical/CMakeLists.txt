# Copyright (C) 2011 - 2024 by the authors of the ASPECT code.
#
# This file is part of ASPECT.
#
# ASPECT is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

cmake_minimum_required(VERSION 3.13.4)

find_package(Aspect 2.4.0 QUIET HINTS ${Aspect_DIR} ../ ../../ $ENV{ASPECT_DIR})

if (NOT Aspect_FOUND)
  message(FATAL_ERROR "\n"
    "Could not find a valid ASPECT build/installation directory. "
    "Please specify the directory where you are building ASPECT by passing\n"
    "   -D Aspect_DIR=<path to ASPECT>\n"
    "to cmake or by setting the environment variable ASPECT_DIR in your shell "
    "before calling cmake.")
endif ()

DEAL_II_INITIALIZE_CACHED_VARIABLES()

set(TARGET "adjoint_dynamic_topography_3d_spherical")
project(${TARGET})

add_library(${TARGET} SHARED cellwise_perturbation.cc)
ASPECT_SETUP_PLUGIN(${TARGET})
