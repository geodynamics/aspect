# ------------------------------------------------------
#
# SPDX-License-Identifier: GPL-2.0-or-later
#
# SPDX-FileCopyrightText: Copyright (C) 2023-2024 by the ASPECT authors.
#
# This file is part of ASPECT.
#
# Detailed license information governing the source code and
# contributions can be found in the folder LICENSES and in
# LICENSE and CONTRIBUTING.md at the top level directory.
#
# ------------------------------------------------------

CMAKE_MINIMUM_REQUIRED(VERSION 3.13.4)

FIND_PACKAGE(Aspect REQUIRED HINTS ${Aspect_DIR} ../ ../../ $ENV{ASPECT_DIR})
DEAL_II_INITIALIZE_CACHED_VARIABLES()

SET(TARGET "solubility")
PROJECT(${TARGET})

ADD_LIBRARY(${TARGET} SHARED solubility.cc)
ASPECT_SETUP_PLUGIN(${TARGET})
