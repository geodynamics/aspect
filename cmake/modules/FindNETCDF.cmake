# ------------------------------------------------------
#
# SPDX-License-Identifier: GPL-2.0-or-later
#
# SPDX-FileCopyrightText: Copyright (C) 2020-2024 by the ASPECT authors.
#
# This file is part of ASPECT.
#
# Detailed license information governing the source code and
# contributions can be found in the folder LICENSES and in
# LICENSE and CONTRIBUTING.md at the top level directory.
#
# ------------------------------------------------------

#
# Try to find the NETCDF C library
#
# This module exports
#
#   NETCDF_LIBRARIES
#   NETCDF_INCLUDE_DIRS
#

# TODO:
# - we should check the macros NC_HAS_PARALLEL4 (parallel hdf5 support) and expose that information

set(NETCDF_DIR "" CACHE PATH "An optional hint to a NETCDF installation")
set_if_empty(NETCDF_DIR "$ENV{NETCDF_DIR}")

find_path(NETCDF_INCLUDE_DIR
  NAMES netcdf_meta.h
  HINTS ${NETCDF_DIR}
  PATH_SUFFIXES netcdf include
        )

find_library(NETCDF_LIBRARY
  NAMES netcdf
  HINTS ${NETCDF_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )


set(_header "${NETCDF_INCLUDE_DIR}/netcdf_meta.h")

set(NETCDF_VERSION "unknown")

if(EXISTS ${_header} )
  file(STRINGS "${_header}" _version_line
    REGEX "#define.*NC_VERSION"
    )
  string(REGEX REPLACE ".*\"(.+)\"" "\\1" NETCDF_VERSION
    "${_version_line}"
    )
endif()


if(NETCDF_INCLUDE_DIR AND NETCDF_LIBRARY)
  set(NETCDF_FOUND TRUE)
  set(NETCDF_LIBRARIES ${NETCDF_LIBRARY})
  set(NETCDF_INCLUDE_DIRS ${NETCDF_INCLUDE_DIR})
else()
    set(NETCDF_FOUND FALSE)
endif()
