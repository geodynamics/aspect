# Copyright (C) 2020 by the authors of the ASPECT code.
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
# Try to find the NETCDF C library
#
# This module exports
#
#   NETCDF_LIBRARIES
#   NETCDF_INCLUDE_DIRS
#

# TODO:
# - we should check the macros NC_HAS_PARALLEL4 (parallel hdf5 support) and expose that information

SET(NETCDF_DIR "" CACHE PATH "An optional hint to a NETCDF installation")
SET_IF_EMPTY(NETCDF_DIR "$ENV{NETCDF_DIR}")

FIND_PATH(NETCDF_INCLUDE_DIR
  NAMES netcdf_meta.h
  HINTS ${NETCDF_DIR}
  PATH_SUFFIXES netcdf include
        )

FIND_LIBRARY(NETCDF_LIBRARY
  NAMES netcdf
  HINTS ${NETCDF_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )


SET(_header "${NETCDF_INCLUDE_DIR}/netcdf_meta.h")

SET(NETCDF_VERSION "unknown")

IF(EXISTS ${_header} )
  FILE(STRINGS "${_header}" _version_line
    REGEX "#define.*NC_VERSION"
    )
  STRING(REGEX REPLACE ".*\"(.+)\"" "\\1" NETCDF_VERSION
    "${_version_line}"
    )
ENDIF()


IF(NETCDF_INCLUDE_DIR AND NETCDF_LIBRARY)
  SET(NETCDF_FOUND TRUE)
  SET(NETCDF_LIBRARIES ${NETCDF_LIBRARY})
  SET(NETCDF_INCLUDE_DIRS ${NETCDF_INCLUDE_DIR})
ELSE()
    SET(NETCDF_FOUND FALSE)
ENDIF()
