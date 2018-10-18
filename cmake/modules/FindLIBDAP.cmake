## ---------------------------------------------------------------------
##
## Copyright (C) 2014 - 2015 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE.md at
## the top level directory of deal.II.
##
## ---------------------------------------------------------------------

# Try to find the Libdap library
MESSAGE(STATUS "Testing the FindLIBDAP.cmake")
SET(LIBDAP_DIR "" CACHE PATH "An optional hint to a libdap directory")
SET_IF_EMPTY(LIBDAP_DIR "$ENV{LIBDAP_DIR}")

# Clear the cache of previous values
UNSET (LIBDAP_INCLUDE_DIR CACHE)
UNSET (LIBDAP_LIBRARY CACHE)


FIND_PATH(LIBDAP_INCLUDE_DIR 
  NAMES Connect.h
  HINTS ${LIBDAP_DIR}
  PATHS /Users/kodi/src/hyrax/build/include/libdap
  )

FIND_LIBRARY(LIBDAP_LIBRARY 
  NAMES libdap.a
  HINTS ${LIBDAP_DIR}
  PATHS /Users/kodi/src/hyrax/build/lib64 /Users/kodi/src/hyrax/build/lib
  )

SET(LIBDAP_FOUND TRUE)
SET(LIBDAP_INCLUDE_DIRS ${LIBDAP_INCLUDE_DIR})
SET(LIBDAP_LIBRARIES ${LIBDAP_LIBRARY})
