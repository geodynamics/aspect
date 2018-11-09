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

# Use the environment variables if they are provided during cmake
SET(LIBDAP_LIB "" CACHE PATH "An optional hint to a libdap directory")
SET_IF_EMPTY(LIBDAP_LIB "$ENV{LIBDAP_LIB}")

SET(LIBDAP_DIR "" CACHE PATH "An optional hint to a libdap directory")
SET_IF_EMPTY(LIBDAP_DIR "$ENV{LIBDAP_DIR}")

# Clear the cache of previous values
UNSET (LIBDAP_INCLUDE_DIR CACHE)
UNSET (LIBDAP_LIBRARY CACHE)
UNSET (LIBDAP_CLIENT_LIBRARY CACHE)


FIND_PATH(LIBDAP_INCLUDE_DIR 
  NAMES Connect.h
  HINTS ${LIBDAP_DIR}
  )

FIND_LIBRARY(LIBDAP_LIBRARY 
  NAMES libdap.a 
  HINTS ${LIBDAP_LIB}
  )
  
FIND_LIBRARY(LIBDAP_CLIENT_LIBRARY
  NAMES libdapclient.a
  HINTS ${LIBDAP_LIB}
  )

SET(LIBDAP_FOUND TRUE)
SET(LIBDAP_INCLUDE_DIRS ${LIBDAP_INCLUDE_DIR})
SET(LIBDAP_LIBRARIES ${LIBDAP_LIBRARY} ${LIBDAP_CLIENT_LIBRARY})
