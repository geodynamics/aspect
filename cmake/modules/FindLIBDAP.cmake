## ---------------------------------------------------------------------
##
## Copyright (C) 2018 by Kodi Neumiller, James Gallagher
##
## ---------------------------------------------------------------------

# Try to find the Libdap library

# Use the environment variables if they are provided during cmake
SET(LIBDAP_LIB "" CACHE PATH "An optional hint to the libdap lib directory")
SET_IF_EMPTY(LIBDAP_LIB "$ENV{LIBDAP_LIB}")

SET(LIBDAP_INC "" CACHE PATH "An optional hint to the libdap include directory")
SET_IF_EMPTY(LIBDAP_DIR "$ENV{LIBDAP_INC}")

SET(LIBDAP_DIR "" CACHE PATH "An optional hint to the libdap installation directory")
SET_IF_EMPTY(LIBDAP_DIR "$ENV{LIBDAP_DIR}")

# Clear the cache of previous values
UNSET (LIBDAP_INCLUDE_DIR CACHE)
UNSET (LIBDAP_LIBRARY CACHE)
UNSET (LIBDAP_CLIENT_LIBRARY CACHE)


FIND_PATH(LIBDAP_INCLUDE_DIR 
  NAMES Connect.h
  HINTS ${LIBDAP_INC} ${LIBDAP_DIR}/include/libdap 
  )

FIND_LIBRARY(LIBDAP_LIBRARY 
  NAMES libdap.a 
  HINTS ${LIBDAP_LIB} ${LIBDAP_DIR}/lib
  )
  
FIND_LIBRARY(LIBDAP_CLIENT_LIBRARY
  NAMES libdapclient.a
  HINTS ${LIBDAP_LIB} ${LIBDAP_DIR}/lib
  )

SET(LIBDAP_FOUND TRUE)
SET(LIBDAP_INCLUDE_DIRS ${LIBDAP_INCLUDE_DIR})
SET(LIBDAP_LIBRARIES ${LIBDAP_CLIENT_LIBRARY} ${LIBDAP_LIBRARY})
