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
MESSAGE(STATUS "libdap include dir = ${LIBDAP_INCLUDE_DIR}")

# Get the remaining libraries necesary for libdap to run.
# Do this by executing a command line command and storing the result
EXECUTE_PROCESS(COMMAND dap-config --libs OUTPUT_VARIABLE result)
MESSAGE(STATUS "testing the dap-config ----- ${result}")

FIND_LIBRARY(LIBDAP_LIBRARY 
  NAMES libdap.a
  HINTS ${LIBDAP_DIR}
  PATHS /Users/kodi/src/hyrax/build/lib64 /Users/kodi/src/hyrax/build/lib
  )
MESSAGE(STATUS "libdap lib = ${LIBDAP_LIBRARY}")

#FIND_LIBRARY(LIBXML_INCLUDE_DIR
#  NAMES xlink.h
#  PATHS /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.14.sdk/usr/include/libxml2/libxml
#  )
#MESSAGE(STATUS "libxml lib = ${LIBXML_INCLUDE_DIR}")

#SET(TESTVAR -L/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.13.sdk/usr/lib -lxml2 -lz -lpthread -licucore -lm -lpthread)
#MESSAGE(STATUS "TESTVAR = ${TESTVAR}")

SET(LIBDAP_FOUND TRUE)
SET(LIBDAP_INCLUDE_DIRS ${LIBDAP_INCLUDE_DIR})
SET(LIBDAP_LIBRARIES ${LIBDAP_LIBRARY})
