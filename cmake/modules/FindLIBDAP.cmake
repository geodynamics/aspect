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

# Try to find the Libdap library

# Use the environment variables if they are provided during cmake
set(LIBDAP_LIB "" CACHE PATH "An optional hint to the libdap lib directory")
set_if_empty(LIBDAP_LIB "$ENV{LIBDAP_LIB}")

set(LIBDAP_INC "" CACHE PATH "An optional hint to the libdap include directory")
set_if_empty(LIBDAP_DIR "$ENV{LIBDAP_INC}")

set(LIBDAP_DIR "" CACHE PATH "An optional hint to the libdap installation directory")
set_if_empty(LIBDAP_DIR "$ENV{LIBDAP_DIR}")

# Clear the cache of previous values
UNSET (LIBDAP_INCLUDE_DIR CACHE)
UNSET (LIBDAP_LIBRARY CACHE)
UNSET (LIBDAP_CLIENT_LIBRARY CACHE)


find_path(LIBDAP_INCLUDE_DIR
        NAMES Connect.h
        HINTS ${LIBDAP_INC} ${LIBDAP_DIR}/include/libdap
        )

find_library(LIBDAP_LIBRARY
        NAMES libdap.so libdap.dylib
        HINTS ${LIBDAP_LIB} ${LIBDAP_DIR}/lib
        )

find_library(LIBDAP_CLIENT_LIBRARY
        NAMES libdapclient.so libdapclient.dylib
        HINTS ${LIBDAP_LIB} ${LIBDAP_DIR}/lib
        )

find_program(LIBDAP_CONFIG_EXECUTABLE
        NAMES dap-config
        HINTS ${LIBDAP_DIR}/bin)

#Lookup and add the CURL and libxml2 libraries
if(LIBDAP_CONFIG_EXECUTABLE)
    message(STATUS "probing ${LIBDAP_CONFIG_EXECUTABLE} for libdap configuration:")
    execute_process(COMMAND ${LIBDAP_CONFIG_EXECUTABLE} --libs
            OUTPUT_VARIABLE _libs
            ERROR_QUIET
            OUTPUT_STRIP_TRAILING_WHITESPACE)
    string(REPLACE " " ";" _libs "${_libs}")
    set(_path "")
    foreach(_lib ${_libs})
        IF (${_lib} MATCHES "^-L")
            string(SUBSTRING ${_lib} 2 -1 _path)
        endif()
        IF (${_lib} MATCHES "^-lcurl")
            set(LIBCURL_PATH "${_path}")
            set(LIBCURL_FOUND TRUE)
        endif()
        IF (${_lib} MATCHES "^-lxml2")
            set(LIBXML2_PATH "${_path}")
            set(LIBXML2_FOUND TRUE)
        endif()
    endforeach()

    execute_process(COMMAND ${LIBDAP_CONFIG_EXECUTABLE} --cflags
            OUTPUT_VARIABLE _flags
            ERROR_QUIET
            OUTPUT_STRIP_TRAILING_WHITESPACE)
    string(REPLACE " " ";" _flags "${_flags}")
    foreach(_flag ${_flags})
        IF (${_flag} MATCHES "^-I")
        string(SUBSTRING ${_flag} 2 -1 _path)
        set(LIBDAP_INCLUDE_DIR ${LIBDAP_INCLUDE_DIR} ${_path})
        endif()
    endforeach()
endif()

find_library(LIBCURL_LIBRARIES
        NAMES libcurl.so libcurl.dylib
        HINTS ${LIBCURL_PATH}
        )
find_library(LIBXML2_LIBRARIES
        NAMES libxml2.so libxml2.dylib
        HINTS ${LIBXML2_PATH}
        )

message(STATUS "  LIBXML2_LIBRARIES: ${LIBXML2_LIBRARIES}")
message(STATUS "  LIBCURL_LIBRARIES: ${LIBCURL_LIBRARIES}")
message(STATUS "  LIBDAP_INCLUDE_DIRS: ${LIBDAP_INCLUDE_DIR}")

if(LIBDAP_LIBRARY AND LIBDAP_CLIENT_LIBRARY AND LIBDAP_CONFIG_EXECUTABLE AND LIBCURL_LIBRARIES AND LIBXML2_LIBRARIES)
    set(LIBDAP_FOUND TRUE)
    set(LIBDAP_INCLUDE_DIRS ${LIBDAP_INCLUDE_DIR})
    set(LIBDAP_LIBRARY ${LIBDAP_LIBRARY})
    set(LIBDAP_LIBRARIES ${LIBDAP_CLIENT_LIBRARY} ${LIBDAP_LIBRARY} ${LIBXML2_LIBRARIES} ${LIBCURL_LIBRARIES})
else()
    set(LIBDAP_FOUND FALSE)
endif()
