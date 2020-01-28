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
        NAMES libdap.so libdap.dylib
        HINTS ${LIBDAP_LIB} ${LIBDAP_DIR}/lib
        )

FIND_LIBRARY(LIBDAP_CLIENT_LIBRARY
        NAMES libdapclient.so libdapclient.dylib
        HINTS ${LIBDAP_LIB} ${LIBDAP_DIR}/lib
        )

FIND_PROGRAM(LIBDAP_CONFIG_EXECUTABLE
        NAMES dap-config
        HINTS ${LIBDAP_DIR}/bin)

#Lookup and add the CURL and libxml2 libraries
IF(LIBDAP_CONFIG_EXECUTABLE)
    EXECUTE_PROCESS(COMMAND ${LIBDAP_CONFIG_EXECUTABLE} --libs
            OUTPUT_VARIABLE _libs
            ERROR_QUIET
            OUTPUT_STRIP_TRAILING_WHITESPACE)
    STRING(REPLACE " " ";" _libs "${_libs}")
    SET(_path "")
    FOREACH(_lib ${_libs})
        IF (${_lib} MATCHES "^-L")
            STRING(SUBSTRING ${_lib} 2 -1 _path)
        ENDIF()
        IF (${_lib} MATCHES "^-lcurl")
            SET(LIBCURL_PATH "${_path}")
            SET(LIBCURL_FOUND TRUE)
        ENDIF()
        IF (${_lib} MATCHES "^-lxml2")
            SET(LIBXML2_PATH "${_path}")
            SET(LIBXML2_FOUND TRUE)
        ENDIF()
    ENDFOREACH()

    EXECUTE_PROCESS(COMMAND ${LIBDAP_CONFIG_EXECUTABLE} --cflags
            OUTPUT_VARIABLE _flags
            ERROR_QUIET
            OUTPUT_STRIP_TRAILING_WHITESPACE)
    STRING(REPLACE " " ";" _flags "${_flags}")
    FOREACH(_flag ${_flags})
        IF (${_flag} MATCHES "^-I")
            SET(LIBDAP_INCLUDE_DIR ${LIBDAP_INCLUDE_DIR} ${_path})
        ENDIF()
    ENDFOREACH()
ENDIF()

FIND_LIBRARY(LIBCURL_LIBRARIES
        NAMES libcurl.so libcurl.dylib
        HINTS ${LIBCURL_PATH}
        )
FIND_LIBRARY(LIBXML2_LIBRARIES
        NAMES libxml2.so libxml2.dylib
        HINTS ${LIBXML2_PATH}
        )

MESSAGE(STATUS "-- LIBXML2_LIBRARIES: ${LIBXML2_LIBRARIES}")
MESSAGE(STATUS "-- LIBCURL_LIBRARIES: ${LIBCURL_LIBRARIES}")
#MESSAGE(STATUS "-- LIBDAP_INCLUDE_DIR: ${LIBDAP_INCLUDE_DIR}")


IF(LIBDAP_LIBRARY AND LIBDAP_CLIENT_LIBRARY AND LIBDAP_CONFIG_EXECUTABLE)
    SET(LIBDAP_FOUND TRUE)
    SET(LIBDAP_INCLUDE_DIRS ${LIBDAP_INCLUDE_DIR})
    SET(LIBDAP_LIBRARIES ${LIBDAP_CLIENT_LIBRARY} ${LIBDAP_LIBRARY} ${LIBXML2_LIBRARIES} ${LIBCURL_LIBRARIES})
ELSE()
    SET(LIBDAP_FOUND FALSE)
ENDIF()