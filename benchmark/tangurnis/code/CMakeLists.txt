CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

FIND_PACKAGE(Aspect 1.5.0 QUIET HINTS ${Aspect_DIR} ../ ../../ $ENV{ASPECT_DIR})

IF (NOT Aspect_FOUND)
  MESSAGE(FATAL_ERROR "\n"
	"Could not find a valid ASPECT build/installation directory. "
	"Please specify the directory where you are building ASPECT by passing\n"
	"   -D Aspect_DIR=<path to ASPECT>\n"
	"to cmake or by setting the environment variable ASPECT_DIR in your shell "
	"before calling cmake. See the section 'How to write a plugin' in the "
        "manual for more information.")
ENDIF ()

DEAL_II_INITIALIZE_CACHED_VARIABLES()

SET(TARGET "tangurnis")
PROJECT(${TARGET})

ADD_LIBRARY(${TARGET} SHARED tangurnis.cc)
ASPECT_SETUP_PLUGIN(${TARGET})
