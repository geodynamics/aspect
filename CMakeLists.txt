MESSAGE("====================================================")
MESSAGE("============ Configuring ASPECT ====================")
MESSAGE("====================================================")

# Set the name of the project and target:
SET(TARGET "aspect")

FILE(GLOB_RECURSE TARGET_SRC  "source/*.cc")
INCLUDE_DIRECTORIES(include)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

FIND_PACKAGE(deal.II 8.0 QUIET
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()

IF(NOT DEAL_II_WITH_TRILINOS OR NOT DEAL_II_WITH_P4EST)
  MESSAGE(FATAL_ERROR
    "\nAspect requires a deal.II installation built with support for Trilinos and p4est but one or both of these appears to be missing!\n"
    )
  IF(NOT DEAL_II_WITH_TRILINOS)
    MESSAGE(FATAL_ERROR
      "\n-- deal.II was built without support for Trilinos!\n"
      )
  ENDIF()

  IF(NOT DEAL_II_WITH_P4EST)
    MESSAGE(FATAL_ERROR
      "\n-- deal.II was built without support for p4est!\n"
      )
  ENDIF()
ENDIF()

DEAL_II_INITIALIZE_CACHED_VARIABLES()
PROJECT(${TARGET})


IF(EXISTS ${CMAKE_SOURCE_DIR}/tests/CMakeLists.txt)
  ENABLE_TESTING()
  ADD_SUBDIRECTORY(tests)
ENDIF()


MESSAGE("====================================================")

# Depending on whether we link statically or allow for shared libs,
# we can or can not load plugins via external shared libs. Pass this
# down during compilation so we can disable it in the code
IF (DEAL_II_STATIC_EXECUTABLE STREQUAL "ON")
  MESSAGE(STATUS "Creating a statically linked executable")
  MESSAGE(STATUS "Disabling dynamic loading of plugins from the input file")
  FOREACH(_source_file ${TARGET_SRC})
    SET_PROPERTY(SOURCE ${_source_file}
      APPEND PROPERTY COMPILE_DEFINITIONS ASPECT_USE_SHARED_LIBS=0)
  ENDFOREACH()
  SET_PROPERTY(TARGET aspect
     APPEND PROPERTY LINK_FLAGS "-static")
ELSE()
  MESSAGE(STATUS "Enabling dynamic loading of plugins from the input file")
  FOREACH(_source_file ${TARGET_SRC})
    SET_PROPERTY(SOURCE ${_source_file}
      APPEND PROPERTY COMPILE_DEFINITIONS ASPECT_USE_SHARED_LIBS=1)
  ENDFOREACH()
ENDIF()

DEAL_II_INVOKE_AUTOPILOT()

