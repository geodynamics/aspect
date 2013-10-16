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

DEAL_II_INVOKE_AUTOPILOT()
