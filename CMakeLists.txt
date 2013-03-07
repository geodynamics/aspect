CMAKE_MINIMUM_REQUIRED (VERSION 2.8)

FIND_PACKAGE(deal.II 8.0 REQUIRED
  HINTS
    ${DEAL_II_DIR}
    $ENV{DEAL_II_DIR}
    "../../"
    "~/workspace/install-deal.II"
  )
MESSAGE(STATUS
  "Using the deal.II installation found at ${DEAL_II_PATH}"
  )

IF(NOT DEAL_II_WITH_TRILINOS)
  MESSAGE(FATAL_ERROR
    "\nAspect requires a deal.II installation built with support for Trilinos!\n"
    )
ENDIF()

IF(NOT DEAL_II_WITH_P4EST)
  MESSAGE(FATAL_ERROR
    "\naspect requires deal.II build with support for p4est!\n"
    )
ENDIF()

DEAL_II_INITIALIZE_CACHED_VARIABLES()

SET (VERSION "0.1")
PROJECT (aspect)

INCLUDE_DIRECTORIES(include)

ADD_SUBDIRECTORY(source)
#ADD_SUBDIRECTORY(run)
