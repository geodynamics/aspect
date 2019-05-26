
cmake_minimum_required (VERSION 2.6)
project(WorldBuilder C CXX Fortran)

# load in version info and export it
SET(WORLD_BUILDER_SOURCE_DIR ${CMAKE_SOURCE_DIR})
INCLUDE("${CMAKE_SOURCE_DIR}/cmake/version.cmake")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

message(STATUS "Building the Geodynamic World Builder version ${WORLD_BUILDER_VERSION} in ${CMAKE_BUILD_TYPE} mode")

set (WORLD_BUILDER_SOURCE_DIR ${PROJECT_SOURCE_DIR})

# generate version.cc
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/source/config.cc.in" "${CMAKE_BINARY_DIR}/source/config.cc" @ONLY)


# Add include directory
include_directories("include/" "tests/") 

# Add source directory
file(GLOB_RECURSE SOURCES "source/*.cc" "source/*.f90" "${CMAKE_BINARY_DIR}/source/*.cc")

# Provide "indent" target for indenting all headers and source files
ADD_CUSTOM_TARGET(indent
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  COMMAND ./doc/indent
  COMMENT "Indenting all the World Builder header and source files..."
  )

if (NOT MSVC AND NOT APPLE)
   if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER MATCHES "[Cc]lang")
       SET(CMAKE_CXX_FLAGS "-std=c++11 -fopenmp -pedantic -fPIC -Wall -Wextra -Wpointer-arith -Wwrite-strings -Wsynth -Wsign-compare -Wswitch -Woverloaded-virtual -Wno-literal-range -Wno-parentheses -Wno-unused-local-typedefs -Wcast-qual -fstrict-aliasing -Wuninitialized -Werror=uninitialized")
   elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
       SET(CMAKE_CXX_FLAGS "-std=c++11 -fopenmp -pedantic -fPIC -Wall -Wextra -Wpointer-arith -Wwrite-strings -Wsign-compare -Woverloaded-virtual -Wno-parentheses -Wcast-qual -fstrict-aliasing -Wmaybe-uninitialized -Werror=maybe-uninitialized")
   else()
       SET(CMAKE_CXX_FLAGS "-std=c++11 -fopenmp -pedantic -fPIC -Wall -Wextra -Wpointer-arith -Wwrite-strings -Wsynth -Wsign-compare -Woverloaded-virtual -Wno-placement-new -Wno-literal-suffix -Wno-parentheses -Wno-unused-local-typedefs -Wcast-qual -fstrict-aliasing -Wmaybe-uninitialized -Werror=maybe-uninitialized")
   endif()
   SET(CMAKE_EXE_LINKER_FLAGS "-lstdc++")
elseif(APPLE)
   SET(CMAKE_CXX_FLAGS "-std=c++11 -Xpreprocessor -fopenmp ${OpenMP_CXX_FLAGS} -pedantic -fPIC -Wall -Wextra -Wpointer-arith -Wwrite-strings -Wsynth -Wsign-compare -Woverloaded-virtual -Wno-literal-range -Wno-parentheses -Wno-unused-local-typedefs -Wcast-qual -fstrict-aliasing -stdlib=libc++ -Wuninitialized -Werror=uninitialized")
   SET(CMAKE_EXE_LINKER_FLAGS "-lc++${OpenMP_EXE_LINKER_FLAGS}") #-stdlib=libc++")
else()
    SET(CMAKE_CXX_FLAGS "/W3")
endif()


IF(NOT CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE Debug
       CACHE STRING "Choose the type of build : None Debug Release RelWithDebInfo MinSizeRel Coverage."
       FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

IF ( CMAKE_BUILD_TYPE STREQUAL Coverage )
  SET( CMAKE_VERBOSE_MAKEFILE on )
  SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} ")
  SET(CMAKE_CXX_FLAGS_COVERAGE "${CMAKE_CXX_FLAGS_DEBUG} --coverage -fprofile-arcs -ftest-coverage ")
  SET(CMAKE_Fortran_FLAGS_COVERAGE --coverage)
  SET(CMAKE_EXE_LINKER_FLAGS_COVERAGE "${CMAKE_EXE_LINKER_FLAGS_DEBUG} --coverage -fprofile-arcs -ftest-coverage ")
  SET(CMAKE_SHARED_LINKER_FLAGS_COVERAGE "${CMAKE_SHARED_LINKER_FLAGS_DEBUG} --coverage -fprofile-arcs -ftest-coverage ")

MARK_AS_ADVANCED(
    CMAKE_CXX_FLAGS_COVERAGE
    CMAKE_C_FLAGS_COVERAGE
    CMAKE_Fortran_FLAGS_COVERAGE
    CMAKE_EXE_LINKER_FLAGS_COVERAGE
    CMAKE_SHARED_LINKER_FLAGS_COVERAGE )
ENDIF()

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

add_library(WorldBuilder ${SOURCES})
add_executable(WorldBuilderApp "${CMAKE_CURRENT_SOURCE_DIR}/app/main.cc")
add_executable(WorldBuilderVisualization "${CMAKE_CURRENT_SOURCE_DIR}/visualization/main.cc")

# Make sure that the whole library is loaded, so the registration is done correctly.
if(NOT APPLE)
  SET(GWB_LIBRARY_WHOLE -Wl,--whole-archive WorldBuilder -Wl,--no-whole-archive)
else()
  SET(GWB_LIBRARY_WHOLE -Wl,-force_load WorldBuilder)
endif()
target_link_libraries (WorldBuilderApp ${GWB_LIBRARY_WHOLE})
if(NOT APPLE)
  target_link_libraries (WorldBuilderVisualization ${GWB_LIBRARY_WHOLE}) 
else()
  target_link_libraries (WorldBuilderVisualization ${GWB_LIBRARY_WHOLE} omp) 
endif()

# binary:
install (TARGETS WorldBuilder EXPORT WorldBuilder  DESTINATION bin)
install (DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION include)
install (EXPORT WorldBuilder  DESTINATION bin)

# headers:
INSTALL(DIRECTORY include/
  DESTINATION include
  COMPONENT includes
  FILES_MATCHING PATTERN "*.h")

enable_testing()
add_subdirectory(tests)
