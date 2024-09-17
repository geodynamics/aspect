cmake_minimum_required(VERSION 3.13.4)

project(WorldBuilder C CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# SWIG: use standard target name.
if(POLICY CMP0078)
  cmake_policy(SET CMP0078 NEW)
endif()
## SWIG: use SWIG_MODULE_NAME property.
if(POLICY CMP0086)
  cmake_policy(SET CMP0086 NEW)
endif()

if(POLICY CMP0058)
  cmake_policy(SET CMP0058 NEW)
endif()

# load in version info and export it
SET(WORLD_BUILDER_SOURCE_DIR ${PROJECT_SOURCE_DIR})
INCLUDE("${PROJECT_SOURCE_DIR}/cmake/version.cmake")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

message(STATUS "Building the Geodynamic World Builder version ${WORLD_BUILDER_VERSION} in ${CMAKE_BUILD_TYPE} mode")


include(CheckCXXSymbolExists)
check_cxx_symbol_exists(feenableexcept fenv.h HAVE_DECL_FEENABLEEXCEPT)
check_cxx_symbol_exists(feclearexcept fenv.h HAVE_DECL_FECLEAREXCEPT)

IF(HAVE_DECL_FEENABLEEXCEPT AND HAVE_DECL_FECLEAREXCEPT)
  SET(WB_USE_FP_EXCEPTIONS TRUE CACHE BOOL "Whether to turn on floating point exceptions in debug mode.")
  MESSAGE(STATUS "Using Floating Point Exceptions if in Debug mode.")
ELSE()
  MESSAGE(STATUS "Not using Floating Point Exceptions.")
ENDIF()

set (WORLD_BUILDER_SOURCE_DIR ${PROJECT_SOURCE_DIR})

# define some config options:
set(WB_ENABLE_APPS ON CACHE BOOL "If enabled, applications gwb-grid and gwb-dat are built.")
set(WB_ENABLE_HELPER_TARGETS ON CACHE BOOL "If enabled, targets to indenting and switching to debug/release mode are configured.")
set(WB_ENABLE_TESTS ON CACHE BOOL "If enabled, tests are configured.")


# generate version.cc
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/include/world_builder/config.h.in" "${PROJECT_BINARY_DIR}/include/world_builder/config.h" @ONLY)


# determine whether to use clang-tidy
set(WB_CHECK_CLANG_TIDY OFF CACHE BOOL "Whether or not to check the code with clang-tidy while compiling the world builder.")

if(WB_CHECK_CLANG_TIDY)
  set(WB_CHECK_CLANG_TIDY_FIX OFF CACHE BOOL "Whether or not to try to automatically fix the clang-tidy issues when found.")
  if(WB_CHECK_CLANG_TIDY_FIX)
    set(CMAKE_CXX_CLANG_TIDY "clang-tidy;--fix;-header-filter=world_builder/.*")
  else()
    set(CMAKE_CXX_CLANG_TIDY "clang-tidy;-header-filter=world_builder/.*}")
  endif()
endif()

# finding support for different languages
# finding a fortran compiler
set(WB_MAKE_FORTRAN_WRAPPER ON CACHE BOOL "Whether or not to create a Fortran wrapper for the world builder.")
if(WB_MAKE_FORTRAN_WRAPPER)
  include(CheckLanguage)
  message(STATUS "Looking for Fortran support.")
  check_language(Fortran)
  if(CMAKE_Fortran_COMPILER)
    enable_language(Fortran)
  else()
    set(WB_MAKE_FORTRAN_WRAPPER OFF CACHE BOOL "Whether or not to create a Fortran wrapper for the world builder." FORCE)
  endif()
endif()

if(WB_MAKE_FORTRAN_WRAPPER)
  message(STATUS "Found Fortran support. Enabling Fortran wrapper and tests.")
else()
  message(STATUS "Did not find Fortran support or it has been requested to be disabled. Disabling Fortran wrapper and tests.")
endif()

#finding support for making a python wrapper
SET(MAKE_PYTHON_WRAPPER FALSE)
#find swig
FIND_PACKAGE(SWIG)
if(SWIG_FOUND)
  INCLUDE(${SWIG_USE_FILE})
SET(CMAKE_SWIG_FLAGS "")
else()
  message(STATUS "SWIG was not found, disabling python wrapper compilation.")
endif()


# Find MPI if available.
find_package(MPI)

if(MPI_FOUND AND ${CMAKE_VERSION} VERSION_GREATER "3.9.0" )
  message(STATUS "MPI found.")
  set(USE_MPI ON CACHE BOOL "Whether or not to enable MPI when compiling the world builder.")
else()
  message(STATUS "MPI not found or cmake version lower than 3.9.0.")
  set(USE_MPI OFF CACHE BOOL "Whether or not to enable MPI when compiling the world builder." FORCE)
endif()

if(USE_MPI)
  message(STATUS "Using MPI found in: ${MPI_INCLUDE_PATH}")
  include_directories(${MPI_INCLUDE_PATH})
else()
  message(STATUS "Not using MPI.")
endif()


#find python libs
SET(WB_ENABLE_PYTHON ON CACHE BOOL "Whether or not to enable Python when compiling the world builder.") 
if(${WB_ENABLE_PYTHON})
  if(${CMAKE_VERSION} VERSION_LESS "3.12.0") 
    FIND_PACKAGE(PythonLibs 3.0)
    FIND_PACKAGE(PythonInterp 3.0)

      if(PYTHONLIBS_FOUND AND PYTHONINTERP_FOUND)
        SET(Python_FOUND TRUE)
        SET(Python_EXECUTABLE ${PYTHON_EXECUTABLE})
        SET(Python_INCLUDE_DIRS ${PYTHON_INCLUDE_DIRS})
        SET(Python_LIBRARIES ${PYTHON_LIBRARIES})
        SET(Python_LIBRARIES_RELEASE ${PYTHON_LIBRARIES_RELEASE})
        include_directories(${PYTHON_INCLUDE_DIRS})
      else()
        SET(Python_FOUND FALSE)
        message(STATUS "Python was not found. Disabling python wrapper compilation.")
      endif()

  else()
    FIND_PACKAGE(Python COMPONENTS Interpreter Development)
    if(Python_FOUND)
      include_directories(${Python_INCLUDE_DIRS})
    else()
      message(STATUS "Python was not found, disabling python wrapper compilation.")
    endif()
  endif()
else()
  message(STATUS "Python was disabled by cmake options.")
endif()


if(SWIG_FOUND AND Python_FOUND)
  SET(MAKE_PYTHON_WRAPPER TRUE)
  message(STATUS "Both SWIG and Python are found, enableling python wrapper compilation.")
endif()

find_package(ZLIB)

SET(WB_USE_ZLIB ON CACHE BOOL "Whether or not to use zlib for compiling the world builder visualizer.")
if(ZLIB_FOUND AND WB_USE_ZLIB)
  message(STATUS "Enabling the use of ZLIB.")
  include_directories(${ZLIB_INCLUDE_DIRS})
else()
  message(STATUS "Disabling the use of ZLIB.")
  SET(WB_USE_ZLIB OFF CACHE BOOL "Whether or not to use zlib for compiling the world builder visualizer." FORCE)
endif()


# Add dependencies are found, now set compiler output directories.
# This needs to be done before the add library and add executable functions.

IF(NOT CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE Debug
       CACHE STRING "Choose the type of build : None Debug Release RelWithDebInfo MinSizeRel Coverage."
       FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

IF (CMAKE_BUILD_TYPE STREQUAL DebugRelease)
    MESSAGE(FATAL_ERROR "DebugRelease build mode is not supported!")
ENDIF()

IF(NOT WB_RUN_APP_TESTS)
   SET(WB_RUN_APP_TESTS ON
       CACHE BOOL "Whether to run or not to run the app tests."
       )
ENDIF()

IF ( CMAKE_BUILD_TYPE STREQUAL Coverage )
MARK_AS_ADVANCED(
    CMAKE_CXX_FLAGS_COVERAGE
    CMAKE_C_FLAGS_COVERAGE
    CMAKE_Fortran_FLAGS_COVERAGE
    CMAKE_EXE_LINKER_FLAGS_COVERAGE
    CMAKE_SHARED_LINKER_FLAGS_COVERAGE )
ENDIF()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
foreach( OUTPUTCONFIG ${CMAKE_CONFIGURATION_TYPES} )
  string( TOUPPER ${OUTPUTCONFIG} OUTPUTCONFIG )
  set( CMAKE_RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${PROJECT_BINARY_DIR}/bin)
  set( CMAKE_LIBRARY_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${PROJECT_BINARY_DIR}/lib)
  set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${PROJECT_BINARY_DIR}/lib)
endforeach( OUTPUTCONFIG CMAKE_CONFIGURATION_TYPES )

if(MAKE_PYTHON_WRAPPER)
  set(CMAKE_SWIG_OUTDIR ${PROJECT_BINARY_DIR}/lib)
endif()
if(WB_MAKE_FORTRAN_WRAPPER)
  set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/mod)
endif()

if(MAKE_PYTHON_WRAPPER)
  # Files to install with Python
  if(WIN32)
    set(PYTHON_INSTALL_FILES
            ${CMAKE_CURRENT_BINARY_DIR}/lib/gwb.py
            ${CMAKE_CURRENT_BINARY_DIR}/lib/_gwb.pyd)
  else()
    set(PYTHON_INSTALL_FILES
            ${CMAKE_CURRENT_BINARY_DIR}/lib/gwb.py
            ${CMAKE_CURRENT_BINARY_DIR}/lib/_gwb.so)
  endif()
  # Configure setup.py and copy to output directory
  set(SETUP_PY_IN ${CMAKE_CURRENT_SOURCE_DIR}/cmake/setup.py.in)
  set(SETUP_PY_OUT ${CMAKE_CURRENT_BINARY_DIR}/lib/setup.py)
  configure_file(${SETUP_PY_IN} ${SETUP_PY_OUT})
endif()

# Add include directory
include_directories("include/" "tests/" "${CMAKE_CURRENT_BINARY_DIR}/include")

# Add source directory
file(GLOB_RECURSE SOURCES_CXX "source/world_builder/*.cc" "${PROJECT_BINARY_DIR}/source/world_builder/*.cc")

if(WB_MAKE_FORTRAN_WRAPPER)
  file(GLOB_RECURSE SOURCES_FORTRAN "source/world_builder/*.f90")
endif()

set(SOURCES ${SOURCES_CXX} ${SOURCES_FORTRAN})

set(CMAKE_UNITY_BUILD_BATCH_SIZE  5000)

# Special treatment of some files for unity builds. We disable the unity
# build for them.
SET(UNITY_DISABLE_FILES
  "source/world_builder/parameters.cc;source/world_builder/point.cc;")

FOREACH(_source_file ${UNITY_DISABLE_FILES})
  SET(_full_name "${PROJECT_SOURCE_DIR}/${_source_file}")
  LIST(FIND SOURCES ${_full_name} _index)
  IF(_index EQUAL -1)
    MESSAGE(FATAL_ERROR "could not find ${_full_name}.")
  ENDIF()

  SET_PROPERTY(SOURCE ${_source_file} PROPERTY SKIP_UNITY_BUILD_INCLUSION TRUE)
ENDFOREACH()

if(NOT DEFINED WB_TARGET)
  set(WB_TARGET "WorldBuilder" CACHE STRING "")
endif()
add_library(${WB_TARGET} ${SOURCES})
set_property(TARGET ${WB_TARGET} PROPERTY POSITION_INDEPENDENT_CODE ON)

IF(NOT CMAKE_VERSION VERSION_LESS 3.16)
  SET(WB_UNITY_BUILD ON CACHE BOOL "Combine source files into less compile targets to speedup compile time. Currently only supported for cmake 3.16 and newer versions.")
ELSE()
  SET(WB_UNITY_BUILD OFF CACHE BOOL "Combine source files into less compile targets to speedup compile time. Currently only supported for cmake 3.16 and newer versions." FORCE)
ENDIF()

IF (WB_UNITY_BUILD)
  IF (NOT CMAKE_VERSION VERSION_LESS 3.16)
    SET_PROPERTY(TARGET ${WB_TARGET} PROPERTY UNITY_BUILD TRUE)
    MESSAGE(STATUS "Combining source files into unity build.")
  ELSEIF(CMAKE_VERSION VERSION_LESS 3.16)
    MESSAGE(FATAL_ERROR "ASPECT_UNITY_BUILD is currently only supported for CMake 3.16 and newer versions.")
  ENDIF()
ELSE()
  IF (NOT CMAKE_VERSION VERSION_LESS 3.16)
    SET_PROPERTY(TARGET ${WB_TARGET} PROPERTY UNITY_BUILD FALSE)
  ENDIF()

  MESSAGE(STATUS "Disabling unity build.")
ENDIF()

if(WB_ENABLE_APPS)
  add_executable(gwb-dat "${CMAKE_CURRENT_SOURCE_DIR}/source/gwb-dat/main.cc")
  add_executable(gwb-grid "${CMAKE_CURRENT_SOURCE_DIR}/source/gwb-grid/main.cc")
endif()

if(MAKE_PYTHON_WRAPPER)

  if(POLICY CMP0078)
    SET(TARGET gwb)
  else()
    SET(TARGET _gwb)
  endif()

  # Add swig module
  SET_SOURCE_FILES_PROPERTIES(SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/source/world_builder/wrapper_cpp.i PROPERTIES CPLUSPLUS ON)
  if(${CMAKE_VERSION} VERSION_LESS "3.8.0") 
    swig_add_module(gwb python ${CMAKE_CURRENT_SOURCE_DIR}/source/world_builder/wrapper_cpp.i)
  else()
    swig_add_library(gwb LANGUAGE python SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/source/world_builder/wrapper_cpp.i)
  endif()

endif()

# Find if we are using make or ninja. 
IF(CMAKE_GENERATOR MATCHES "Ninja")
  SET(_make_command "ninja")
ELSE()
  SET(_make_command "make")
ENDIF()

# Provide "indent" target for indenting all headers and source files
IF(CMAKE_SOURCE_DIR MATCHES "${PROJECT_SOURCE_DIR}")
  SET(TARGET_WB_PREFIX "")
ELSE()
  SET(TARGET_WB_PREFIX "wb_")
ENDIF()

if(WB_ENABLE_HELPER_TARGETS)
  # Provide "indent" target for indenting all headers and source files
  ADD_CUSTOM_TARGET(${TARGET_WB_PREFIX}indent
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    COMMAND ./doc/indent
    COMMENT "Indenting all the World Builder header and source files..."
    )

  # Provide "release" target for switching to release mode
  ADD_CUSTOM_TARGET(${TARGET_WB_PREFIX}release
    COMMAND ${CMAKE_COMMAND} -D CMAKE_BUILD_TYPE=Release .
    COMMAND ${CMAKE_COMMAND} -E echo "Switched to Release mode. Now recompile with: ${_make_command}"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    VERBATIM
    COMMENT "switching to RELEASE mode..."
    )


  # Provide "debug" target for switching to release mode
  ADD_CUSTOM_TARGET(${TARGET_WB_PREFIX}debug
    COMMAND ${CMAKE_COMMAND} -D CMAKE_BUILD_TYPE=Debug .
    COMMAND ${CMAKE_COMMAND} -E echo "Switched to Debug mode. Now recompile with: ${_make_command}"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    VERBATIM
    COMMENT "switching to DEBUG mode..."
    )
endif()

SET(FORCE_COLORED_OUTPUT ON CACHE BOOL "Forces colored output when compiling with gcc and clang.")


IF ( CMAKE_BUILD_TYPE STREQUAL Coverage )
  SET(CMAKE_VERBOSE_MAKEFILE on )
  SET(WB_COMPILER_OPTIONS_PRIVATE_COVERAGE_NEW --coverage -fprofile-arcs -ftest-coverage -g)
  SET(WB_COMPILER_OPTIONS_PRIVATE_COVERAGE_OLD "--coverage -fprofile-arcs -ftest-coverage -g")

  ## for old versions of cmake
  SET(WB_FORTRAN_COMPILER_FLAGS_COVERAGE "--coverage")
endif()

  SET(WB_LINKER_OPTIONS "")
  SET(WB_VISU_LINKER_OPTIONS "")
  SET(WB_COMPILER_OPTIONS_INTERFACE "")
  SET(WB_COMPILER_OPTIONS_PRIVATE "")
  SET(WB_VISU_COMPILER_OPTIONS_PRIVATE "")

if (NOT MSVC AND NOT APPLE)
   #TODO: if other compiles need to be threaded, see https://computing.llnl.gov/tutorials/pthreads/#Compiling to add the correct cases.
   if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER MATCHES "AppleClang")
      # Preventing issues with older cmake compilers (<3.7) which do not support VERSION_GREATER_EQUAL
      # cmake vesrion 3.12.0 introduces COMPILE_LANGUAGE:FORTRAN, otherwise this would be >=2.8.12
      if(NOT ${CMAKE_VERSION} VERSION_LESS "3.9.0")
        SET(WB_COMPILER_OPTIONS_PRIVATE -pedantic -Wall -Wextra  
        $<$<COMPILE_LANGUAGE:CXX>:-Wmost  -Wunreachable-code -Wuninitialized -Wmissing-braces -Wunused-parameter -Wold-style-cast -Wshadow -Wfloat-equal -Wpointer-arith -Wwrite-strings
        -Wsynth -Wsign-compare -Woverloaded-virtual -Wliteral-range -Wparentheses -Wunused-local-typedefs -Wcast-qual -fstrict-aliasing -Werror=uninitialized -Wundef 
        -Wcast-align -Wmissing-declarations -Wredundant-decls -Wdiv-by-zero -Wdisabled-optimization -Wswitch-default -Wunused>)

        if (${FORCE_COLORED_OUTPUT})
          SET(WB_COMPILER_OPTIONS_PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-fcolor-diagnostics> ${WB_COMPILER_OPTIONS_PRIVATE})
        endif ()

      else()
        SET(WB_COMPILER_OPTIONS_PRIVATE "-pedantic -Wall -Wextra -Wmost -Wunreachable-code -Wuninitialized -Wmissing-braces -Wunused-parameter -Wold-style-cast -Wshadow -Wfloat-equal -Wpointer-arith -Wwrite-strings -Wsynth -Wsign-compare -Woverloaded-virtual -Wliteral-range -Wparentheses -Wunused-local-typedefs -Wcast-qual -fstrict-aliasing -Werror=uninitialized -Wundef -Wcast-align -Wmissing-declarations -Wredundant-decls -Wdiv-by-zero -Wdisabled-optimization -Wswitch-default -Wunused")

        if (${FORCE_COLORED_OUTPUT})
          SET(WB_COMPILER_OPTIONS_PRIVATE "-fcolor-diagnostics ${WB_COMPILER_OPTIONS_PRIVATE}")
        endif ()

      endif()

   elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
     SET(WB_COMPILER_OPTIONS_PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-pthread -pedantic -Wall -Wextra -Wpointer-arith -Wwrite-strings -Wsign-compare 
     -Woverloaded-virtual -Wno-parentheses -Wcast-qual -fstrict-aliasing -Wmaybe-uninitialized -Werror=maybe-uninitialized>)
   else()
      # gcc linux

      # Preventing issues with older cmake compilers (<3.7) which do not support VERSION_GREATER_EQUAL
      # cmake vesrion 3.12.0 introduces COMPILE_LANGUAGE:FORTRAN, otherwise this would be >=2.8.12
      # -Warith-conversion was added in gcc 10, add back in when the minimum becomes gcc 10.
      if(NOT ${CMAKE_VERSION} VERSION_LESS "3.9.0")
        SET(WB_COMPILER_OPTIONS_PRIVATE -pedantic -Wall -Wextra  $<$<COMPILE_LANGUAGE:CXX>:-Wunused-variable -Wold-style-cast -Wconversion -Wsign-conversion -Wmissing-braces -Wunused-parameter -Wpointer-arith -Wwrite-strings -Wsynth -Wsign-compare -Woverloaded-virtual -Wno-placement-new -Wno-literal-suffix -Wno-parentheses -Wno-unused-local-typedefs -Wcast-qual -fstrict-aliasing -Wmaybe-uninitialized -Werror=maybe-uninitialized -Wparentheses -Wfloat-equal -Wundef -Wcast-align -Wlogical-op -Wmissing-declarations -Wredundant-decls -Wdiv-by-zero -Wdisabled-optimization -Wswitch-default -Wno-unused>)

        if (${FORCE_COLORED_OUTPUT})
          SET(WB_COMPILER_OPTIONS_PRIVATE -fdiagnostics-color=always ${WB_COMPILER_OPTIONS_PRIVATE})
        endif()

      else()
        SET(WB_COMPILER_OPTIONS_PRIVATE "-pedantic -Wall -Wextra -Wold-style-cast -Wconversion -Wsign-conversion -Wunused-variable -Wmissing-braces -Wunused-parameter -Wpointer-arith -Wwrite-strings -Wsynth -Wsign-compare -Woverloaded-virtual -Wno-placement-new -Wno-literal-suffix -Wno-parentheses -Wno-unused-local-typedefs -Wcast-qual -fstrict-aliasing -Wmaybe-uninitialized -Werror=maybe-uninitialized -Wparentheses -Wfloat-equal -Wundef -Wcast-align -Wlogical-op -Wmissing-declarations -Wredundant-decls -Wdiv-by-zero -Wdisabled-optimization -Wswitch-default -Wno-unused>")

        if (${FORCE_COLORED_OUTPUT})
          SET(WB_COMPILER_OPTIONS_PRIVATE "-fdiagnostics-color=always ${WB_COMPILER_OPTIONS_PRIVATE}")
        endif()

      endif()

      # adding flags is a hassle before cmake 3.12.0, and these flags are mostly useful for development, so just ignore those flags.
      if(NOT ${CMAKE_VERSION} VERSION_LESS "3.9.0") # Preventing issues with older cmake compilers which do not support VERSION_GREATER_EQUAL
        include(CheckCXXCompilerFlag)

        check_cxx_compiler_flag(-Wsuggest-override SUGGEST_OVERWRITE_CXX_COMPILE_FLAG)
        if(SUGGEST_OVERWRITE_CXX_COMPILE_FLAG)
          SET(WB_COMPILER_OPTIONS_PRIVATE ${WB_COMPILER_OPTIONS_PRIVATE} $<$<COMPILE_LANGUAGE:CXX>:-Wsuggest-override>)
        endif()

        check_cxx_compiler_flag(-Wno-placement-new NO_PLACEMENT_NEW_CXX_COMPILE_FLAG)
        if(NO_PLACEMENT_NEW_CXX_COMPILE_FLAG)
          SET(WB_COMPILER_OPTIONS_PRIVATE ${WB_COMPILER_OPTIONS_PRIVATE} $<$<COMPILE_LANGUAGE:CXX>:-Wno-placement-new>)
        endif()

        # until 4.9.0 Wshadow is too strict: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=57709,
        # but still some issues in 4.9.0, so only add it after 5.0.0.
        if(CMAKE_COMPILER_IS_GNUCC AND NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5.0.0)
          SET(WB_COMPILER_OPTIONS_PRIVATE ${WB_COMPILER_OPTIONS_PRIVATE} $<$<COMPILE_LANGUAGE:CXX>:-Wshadow>)
        endif()
      endif()
   endif()
   IF ( CMAKE_BUILD_TYPE STREQUAL Coverage )
     if(NOT ${CMAKE_VERSION} VERSION_LESS "3.13.0") # Preventing issues with older cmake compilers which do not support VERSION_GREATER_EQUAL
     SET(WB_LINKER_OPTIONS --coverage -fprofile-arcs -ftest-coverage)
     else()
      SET(WB_LINKER_OPTIONS "--coverage -fprofile-arcs -ftest-coverage")
     endif()
   endif()

   SET(WB_VISU_LINKER_OPTIONS "-pthread")
elseif(APPLE)
    SET(WB_COMPILER_OPTIONS_PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-Wpointer-arith -Wwrite-strings -Wsynth -Wsign-compare -Woverloaded-virtual -Wno-literal-range -Wno-parentheses -Wno-unused-local-typedefs -Wcast-qual -fstrict-aliasing -Wuninitialized -Werror=uninitialized -Wdangling-else>)
else()
    #MSVS
    SET(WB_COMPILER_OPTIONS_PRIVATE $<$<COMPILE_LANGUAGE:CXX>:/W3 /EHsc>)
    SET(WB_LINKER_OPTIONS  /WHOLEARCHIVE:${CMAKE_CURRENT_BINARY_DIR}/lib/WorldBuilder.lib)
endif()



if(${CMAKE_VERSION} VERSION_LESS "3.9.0")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${WB_COMPILER_OPTIONS_INTERFACE} ${WB_COMPILER_OPTIONS_PRIVATE} ${WB_COMPILER_OPTIONS_PRIVATE_COVERAGE_OLD}")
else()
  target_compile_options(${WB_TARGET} INTERFACE ${WB_COMPILER_OPTIONS_INTERFACE} PRIVATE ${WB_COMPILER_OPTIONS_PRIVATE} ${WB_COMPILER_OPTIONS_PRIVATE_COVERAGE_NEW})
  if (WB_ENABLE_APPS)
    target_compile_options(gwb-dat INTERFACE ${WB_COMPILER_OPTIONS_INTERFACE} PRIVATE ${WB_COMPILER_OPTIONS_PRIVATE} ${WB_COMPILER_OPTIONS_PRIVATE_COVERAGE_NEW})
    target_compile_options(gwb-grid INTERFACE ${WB_COMPILER_OPTIONS_INTERFACE} PRIVATE ${WB_COMPILER_OPTIONS_PRIVATE} ${WB_COMPILER_OPTIONS_PRIVATE_COVERAGE_NEW})
  endif()
endif()

if(${CMAKE_VERSION} VERSION_LESS "3.13.0") 
    SET(CMAKE_EXE_LINKER_FLAGS "${WB_LINKER_OPTIONS} ${WB_VISU_LINKER_OPTIONS} ${WB_COMPILER_OPTIONS_PRIVATE_COVERAGE_OLD}")
else()
    target_link_options(${WB_TARGET} INTERFACE ${WB_LINKER_OPTIONS})
    if (WB_ENABLE_APPS)
      target_link_options(gwb-dat INTERFACE ${WB_LINKER_OPTIONS} )
      target_link_options(gwb-grid INTERFACE ${WB_LINKER_OPTIONS})
      target_link_options(gwb-grid PRIVATE  ${WB_VISU_LINKER_OPTIONS})
    endif()
endif()

if(${WB_USE_FP_EXCEPTIONS})
  add_definitions(-DWB_USE_FP_EXCEPTIONS)
endif()

if(${USE_MPI})
  add_definitions(-DWB_WITH_MPI)
endif()

if(${WB_USE_ZLIB})
  add_definitions(-DWB_WITH_ZLIB -DVTU11_ENABLE_ZLIB)
endif()


# Make sure that the whole library is loaded, so the registration is done correctly.
if(NOT APPLE AND NOT MSVC)
  SET(GWB_LIBRARY_WHOLE -Wl,--whole-archive ${WB_TARGET} -Wl,--no-whole-archive)
elseif(MSVC)
  SET(GWB_LIBRARY_WHOLE ${WB_TARGET})
else()
  SET(GWB_LIBRARY_WHOLE -Wl,-force_load ${WB_TARGET})
endif()

if(${USE_MPI})
  target_link_libraries (${WB_TARGET} PUBLIC MPI::MPI_CXX MPI::MPI_C)
endif()

if(WB_ENABLE_APPS)
  set(_apps "gwb-dat" "gwb-grid")
  set(_libs ${GWB_LIBRARY_WHOLE})
  if(${WB_USE_ZLIB})
    set(_libs ${_libs} ${ZLIB_LIBRARIES})
  endif()
  if (USE_MPI)
    set(_mpi_dep "PUBLIC MPI::MPI_CXX MPI::MPI_C")
  else()
    set(_mpi_dep "")
  endif()

  foreach (_app ${_apps})
    target_link_libraries (${_app} ${_mpi_deb} ${_libs})
  endforeach()
endif()

if(MAKE_PYTHON_WRAPPER)

  # see https://github.com/swig/swig/issues/1259 for more info.
  if(SWIG_VERSION VERSION_LESS "4.0.1" AND NOT MSVC)
    set_target_properties(${TARGET} PROPERTIES COMPILE_FLAGS "-Wno-cast-qual")
  elseif(SWIG_VERSION VERSION_GREATER_EQUAL "4.0.1" AND APPLE)
    # see https://github.com/swig/swig/issues/1790 for more info.
    # this could be made more specific, but decrecation for the rest will be
    # testers will still be tested. TODO: only do this for the file which needs it
    set_target_properties(${TARGET} PROPERTIES COMPILE_FLAGS "-Wno-error=deprecated-declarations")
  endif()

  IF ( CMAKE_BUILD_TYPE STREQUAL Coverage )
    set_target_properties(${TARGET} PROPERTIES COMPILE_FLAGS "--coverage -fprofile-arcs -ftest-coverage")
  endif()

  # visual studio has a problem with finding the debug version of python (pythonxx_d.dll),
  # so force release usage.
  if(WIN32 AND MSVC)
    SET(CMAKE_SWIG_FLAGS "-D_SWIG_WIN32 /WHOLEARCHIVE:${CMAKE_CURRENT_BINARY_DIR}/lib/WorldBuilder.lib")
    swig_link_libraries(gwb ${GWB_LIBRARY_WHOLE} ${Python_LIBRARY_RELEASE})

  else()
    IF ( CMAKE_BUILD_TYPE STREQUAL Coverage )
      swig_link_libraries(gwb ${GWB_LIBRARY_WHOLE} ${Python_LIBRARIES} --coverage -fprofile-arcs -ftest-coverage)
    else()
      swig_link_libraries(gwb ${GWB_LIBRARY_WHOLE} ${Python_LIBRARIES})
    endif()
  endif()
  
  set_target_properties(${TARGET} PROPERTIES PREFIX "_" OUTPUT_NAME "gwb")


  # Install target to call setup.py
  add_custom_target(install-python
          DEPENDS _gwb
          BYPRODUCTS _gwb
          COMMAND ${Python_EXECUTABLE} ${SETUP_PY_OUT} install)
endif()

# binary:
install(CODE "message(Installing...)")
install(TARGETS ${WB_TARGET} EXPORT WorldBuilder DESTINATION bin)
install(EXPORT WorldBuilder  DESTINATION bin)

if(WB_MAKE_FORTRAN_WRAPPER)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/mod/worldbuilder.mod DESTINATION include)
endif()

if(MAKE_PYTHON_WRAPPER)
  install(CODE "execute_process(COMMAND ${Python_EXECUTABLE} ${SETUP_PY_OUT} install)")
endif()

# headers:
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/
  DESTINATION include
  COMPONENT includes
  FILES_MATCHING PATTERN "*.h")

if(WB_ENABLE_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

# Add the cmake folder so the FindSphinx module is found
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

add_subdirectory (doc)
