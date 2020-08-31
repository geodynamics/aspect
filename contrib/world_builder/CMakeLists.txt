cmake_minimum_required (VERSION 2.8.4)
project(WorldBuilder C CXX)

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
SET(WORLD_BUILDER_SOURCE_DIR ${CMAKE_SOURCE_DIR})
INCLUDE("${CMAKE_SOURCE_DIR}/cmake/version.cmake")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

message(STATUS "Building the Geodynamic World Builder version ${WORLD_BUILDER_VERSION} in ${CMAKE_BUILD_TYPE} mode")

set (WORLD_BUILDER_SOURCE_DIR ${PROJECT_SOURCE_DIR})

# generate version.cc
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/source/config.cc.in" "${CMAKE_BINARY_DIR}/source/config.cc" @ONLY)

# finding support for different languages
# finding a fortran compiler
if(${CMAKE_VERSION} VERSION_LESS "2.8.8") 
  message(STATUS "Fortran not supported with cmake version lower than  2.8.8. Please upgrade your cmake version if you need it. Disabling Fortran wrapper and tests.")
  SET(CMAKE_Fortran_COMPILER FALSE)
else()
  include(CheckLanguage)
  check_language(Fortran)
  if(CMAKE_Fortran_COMPILER)
    enable_language(Fortran)
    message(STATUS "Found Fortran support. Enabling Fortran wrapper and tests.")
  else()
    message(STATUS "Did not find Fortran support. Disabling Fortran wrapper and tests.")
  endif()
endif()

#finding support for makeing a python wrapper
SET(MAKE_PYTHON_WRAPPER FALSE)
#find swig
FIND_PACKAGE(SWIG)
if(SWIG_FOUND)
  INCLUDE(${SWIG_USE_FILE})
SET(CMAKE_SWIG_FLAGS "")
else()
  message(STATUS "SWIG was not found, disableling python wrapper compilation.")
endif()

#find python libs
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
      message(STATUS "Python was not found. Disableling python wrapper compilation.")
    endif()

else()
  FIND_PACKAGE(Python COMPONENTS Interpreter Development)
  if(Python_FOUND)
    include_directories(${Python_INCLUDE_DIRS})
  else()
    message(STATUS "Python was not found, disableling python wrapper compilation.")
  endif()
endif()


if(SWIG_FOUND AND Python_FOUND)
  SET(MAKE_PYTHON_WRAPPER TRUE)
  message(STATUS "Both SWIG and Python are found, enableling python wrapper compilation.")
endif()

# Add dependencies are found, now set compiler output directories.
# This needs to be done before the add library and add executable functions.

IF(NOT CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE Debug
       CACHE STRING "Choose the type of build : None Debug Release RelWithDebInfo MinSizeRel Coverage."
       FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

IF(NOT CMAKE_RUN_APP_TESTS)
   SET(CMAKE_RUN_APP_TESTS TRUE
       CACHE STRING "Whether to run or not to run the app tests."
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

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
foreach( OUTPUTCONFIG ${CMAKE_CONFIGURATION_TYPES} )
  string( TOUPPER ${OUTPUTCONFIG} OUTPUTCONFIG )
  set( CMAKE_RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_BINARY_DIR}/bin)
  set( CMAKE_LIBRARY_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_BINARY_DIR}/lib)
  set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_BINARY_DIR}/lib)
endforeach( OUTPUTCONFIG CMAKE_CONFIGURATION_TYPES )

if(MAKE_PYTHON_WRAPPER)
  set(CMAKE_SWIG_OUTDIR ${CMAKE_BINARY_DIR}/lib)
endif()
if(CMAKE_Fortran_COMPILER)
  set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)
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
include_directories("include/" "tests/") 

# Add source directory
file(GLOB_RECURSE SOURCES_CXX "source/*.cc" "${CMAKE_BINARY_DIR}/source/*.cc")

if(CMAKE_Fortran_COMPILER)
  file(GLOB_RECURSE SOURCES_FORTAN "source/*.f90")
endif()

set(SOURCES ${SOURCES_CXX} ${SOURCES_FORTAN})

add_library(WorldBuilder SHARED ${SOURCES})

add_executable(WorldBuilderApp "${CMAKE_CURRENT_SOURCE_DIR}/app/main.cc")
add_executable(WorldBuilderVisualization "${CMAKE_CURRENT_SOURCE_DIR}/visualization/main.cc")

if(MAKE_PYTHON_WRAPPER)

  if(POLICY CMP0078)
    SET(TARGET gwb)
  else()
    SET(TARGET _gwb)
  endif()

  # Add swig module
  SET_SOURCE_FILES_PROPERTIES(SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/source/wrapper_cpp.i PROPERTIES CPLUSPLUS ON)
  if(${CMAKE_VERSION} VERSION_LESS "3.8.0") 
    swig_add_module(gwb python ${CMAKE_CURRENT_SOURCE_DIR}/source/wrapper_cpp.i)
  else()
    swig_add_library(gwb LANGUAGE python SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/source/wrapper_cpp.i)
  endif()

endif()


# Provide "indent" target for indenting all headers and source files
ADD_CUSTOM_TARGET(indent
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  COMMAND ./doc/indent
  COMMENT "Indenting all the World Builder header and source files..."
  )



IF ( CMAKE_BUILD_TYPE STREQUAL Coverage )
  SET(CMAKE_VERBOSE_MAKEFILE on )
  SET(WB_COMPILER_OPTIONS_PRIVATE_COVERAGE_NEW --coverage -fprofile-arcs -ftest-coverage)
  SET(WB_COMPILER_OPTIONS_PRIVATE_COVERAGE_OLD "--coverage -fprofile-arcs -ftest-coverage")

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
   if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER MATCHES "[Cc]lang")
      # Preventing issues with older cmake compilers (<2.8.11) which do not support VERSION_GREATER_EQUAL
      # cmake vesrion 3.12.0 introduces COMPILE_LANGUAGE:FORTAN, otherwise this would be >=2.8.11
      if(NOT ${CMAKE_VERSION} VERSION_LESS "3.12.0") 
        SET(WB_COMPILER_OPTIONS_PRIVATE $<$<COMPILE_LANGUAGE:FORTAN>:-pedantic -fPIC -Wall -Wextra>  
        $<$<COMPILE_LANGUAGE:CXX>:-std=c++11 -pedantic -fPIC -Wall -Wextra -Wmost -Wconversion -Wunreachable-code -Wuninitialized -Wold-style-cast -Wshadow -Wfloat-equal -Wpointer-arith -Wwrite-strings 
        -Wsynth -Wsign-compare -Woverloaded-virtual -Wliteral-range -Wparentheses -Wunused-local-typedefs -Wcast-qual -fstrict-aliasing -Werror=uninitialized -Wundef 
        -Wcast-align -Wmissing-declarations -Wredundant-decls -Wdiv-by-zero -Wdisabled-optimization -Wswitch-default -Wunused>)
      else()
        SET(WB_COMPILER_OPTIONS_PRIVATE "-std=c++11 -pedantic -fPIC -Wall -Wextra -Wmost -Wconversion -Wunreachable-code -Wuninitialized -Wold-style-cast -Wshadow -Wfloat-equal -Wpointer-arith -Wwrite-strings -Wsynth -Wsign-compare -Woverloaded-virtual -Wliteral-range -Wparentheses -Wunused-local-typedefs -Wcast-qual -fstrict-aliasing -Werror=uninitialized -Wundef -Wcast-align -Wmissing-declarations -Wredundant-decls -Wdiv-by-zero -Wdisabled-optimization -Wswitch-default -Wunused")
      endif()
   elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
   SET(WB_COMPILER_OPTIONS_PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-std=c++11 -pthread -pedantic -fPIC -Wall -Wextra -Wpointer-arith -Wwrite-strings -Wsign-compare -Woverloaded-virtual -Wno-parentheses -Wcast-qual -fstrict-aliasing -Wmaybe-uninitialized -Werror=maybe-uninitialized>)
   else()
      # gcc linux

      # Preventing issues with older cmake compilers (<2.8.11) which do not support VERSION_GREATER_EQUAL
      # cmake vesrion 3.12.0 introduces COMPILE_LANGUAGE:FORTAN, otherwise this would be >=2.8.11
      if(NOT ${CMAKE_VERSION} VERSION_LESS "3.12.0") 
        SET(WB_COMPILER_OPTIONS_PRIVATE $<$<COMPILE_LANGUAGE:FORTAN>:-pedantic -fPIC -Wall -Wextra>  $<$<COMPILE_LANGUAGE:CXX>:-std=c++11 -pedantic -fPIC -Wall -Wextra -Wpointer-arith -Wwrite-strings -Wsynth -Wsign-compare -Woverloaded-virtual -Wno-placement-new -Wno-literal-suffix -Wno-parentheses -Wno-unused-local-typedefs -Wcast-qual -fstrict-aliasing -Wmaybe-uninitialized -Werror=maybe-uninitialized -Wparentheses -Wfloat-equal -Wundef -Wcast-align -Wlogical-op -Wmissing-declarations -Wredundant-decls -Wdiv-by-zero -Wdisabled-optimization -Wswitch-default -Wno-unused>)
      else()
        SET(WB_COMPILER_OPTIONS_PRIVATE "-std=c++11 -pedantic -fPIC -Wall -Wextra -Wpointer-arith -Wwrite-strings -Wsynth -Wsign-compare -Woverloaded-virtual -Wno-placement-new -Wno-literal-suffix -Wno-parentheses -Wno-unused-local-typedefs -Wcast-qual -fstrict-aliasing -Wmaybe-uninitialized -Werror=maybe-uninitialized -Wparentheses -Wfloat-equal -Wundef -Wcast-align -Wlogical-op -Wmissing-declarations -Wredundant-decls -Wdiv-by-zero -Wdisabled-optimization -Wswitch-default -Wno-unused>")
      endif()

      # adding flags is a hassle before cmake 3.12.0, and these flags are mostly useful for development, so just ignore those flags.
      if(NOT ${CMAKE_VERSION} VERSION_LESS "3.12.0") # Preventing issues with older cmake compilers which do not support VERSION_GREATER_EQUAL
        include(CheckCXXCompilerFlag)

        check_cxx_compiler_flag(-Wsuggest-override SUGGEST_OVERWRITE_CXX_COMPILE_FLAG)
        if(SUGGEST_OVERWRITE_CXX_COMPILE_FLAG)
          SET(WB_COMPILER_OPTIONS_PRIVATE ${WB_COMPILER_OPTIONS_PRIVATE} $<$<COMPILE_LANGUAGE:CXX>:-Wsuggest-override>)
        endif()

        check_cxx_compiler_flag(-Wno-placement-new NO_PLACEMENT_NEW_CXX_COMPILE_FLAG)
        if(NO_PLACEMENT_NEW_CXX_COMPILE_FLAG)
          SET(WB_COMPILER_OPTIONS_PRIVATE ${WB_COMPILER_OPTIONS_PRIVATE} $<$<COMPILE_LANGUAGE:CXX>:-Wno-placement-new>)
        endif()

        # untill 4.9.0 Wshadow is too strict: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=57709,
        # but still some issues in 4.9.0, so only add it after 5.0.0.
        if(CMAKE_COMPILER_IS_GNUCC AND NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5.0.0)
          SET(WB_COMPILER_OPTIONS_PRIVATE ${WB_COMPILER_OPTIONS_PRIVATE} $<$<COMPILE_LANGUAGE:CXX>:-Wshadow>)
        endif()
      endif()
   endif()
   IF ( CMAKE_BUILD_TYPE STREQUAL Coverage )
     if(NOT ${CMAKE_VERSION} VERSION_LESS "3.13.0") # Preventing issues with older cmake compilers which do not support VERSION_GREATER_EQUAL
     SET(WB_LINKER_OPTIONS -lstdc++ --coverage -fprofile-arcs -ftest-coverage)
     else()
      SET(WB_LINKER_OPTIONS "-lstdc++ --coverage -fprofile-arcs -ftest-coverage")
     endif()
   else()
     SET(WB_LINKER_OPTIONS -lstdc++)
   endif()

   SET(WB_VISU_LINKER_OPTIONS "-pthread")
elseif(APPLE)
    SET(WB_COMPILER_OPTIONS_PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-std=c++11 -Wpointer-arith -Wwrite-strings -Wsynth -Wsign-compare -Woverloaded-virtual -Wno-literal-range -Wno-parentheses -Wno-unused-local-typedefs -Wcast-qual -fstrict-aliasing -stdlib=libc++ -Wuninitialized -Werror=uninitialized -Wdangling-else>)
    SET(WB_LINKER_OPTIONS -lc++)
else()
    #MSVS
    SET(WB_COMPILER_OPTIONS_PRIVATE $<$<COMPILE_LANGUAGE:CXX>:/W3 /EHsc -std=c++11>)
    SET(WB_LINKER_OPTIONS  /WHOLEARCHIVE:${CMAKE_CURRENT_BINARY_DIR}/lib/WorldBuilder.lib)
endif()



if(${CMAKE_VERSION} VERSION_LESS "3.12.0") 
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${WB_COMPILER_OPTIONS_INTERFACE} ${WB_COMPILER_OPTIONS_PRIVATE} ${WB_COMPILER_OPTIONS_PRIVATE_COVERAGE_OLD}")
else()
  target_compile_options(WorldBuilder INTERFACE ${WB_COMPILER_OPTIONS_INTERFACE} PRIVATE ${WB_COMPILER_OPTIONS_PRIVATE} ${WB_COMPILER_OPTIONS_PRIVATE_COVERAGE_NEW})
  target_compile_options(WorldBuilderApp INTERFACE ${WB_COMPILER_OPTIONS_INTERFACE} PRIVATE ${WB_COMPILER_OPTIONS_PRIVATE} ${WB_COMPILER_OPTIONS_PRIVATE_COVERAGE_NEW})
  target_compile_options(WorldBuilderVisualization INTERFACE ${WB_COMPILER_OPTIONS_INTERFACE} PRIVATE ${WB_COMPILER_OPTIONS_PRIVATE} ${WB_COMPILER_OPTIONS_PRIVATE_COVERAGE_NEW})
endif()

if(${CMAKE_VERSION} VERSION_LESS "3.13.0") 
    SET(CMAKE_EXE_LINKER_FLAGS "${WB_LINKER_OPTIONS} ${WB_VISU_LINKER_OPTIONS} ${WB_COMPILER_OPTIONS_PRIVATE_COVERAGE_OLD}")
else()
    target_link_options(WorldBuilder INTERFACE ${WB_LINKER_OPTIONS} ${WB_COMPILER_OPTIONS_PRIVATE_COVERAGE_NEW})
    target_link_options(WorldBuilderApp INTERFACE ${WB_LINKER_OPTIONS} ${WB_COMPILER_OPTIONS_PRIVATE_COVERAGE_NEW})
    target_link_options(WorldBuilderVisualization INTERFACE ${WB_LINKER_OPTIONS} ${WB_COMPILER_OPTIONS_PRIVATE_COVERAGE_NEW})
    target_link_options(WorldBuilderVisualization PRIVATE  ${WB_VISU_LINKER_OPTIONS} ${WB_COMPILER_OPTIONS_PRIVATE_COVERAGE_NEW})
endif()




# Make sure that the whole library is loaded, so the registration is done correctly.
if(NOT APPLE AND NOT MSVC)
  SET(GWB_LIBRARY_WHOLE -Wl,--whole-archive WorldBuilder -Wl,--no-whole-archive)
elseif(MSVC)
  SET(GWB_LIBRARY_WHOLE WorldBuilder)
else()
  SET(GWB_LIBRARY_WHOLE -Wl,-force_load WorldBuilder)
endif()
target_link_libraries (WorldBuilderApp ${GWB_LIBRARY_WHOLE})
target_link_libraries (WorldBuilderVisualization ${GWB_LIBRARY_WHOLE}) 


if(MAKE_PYTHON_WRAPPER)

  # see https://github.com/swig/swig/issues/1259 for more info.
  if(SWIG_VERSION VERSION_LESS "4.0.1" AND NOT MSVC)
    set_target_properties(${TARGET} PROPERTIES COMPILE_FLAGS "-Wno-cast-qual")
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
          COMMAND ${Python_EXECUTABLE} ${SETUP_PY_OUT} install)
endif()

# binary:
install(CODE "message(Installing...)")
install (TARGETS WorldBuilder EXPORT WorldBuilder  DESTINATION bin)
install (EXPORT WorldBuilder  DESTINATION bin)

if(CMAKE_Fortran_COMPILER)
  install (FILES ${CMAKE_CURRENT_BINARY_DIR}/mod/worldbuilder.mod DESTINATION include)
endif()

if(MAKE_PYTHON_WRAPPER)
  install(CODE "execute_process(COMMAND ${Python_EXECUTABLE} ${SETUP_PY_OUT} install)")
endif()

# headers:
INSTALL(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/
  DESTINATION include
  COMPONENT includes
  FILES_MATCHING PATTERN "*.h")

enable_testing()
add_subdirectory(tests)
