# Copyright (C) 2011 - 2024 by the authors of the ASPECT code.
#
# This file is part of ASPECT.
#
# ASPECT is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# ASPECT is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ASPECT; see the file LICENSE.  If not see
# <http://www.gnu.org/licenses/>.

################### top matter ################
cmake_minimum_required(VERSION 3.13.4)

if (CMP0058)
 cmake_policy(CMP0058 NEW)
endif()

project(testsuite CXX)

set(ASPECT_BINARY "${ASPECT_BINARY}" CACHE STRING "" FORCE)
set(Aspect_DIR "${Aspect_DIR}" CACHE STRING "" FORCE)
set(ASPECT_COMPARE_TEST_RESULTS "${ASPECT_COMPARE_TEST_RESULTS}" CACHE BOOL "" FORCE)
set(ASPECT_RUN_ALL_TESTS OFF CACHE BOOL "")

find_package(Aspect 3.1.0 QUIET REQUIRED HINTS ${Aspect_DIR})
DEAL_II_INITIALIZE_CACHED_VARIABLES()

find_package(Perl)

message(STATUS "Setting up tests with CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE}")


enable_testing()

macro(SET_IF_EMPTY _variable _value)
  if("${${_variable}}" STREQUAL "")
    set(${_variable} ${_value})
  endif()
endmacro()

# suppress mac os warning
if(APPLE AND NOT DEFINED CMAKE_MACOSX_RPATH)
  set(CMAKE_MACOSX_RPATH 0) # see policy CMP0042
endif()

# A function that extracts from a file (presumably a .prm file)
# the number of MPI processes this test is to be invoked as.
# This is encoded in .prm files through lines of the form
#    '# MPI: 4'
# The result is returned in a variable _mpi_count in the
# caller's scope.
function(get_mpi_count _filename)
  file(STRINGS ${_filename} _input_lines
       REGEX "MPI:")
  if("${_input_lines}" STREQUAL "")
    set(_mpi_count 1 PARENT_SCOPE)
  else()
    # go over the (possibly multiple) lines with MPI markers and choose the last
    foreach(_input_line ${_input_lines})
     set(_last_line ${_input_line})
    endforeach()
    string(REGEX REPLACE "^ *# *MPI: *([0-9]+) *$" "\\1"
           _mpi_count ${_last_line})
    set(_mpi_count "${_mpi_count}" PARENT_SCOPE)
  endif()
endfunction()

# A function that checks if a test depends on another test.
# This is encoded in .prm files through lines of the form
#    '# DEPENDS-ON: testname'
# The result is returned in a variable _depends_on in the
# caller's scope. Having test dependencies is helpful if
# one test requires another test to finish first, for example
# because the result of the first test is used by the later test.
function(get_depends_on _filename)
  file(STRINGS ${_filename} _input_lines
       REGEX "DEPENDS-ON:")
  if("${_input_lines}" STREQUAL "")
    set(_depends_on "" PARENT_SCOPE)
  else()
    # go over the (possibly multiple) lines with DEPENDS-ON markers and choose the last
    foreach(_input_line ${_input_lines})
     set(_last_line ${_input_line})
    endforeach()
    string(REGEX REPLACE "^ *# *DEPENDS-ON: *(.*) *$" "\\1"
           _depends_on ${_last_line})
    set(_depends_on "${_depends_on}" PARENT_SCOPE)
  endif()
endfunction()

# Analyze the .prm to decide if this test should be run or not. We are
# checking for the following strings in the .prm:
# 1. "Enable if: ASPECT_WITH_WORLD_BUILDER"
# 2. "Enable if: ASPECT_WITH_FASTSCAPE"
# 3. "Enable if: ASPECT_WITH_PYTHON"
# 4. "Enable if: QUICK_TEST" - enable the test even if RUN_ALL_TESTS is false
function(SHOULD_ENABLE_TEST _filename)

  file(STRINGS ${_filename} _input_lines
       REGEX "Enable if: ASPECT_WITH_WORLD_BUILDER")
  if(NOT "${_input_lines}" STREQUAL "")
    if (${_input_lines} MATCHES "Enable if: ASPECT_WITH_WORLD_BUILDER"
        AND NOT ASPECT_WITH_WORLD_BUILDER)
      set(_use_test OFF PARENT_SCOPE)
    endif()
  endif()

  file(STRINGS ${_filename} _input_lines
       REGEX "Enable if: ASPECT_WITH_FASTSCAPE")
  if(NOT "${_input_lines}" STREQUAL "")
    if (${_input_lines} MATCHES "Enable if: ASPECT_WITH_FASTSCAPE"
        AND NOT ASPECT_WITH_FASTSCAPE)
      set(_use_test OFF PARENT_SCOPE)
    endif()
  endif()

  file(STRINGS ${_filename} _input_lines
       REGEX "Enable if: ASPECT_WITH_PYTHON")
  if(NOT "${_input_lines}" STREQUAL "")
    if (${_input_lines} MATCHES "Enable if: ASPECT_WITH_PYTHON"
        AND NOT ASPECT_WITH_PYTHON)
      set(_use_test OFF PARENT_SCOPE)
    endif()
  endif()

  file(STRINGS ${_filename} _input_lines
       REGEX "Enable if: QUICK_TEST")
  if(NOT ASPECT_RUN_ALL_TESTS AND "${_input_lines}" STREQUAL "")
    set(_use_test OFF PARENT_SCOPE)
  endif()
endfunction()


# set a time limit of 10 minutes per test. this should be long
# enough even for long-running tests, and short enough to not
# get into trouble should we have an infinite loop.
SET_IF_EMPTY(TEST_TIME_LIMIT 600)

############################3

add_custom_target(tests)

#
# We need a diff tool, preferably numdiff otherwise diff. Let the user
# override it by specifying TEST_DIFF.
#
find_program(DIFF_EXECUTABLE
  NAMES diff
  HINTS ${DIFF_DIR}
  PATH_SUFFIXES bin
  )

find_program(NUMDIFF_EXECUTABLE
  NAMES numdiff
  HINTS ${NUMDIFF_DIR}
  PATH_SUFFIXES bin
  )

mark_as_advanced(DIFF_EXECUTABLE NUMDIFF_EXECUTABLE)

if("${TEST_DIFF}" STREQUAL "")
  #

  if(NOT NUMDIFF_EXECUTABLE MATCHES "-NOTFOUND")
    set(TEST_DIFF ${NUMDIFF_EXECUTABLE})
  elseif(NOT DIFF_EXECUTABLE MATCHES "-NOTFOUND")
    set(TEST_DIFF ${DIFF_EXECUTABLE})
  else()
    message(FATAL_ERROR
      "Could not find diff or numdiff. One of those are required for running the testsuite.\n"
      "Please specify TEST_DIFF by hand."
      )
  endif()
endif()

file(GLOB _tests *.prm)

set(_n_tests "0")
list(SORT _tests)

foreach(_test ${_tests})
  set(_test_full ${_test})
  get_filename_component(_testname ${_test} NAME_WE)

  set(_use_test ON)

  if ("${_testname}" STREQUAL "")
    message("Ignoring invalid .prm '${_test_full}'...")
    set(_use_test OFF)
  endif()
  if (${_test_full} MATCHES "\\.x\\.prm$")
    # Skip files generated by in source builds:
    set(_use_test OFF)
  endif()

  SHOULD_ENABLE_TEST(${CMAKE_CURRENT_SOURCE_DIR}/${_testname}.prm)

  if(_use_test)

    message(STATUS "Processing test ${_testname}:")

    # Create main target for this test. We let it depend on the screen output
    # even if the folder with reference data for this test doesn't contain any
    # files. This is useful when you are constructing a new test.
    add_custom_target(tests.${_testname}
      DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/output-${_testname}/screen-output)

    math(EXPR _n_tests "${_n_tests} + 1")

    # create the target ${_testname} which creates the library
    # lib${_testname}.so if we have a .cc file
    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${_testname}.cc)
      add_library(${_testname} SHARED EXCLUDE_FROM_ALL ${_testname}.cc)
      ASPECT_SETUP_PLUGIN(${_testname})

      set(_testlib 'set Additional shared libraries = ./lib${_testname}.so')
      set(_testdepends ${_testname})
    else()
      set(_testlib)
      set(_testdepends)
    endif()

    # Create the output directory and subdirectories as needed. To do this, we
    # recursively glob all files (and only files, not subdirectories) located
    # in this test's directory.  For each file get the directory part relative
    # to the test directory and create a corresponding subdirectory in the
    # output directory of this test. MAKE_DIRECTORY then does a `mkdir -p`.
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/output-${_testname})

    file(GLOB_RECURSE _outputs RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/${_testname}/
  ${CMAKE_CURRENT_SOURCE_DIR}/${_testname}/[a-zA-Z0-9]*)
    foreach(_output ${_outputs})

      # get the directory name we should work on; for cmake <2.8.11, the
      # subcommand was PATH (legacy), for everything after we can use
      # DIRECTORY
      if (${CMAKE_VERSION} VERSION_LESS "2.8.12")
        get_filename_component(_dir ${_output} PATH)
      else()
        get_filename_component(_dir ${_output} DIRECTORY)
      endif()
      if(NOT  ${_dir} STREQUAL "" )
        file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/output-${_testname}/${_dir})
      endif()
    endforeach()

    # A rule to generate the input file therein.  This input file is copied
    # from the one found in the source dir, but we make sure it has the
    # correct output directory so that we do not need to specify it by hand (a
    # common mistake when copying a previous test) and we specify any
    # additional plugin shared libraries, should they exist for this test
    #
    # We also replace all occurrences of @SOURCE_DIR@ by
    # ${CMAKE_CURRENT_SOURCE_DIR} so that input files can reference
    # mesh, input files, etc, in the tests/ directory without having
    # to know where exactly they are run.
    add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${_testname}.x.prm
      COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/${_testname}.prm
              ${CMAKE_CURRENT_BINARY_DIR}/${_testname}.x.prm
      COMMAND echo ''
              >> ${CMAKE_CURRENT_BINARY_DIR}/${_testname}.x.prm
      COMMAND echo 'set Output directory = output-${_testname}'
              >> ${CMAKE_CURRENT_BINARY_DIR}/${_testname}.x.prm
      COMMAND echo '${_testlib}'
              >> ${CMAKE_CURRENT_BINARY_DIR}/${_testname}.x.prm
      COMMAND ${PERL_EXECUTABLE} -pi
              -e 's!\\@SOURCE_DIR\\@!${CMAKE_CURRENT_SOURCE_DIR}!g;'
              ${CMAKE_CURRENT_BINARY_DIR}/${_testname}.x.prm
      DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${_testname}.prm
      )

    # then generate a rule that runs aspect and normalizes the output
    # files. Before running aspect, delete prior content of the directory to
    # make sure no dead files are left there (this is one way to trip up the
    # 'terminate_user_request' test of Aspect which terminates the program when
    # a certain file appears). we have to take care of not deleting those
    # files that have been placed there on purpose, however, which are all
    # of the .cmp.notime files.
    #
    # the actual run command is a bit complicated because we have to figure out
    # whether we want the test to run in parallel using MPI or not
    GET_MPI_COUNT(${CMAKE_CURRENT_SOURCE_DIR}/${_testname}.prm)

    # detect the optional script <_testname>.sh that a test can use to process
    # output (currently screen-output only).
    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${_testname}.sh)
      set(_replacement_script "${CMAKE_CURRENT_SOURCE_DIR}/${_testname}.sh")
      set(_testdepends ${_testdepends} ${CMAKE_CURRENT_SOURCE_DIR}/${_testname}.sh)
    else()
      set(_replacement_script "${CMAKE_CURRENT_SOURCE_DIR}/cmake/default")
    endif()

    # If a .prm contains the keyword "EXPECT FAILURE", make sure aspect fails
    # with a non-zero return value. Otherwise make sure aspect terminates
    # without errors. The logic for this is in tests/cmake/run_test.sh,
    # because we can not do this in ADD_CUSTOM_COMMAND directly.
    file(STRINGS ${CMAKE_CURRENT_SOURCE_DIR}/${_testname}.prm _input_lines REGEX "EXPECT FAILURE")
    if("${_input_lines}" STREQUAL "")
      set(EXPECT "0")
    else()
      set(EXPECT "1")
      message(STATUS "Test ${_testname} is expected to fail.")
      if (NOT "${_mpi_count}" STREQUAL "1")
        message(FATAL_ERROR "Invalid setup in test '${_testname}.prm': Tests with 'EXPECT FAILURE' are only supported if they use a single MPI rank.")
      endif()
    endif()

    file(GLOB_RECURSE _relative_output_filenames RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/${_testname}/
         ${CMAKE_CURRENT_SOURCE_DIR}/${_testname}/[a-zA-Z0-9]*)
    set(_output_files "")

    foreach(_name ${_relative_output_filenames})
      set(_path_and_name ${CMAKE_CURRENT_BINARY_DIR}/output-${_testname}/${_name})
      set(_output_files ${_output_files} ${_path_and_name})
    endforeach()

    # fail is the user forgot to create reference files instead of silently passing the test
    if("${_output_files}" STREQUAL "")
      message(FATAL_ERROR "test ${_testname}.prm is missing a folder with reference data")
    endif()

    # Create the run command.
    # Note that we are using the timeout.sh script (which is a wrapper around
    # the bash timeout command not available on OSX). This might seem
    # redundant as ctest has an option to limit test time (see below), but we
    # "prebuild" tests on the CI system, which does not use the ctest
    # mechanism. So without this additional timeout command, tests would hang
    # indefinitely.
    # We are creating an empty file "/runs" if the test runs successfully without
    # crashing. This is used later to only diff if necessary.
    if("${_mpi_count}" STREQUAL "1")
      add_custom_command(
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/output-${_testname}/screen-output
        COMMAND rm -f ${CMAKE_CURRENT_BINARY_DIR}/output-${_testname}/runs
        COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/cmake/timeout.sh ${TEST_TIME_LIMIT}s ${CMAKE_CURRENT_SOURCE_DIR}/cmake/run_test.sh ${EXPECT} ${ASPECT_BINARY} ${CMAKE_CURRENT_BINARY_DIR}/${_testname}.x.prm ${CMAKE_CURRENT_BINARY_DIR}/output-${_testname}/screen-output.tmp
        COMMAND ${_replacement_script} screen-output <${CMAKE_CURRENT_BINARY_DIR}/output-${_testname}/screen-output.tmp >${CMAKE_CURRENT_BINARY_DIR}/output-${_testname}/screen-output
        COMMAND touch ${CMAKE_CURRENT_BINARY_DIR}/output-${_testname}/runs
        DEPENDS ${ASPECT_BINARY} ${CMAKE_CURRENT_BINARY_DIR}/${_testname}.x.prm ${_testdepends})
    else()
      add_custom_command(
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/output-${_testname}/screen-output
        COMMAND rm -f ${CMAKE_CURRENT_BINARY_DIR}/output-${_testname}/runs
        COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/cmake/timeout.sh ${TEST_TIME_LIMIT}s mpirun -np ${_mpi_count} ${ASPECT_BINARY} ${CMAKE_CURRENT_BINARY_DIR}/${_testname}.x.prm
                > ${CMAKE_CURRENT_BINARY_DIR}/output-${_testname}/screen-output.tmp
        COMMAND ${_replacement_script} screen-output <${CMAKE_CURRENT_BINARY_DIR}/output-${_testname}/screen-output.tmp >${CMAKE_CURRENT_BINARY_DIR}/output-${_testname}/screen-output
        COMMAND touch ${CMAKE_CURRENT_BINARY_DIR}/output-${_testname}/runs
        DEPENDS ${ASPECT_BINARY} ${CMAKE_CURRENT_BINARY_DIR}/${_testname}.x.prm ${_testdepends})
    endif()


    get_filename_component(ASPECT_BASE_DIR ${CMAKE_CURRENT_SOURCE_DIR} PATH)

    if(ASPECT_COMPARE_TEST_RESULTS)
      # Commands to normalize each output. Note that we depend on screen-output
      # instead of ${_output_name} here, to avoid having to specify BYPRODUCTS
      # for ninja.
      foreach(_output_name ${_output_files})
        add_custom_command(OUTPUT ${_output_name}.notime
          COMMAND ${PERL_EXECUTABLE}
                  ${CMAKE_CURRENT_SOURCE_DIR}/normalize.pl
                  ${_output_name} ${ASPECT_BASE_DIR}
                  >${_output_name}.notime
          DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/output-${_testname}/screen-output)
      endforeach()

      # Create commands to copy and normalize reference output
      foreach(_name ${_relative_output_filenames})
        add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/output-${_testname}/${_name}.cmp.notime
          COMMAND ${PERL_EXECUTABLE}
                  ${CMAKE_CURRENT_SOURCE_DIR}/normalize.pl
                  ${CMAKE_CURRENT_SOURCE_DIR}/${_testname}/${_name} ${ASPECT_BASE_DIR}
                  >${CMAKE_CURRENT_BINARY_DIR}/output-${_testname}/${_name}.cmp.notime
          DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${_testname}/${_name})
      endforeach()

      # Now create commands to diff reference with output file:
      foreach(_name ${_relative_output_filenames})
        # the normalized reference file:
        set(_ref_file ${CMAKE_CURRENT_BINARY_DIR}/output-${_testname}/${_name}.cmp.notime)
        # the normalized output file:
        set(_run_file ${CMAKE_CURRENT_BINARY_DIR}/output-${_testname}/${_name}.notime)

        # generate .diff file, only run diff if "/runs" exists (test ran successfully),
        # otherwise GENERATE_REFERENCE_OUTPUT would overwrite reference data with garbage.
        add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/output-${_testname}/${_name}.diff
          COMMAND test -f ${CMAKE_CURRENT_BINARY_DIR}/output-${_testname}/runs
          COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/cmake/diff_test.sh
                  ${TEST_DIFF}
                  ${_testname}/${_name}
                  ${CMAKE_CURRENT_SOURCE_DIR}
                  ${CMAKE_CURRENT_BINARY_DIR}
                  ${ASPECT_COMPARE_TEST_RESULTS}
          DEPENDS ${_ref_file} ${_run_file})

        # Add the target for this output file to the dependencies of this
        # test. Note that a target may not contain "/" but outputs can be
        # in subdirectories (for example solution/solution...), so replace them.
        string(REPLACE "/" "-" _target_name ${_testname}.${_name}.diff)
        add_custom_target(${_target_name}
          DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/output-${_testname}/${_name}.diff)
        add_dependencies(tests.${_testname} ${_target_name})

      endforeach()
    endif()

    # add this test to the dependencies of the overall 'tests' target
    # and declare it to ctest
    add_dependencies(tests tests.${_testname})
    add_test(NAME ${_testname}
      COMMAND ${CMAKE_COMMAND}
              -DBINARY_DIR=${CMAKE_BINARY_DIR}
              -DTESTNAME=tests.${_testname}
              -DERROR="Test ${_testname} failed"
              -P ${CMAKE_SOURCE_DIR}/run_test.cmake
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set_tests_properties(${_testname} PROPERTIES TIMEOUT ${TEST_TIME_LIMIT})

    # Check if this test depends on another test
    GET_DEPENDS_ON(${CMAKE_CURRENT_SOURCE_DIR}/${_testname}.prm)
    if(NOT "${_depends_on}" STREQUAL "")
      message(STATUS "Test dependency detected for test:" ${_testname} ", depending on:" ${_depends_on})
      set_tests_properties(${_testname} PROPERTIES DEPENDS ${_depends_on})
    endif()
  endif()
endforeach()

message(STATUS "Found ${_n_tests} tests.")
