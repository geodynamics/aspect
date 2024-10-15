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
CMAKE_MINIMUM_REQUIRED(VERSION 3.13.4)

IF (CMP0058)
 cmake_policy(CMP0058 NEW)
ENDIF()

PROJECT(testsuite CXX)

SET(ASPECT_BINARY "${ASPECT_BINARY}" CACHE STRING "" FORCE)
SET(Aspect_DIR "${Aspect_DIR}" CACHE STRING "" FORCE)
SET(ASPECT_COMPARE_TEST_RESULTS "${ASPECT_COMPARE_TEST_RESULTS}" CACHE BOOL "" FORCE)
SET(ASPECT_RUN_ALL_TESTS OFF CACHE BOOL "")

FIND_PACKAGE(Aspect 2.4.0 QUIET REQUIRED HINTS ${Aspect_DIR})
DEAL_II_INITIALIZE_CACHED_VARIABLES()

FIND_PACKAGE(Perl)

MESSAGE(STATUS "Setting up tests with CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE}")


ENABLE_TESTING()

MACRO(SET_IF_EMPTY _variable _value)
  IF("${${_variable}}" STREQUAL "")
    SET(${_variable} ${_value})
  ENDIF()
ENDMACRO()

# suppress mac os warning
IF(APPLE AND NOT DEFINED CMAKE_MACOSX_RPATH)
  SET(CMAKE_MACOSX_RPATH 0) # see policy CMP0042
ENDIF()

# A function that extracts from a file (presumably a .prm file)
# the number of MPI processes this test is to be invoked as.
# This is encoded in .prm files through lines of the form
#    '# MPI: 4'
# The result is returned in a variable _mpi_count in the
# caller's scope.
FUNCTION(get_mpi_count _filename)
  FILE(STRINGS ${_filename} _input_lines
       REGEX "MPI:")
  IF("${_input_lines}" STREQUAL "")
    SET(_mpi_count 1 PARENT_SCOPE)
  ELSE()
    # go over the (possibly multiple) lines with MPI markers and choose the last
    FOREACH(_input_line ${_input_lines})
     SET(_last_line ${_input_line})
    ENDFOREACH()
    STRING(REGEX REPLACE "^ *# *MPI: *([0-9]+) *$" "\\1"
           _mpi_count ${_last_line})
    SET(_mpi_count "${_mpi_count}" PARENT_SCOPE)
  ENDIF()
ENDFUNCTION()

# A function that checks if a test depends on another test.
# This is encoded in .prm files through lines of the form
#    '# DEPENDS-ON: testname'
# The result is returned in a variable _depends_on in the
# caller's scope. Having test dependencies is helpful if
# one test requires another test to finish first, for example
# because the result of the first test is used by the later test.
FUNCTION(get_depends_on _filename)
  FILE(STRINGS ${_filename} _input_lines
       REGEX "DEPENDS-ON:")
  IF("${_input_lines}" STREQUAL "")
    SET(_depends_on "" PARENT_SCOPE)
  ELSE()
    # go over the (possibly multiple) lines with DEPENDS-ON markers and choose the last
    FOREACH(_input_line ${_input_lines})
     SET(_last_line ${_input_line})
    ENDFOREACH()
    STRING(REGEX REPLACE "^ *# *DEPENDS-ON: *(.*) *$" "\\1"
           _depends_on ${_last_line})
    SET(_depends_on "${_depends_on}" PARENT_SCOPE)
  ENDIF()
ENDFUNCTION()

# Analyze the .prm to decide if this test should be run or not. We are
# checking for the following strings in the .prm:
# 1. "Enable if: ASPECT_WITH_WORLD_BUILDER" - only run with World Builder.
# 2. "Enable if: QUICK_TEST" - enable the test even if RUN_ALL_TESTS is false
FUNCTION(SHOULD_ENABLE_TEST _filename)

  FILE(STRINGS ${_filename} _input_lines
       REGEX "Enable if: ASPECT_WITH_WORLD_BUILDER")
  IF(NOT "${_input_lines}" STREQUAL "")
    IF (${_input_lines} MATCHES "Enable if: ASPECT_WITH_WORLD_BUILDER"
        AND NOT ASPECT_WITH_WORLD_BUILDER)
      SET(_use_test OFF PARENT_SCOPE)
    ENDIF()
  ENDIF()

  FILE(STRINGS ${_filename} _input_lines
       REGEX "Enable if: QUICK_TEST")
  IF(NOT ASPECT_RUN_ALL_TESTS AND "${_input_lines}" STREQUAL "")
    SET(_use_test OFF PARENT_SCOPE)
  ENDIF()
ENDFUNCTION()


# set a time limit of 10 minutes per test. this should be long
# enough even for long-running tests, and short enough to not
# get into trouble should we have an infinite loop.
SET_IF_EMPTY(TEST_TIME_LIMIT 600)

############################3

ADD_CUSTOM_TARGET(tests)

#
# We need a diff tool, preferably numdiff otherwise diff. Let the user
# override it by specifying TEST_DIFF.
#
FIND_PROGRAM(DIFF_EXECUTABLE
  NAMES diff
  HINTS ${DIFF_DIR}
  PATH_SUFFIXES bin
  )

FIND_PROGRAM(NUMDIFF_EXECUTABLE
  NAMES numdiff
  HINTS ${NUMDIFF_DIR}
  PATH_SUFFIXES bin
  )

MARK_AS_ADVANCED(DIFF_EXECUTABLE NUMDIFF_EXECUTABLE)

IF("${TEST_DIFF}" STREQUAL "")
  #

  IF(NOT NUMDIFF_EXECUTABLE MATCHES "-NOTFOUND")
    SET(TEST_DIFF ${NUMDIFF_EXECUTABLE})
  ELSEIF(NOT DIFF_EXECUTABLE MATCHES "-NOTFOUND")
    SET(TEST_DIFF ${DIFF_EXECUTABLE})
  ELSE()
    MESSAGE(FATAL_ERROR
      "Could not find diff or numdiff. One of those are required for running the testsuite.\n"
      "Please specify TEST_DIFF by hand."
      )
  ENDIF()
ENDIF()

FILE(GLOB _tests *.prm)

SET(_n_tests "0")
LIST(SORT _tests)

FOREACH(_test ${_tests})

  SET(_test_full ${_test})
  GET_FILENAME_COMPONENT(_testname ${_test} NAME_WE)

  SET(_use_test ON)

  IF ("${_testname}" STREQUAL "")
    MESSAGE("Ignoring invalid .prm '${_test_full}'...")
    SET(_use_test OFF)
  ENDIF()
  IF (${_test_full} MATCHES "\\.x\\.prm$")
    # Skip files generated by in source builds:
    SET(_use_test OFF)
  ENDIF()

  SHOULD_ENABLE_TEST(${CMAKE_CURRENT_SOURCE_DIR}/${_testname}.prm)

  IF(_use_test)

    MESSAGE(STATUS "** processing test ${_testname}:")

    # Create main target for this test. We let it depend on the screen output
    # even if the folder with reference data for this test doesn't contain any
    # files. This is useful when you are constructing a new test.
    ADD_CUSTOM_TARGET(tests.${_testname}
      DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/output-${_testname}/screen-output)

    MATH(EXPR _n_tests "${_n_tests} + 1")

    # create the target ${_testname} which creates the library
    # lib${_testname}.so if we have a .cc file
    IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${_testname}.cc)
      ADD_LIBRARY(${_testname} SHARED EXCLUDE_FROM_ALL ${_testname}.cc)
      ASPECT_SETUP_PLUGIN(${_testname})

      SET(_testlib 'set Additional shared libraries = ./lib${_testname}.so')
      SET(_testdepends ${_testname})
    ELSE()
      SET(_testlib)
      SET(_testdepends)
    ENDIF()

    # Create the output directory and subdirectories as needed. To do this, we
    # recursively glob all files (and only files, not subdirectories) located
    # in this test's directory.  For each file get the directory part relative
    # to the test directory and create a corresponding subdirectory in the
    # output directory of this test. MAKE_DIRECTORY then does a `mkdir -p`.
    FILE(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/output-${_testname})

    FILE(GLOB_RECURSE _outputs RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/${_testname}/
  ${CMAKE_CURRENT_SOURCE_DIR}/${_testname}/[a-zA-Z0-9]*)
    FOREACH(_output ${_outputs})

      # get the directory name we should work on; for cmake <2.8.11, the
      # subcommand was PATH (legacy), for everything after we can use
      # DIRECTORY
      IF (${CMAKE_VERSION} VERSION_LESS "2.8.12")
        GET_FILENAME_COMPONENT(_dir ${_output} PATH)
      ELSE()
        GET_FILENAME_COMPONENT(_dir ${_output} DIRECTORY)
      ENDIF()
      IF(NOT  ${_dir} STREQUAL "" )
        FILE(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/output-${_testname}/${_dir})
      ENDIF()
    ENDFOREACH()

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
    ADD_CUSTOM_COMMAND(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${_testname}.x.prm
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
    IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${_testname}.sh)
      SET(_replacement_script "${CMAKE_CURRENT_SOURCE_DIR}/${_testname}.sh")
      SET(_testdepends ${_testdepends} ${CMAKE_CURRENT_SOURCE_DIR}/${_testname}.sh)
    ELSE()
      SET(_replacement_script "${CMAKE_CURRENT_SOURCE_DIR}/cmake/default")
    ENDIF()

    # If a .prm contains the keyword "EXPECT FAILURE", make sure aspect fails
    # with a non-zero return value. Otherwise make sure aspect terminates
    # without errors. The logic for this is in tests/cmake/run_test.sh,
    # because we can not do this in ADD_CUSTOM_COMMAND directly.
    FILE(STRINGS ${CMAKE_CURRENT_SOURCE_DIR}/${_testname}.prm _input_lines REGEX "EXPECT FAILURE")
    IF("${_input_lines}" STREQUAL "")
      SET(EXPECT "0")
    ELSE()
      SET(EXPECT "1")
      MESSAGE(STATUS "Test ${_testname} is expected to fail.")
      IF (NOT "${_mpi_count}" STREQUAL "1")
        MESSAGE(FATAL_ERROR "Invalid setup in test '${_testname}.prm': Tests with 'EXPECT FAILURE' are only supported if they use a single MPI rank.")
      ENDIF()
    ENDIF()

    FILE(GLOB_RECURSE _relative_output_filenames RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/${_testname}/
         ${CMAKE_CURRENT_SOURCE_DIR}/${_testname}/[a-zA-Z0-9]*)
    SET(_output_files "")

    FOREACH(_name ${_relative_output_filenames})
      SET(_path_and_name ${CMAKE_CURRENT_BINARY_DIR}/output-${_testname}/${_name})
      SET(_output_files ${_output_files} ${_path_and_name})
    ENDFOREACH()

    # fail is the user forgot to create reference files instead of silently passing the test
    IF("${_output_files}" STREQUAL "")
      MESSAGE(FATAL_ERROR "test ${_testname}.prm is missing a folder with reference data")
    ENDIF()

    # Create the run command.
    # Note that we are using the timeout.sh script (which is a wrapper around
    # the bash timeout command not available on OSX). This might seem
    # redundant as ctest has an option to limit test time (see below), but we
    # "prebuild" tests on the CI system, which does not use the ctest
    # mechanism. So without this additional timeout command, tests would hang
    # indefinitely.
    # We are creating an empty file "/runs" if the test runs successfully without
    # crashing. This is used later to only diff if necessary.
    IF("${_mpi_count}" STREQUAL "1")
      ADD_CUSTOM_COMMAND(
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/output-${_testname}/screen-output
        COMMAND rm -f ${CMAKE_CURRENT_BINARY_DIR}/output-${_testname}/runs
        COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/cmake/timeout.sh ${TEST_TIME_LIMIT}s ${CMAKE_CURRENT_SOURCE_DIR}/cmake/run_test.sh ${EXPECT} ${ASPECT_BINARY} ${CMAKE_CURRENT_BINARY_DIR}/${_testname}.x.prm ${CMAKE_CURRENT_BINARY_DIR}/output-${_testname}/screen-output.tmp
        COMMAND ${_replacement_script} screen-output <${CMAKE_CURRENT_BINARY_DIR}/output-${_testname}/screen-output.tmp >${CMAKE_CURRENT_BINARY_DIR}/output-${_testname}/screen-output
        COMMAND touch ${CMAKE_CURRENT_BINARY_DIR}/output-${_testname}/runs
        DEPENDS ${ASPECT_BINARY} ${CMAKE_CURRENT_BINARY_DIR}/${_testname}.x.prm ${_testdepends})
    ELSE()
      ADD_CUSTOM_COMMAND(
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/output-${_testname}/screen-output
        COMMAND rm -f ${CMAKE_CURRENT_BINARY_DIR}/output-${_testname}/runs
        COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/cmake/timeout.sh ${TEST_TIME_LIMIT}s mpirun -np ${_mpi_count} ${ASPECT_BINARY} ${CMAKE_CURRENT_BINARY_DIR}/${_testname}.x.prm
                > ${CMAKE_CURRENT_BINARY_DIR}/output-${_testname}/screen-output.tmp
        COMMAND ${_replacement_script} screen-output <${CMAKE_CURRENT_BINARY_DIR}/output-${_testname}/screen-output.tmp >${CMAKE_CURRENT_BINARY_DIR}/output-${_testname}/screen-output
        COMMAND touch ${CMAKE_CURRENT_BINARY_DIR}/output-${_testname}/runs
        DEPENDS ${ASPECT_BINARY} ${CMAKE_CURRENT_BINARY_DIR}/${_testname}.x.prm ${_testdepends})
    ENDIF()


    GET_FILENAME_COMPONENT(ASPECT_BASE_DIR ${CMAKE_CURRENT_SOURCE_DIR} PATH)

    IF(ASPECT_COMPARE_TEST_RESULTS)
      # Commands to normalize each output. Note that we depend on screen-output
      # instead of ${_output_name} here, to avoid having to specify BYPRODUCTS
      # for ninja.
      FOREACH(_output_name ${_output_files})
        ADD_CUSTOM_COMMAND(OUTPUT ${_output_name}.notime
          COMMAND ${PERL_EXECUTABLE}
                  ${CMAKE_CURRENT_SOURCE_DIR}/normalize.pl
                  ${_output_name} ${ASPECT_BASE_DIR}
                  >${_output_name}.notime
          DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/output-${_testname}/screen-output)
      ENDFOREACH()

      # Create commands to copy and normalize reference output
      FOREACH(_name ${_relative_output_filenames})
        ADD_CUSTOM_COMMAND(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/output-${_testname}/${_name}.cmp.notime
          COMMAND ${PERL_EXECUTABLE}
                  ${CMAKE_CURRENT_SOURCE_DIR}/normalize.pl
                  ${CMAKE_CURRENT_SOURCE_DIR}/${_testname}/${_name} ${ASPECT_BASE_DIR}
                  >${CMAKE_CURRENT_BINARY_DIR}/output-${_testname}/${_name}.cmp.notime
          DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${_testname}/${_name})
      ENDFOREACH()

      # Now create commands to diff reference with output file:
      FOREACH(_name ${_relative_output_filenames})
        # the normalized reference file:
        SET(_ref_file ${CMAKE_CURRENT_BINARY_DIR}/output-${_testname}/${_name}.cmp.notime)
        # the normalized output file:
        SET(_run_file ${CMAKE_CURRENT_BINARY_DIR}/output-${_testname}/${_name}.notime)

        # generate .diff file, only run diff if "/runs" exists (test ran successfully),
        # otherwise GENERATE_REFERENCE_OUTPUT would overwrite reference data with garbage.
        ADD_CUSTOM_COMMAND(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/output-${_testname}/${_name}.diff
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
        STRING(REPLACE "/" "-" _target_name ${_testname}.${_name}.diff)
        ADD_CUSTOM_TARGET(${_target_name}
          DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/output-${_testname}/${_name}.diff)
        ADD_DEPENDENCIES(tests.${_testname} ${_target_name})

      ENDFOREACH()
    ENDIF()

    # add this test to the dependencies of the overall 'tests' target
    # and declare it to ctest
    ADD_DEPENDENCIES(tests tests.${_testname})
    ADD_TEST(NAME ${_testname}
      COMMAND ${CMAKE_COMMAND}
              -DBINARY_DIR=${CMAKE_BINARY_DIR}
              -DTESTNAME=tests.${_testname}
              -DERROR="Test ${_testname} failed"
              -P ${CMAKE_SOURCE_DIR}/run_test.cmake
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    SET_TESTS_PROPERTIES(${_testname} PROPERTIES TIMEOUT ${TEST_TIME_LIMIT})

    # Check if this test depends on another test
    GET_DEPENDS_ON(${CMAKE_CURRENT_SOURCE_DIR}/${_testname}.prm)
    IF(NOT "${_depends_on}" STREQUAL "")
      MESSAGE(STATUS "Test dependency detected for test:" ${_testname} ", depending on:" ${_depends_on})
      SET_TESTS_PROPERTIES(${_testname} PROPERTIES DEPENDS ${_depends_on})
    ENDIF()
  ENDIF()
ENDFOREACH()

MESSAGE("     Found ${_n_tests} tests.")
