################### top matter ################
FIND_PACKAGE(Perl)

MACRO(SET_IF_EMPTY _variable _value)
  IF("${${_variable}}" STREQUAL "")
    SET(${_variable} ${_value})
  ENDIF()
ENDMACRO()

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
  endif()
ENDFUNCTION()

# extract from the prm if the test should be run with trilinos or
# petsc. This is encoded in the .prm in form of line
# '# LINEAR ALGEBRA: PETSC'
# '# LINEAR ALGEBRA: PETSC TRILINOS'
# If no line is there, Trilinos is assumed. 
FUNCTION(check_linear_algebra _filename)
  FILE(STRINGS ${_filename} _input_lines
       REGEX "LINEAR ALGEBRA:")
  IF("${_input_lines}" STREQUAL "")
    IF(ASPECT_USE_PETSC)
      SET(_use_test 0 PARENT_SCOPE)
    ELSE()
      SET(_use_test 1 PARENT_SCOPE)
    ENDIF()
  ELSE()
    SET(_use_test 0 PARENT_SCOPE)
    FOREACH(_input_line ${_input_lines})
      IF (${_input_line} MATCHES "PETSC" AND ASPECT_USE_PETSC)
        SET(_use_test 1 PARENT_SCOPE)
      ENDIF()
      IF (${_input_line} MATCHES "TRILINOS" AND NOT ASPECT_USE_PETSC)
        SET(_use_test 1 PARENT_SCOPE)
      ENDIF()
    ENDFOREACH()
  ENDIF()
ENDFUNCTION()


# set a time limit of 10 minutes per test. this should be long
# enough even for long-running tests, and short enough to not
# get into trouble should we have an infinite loop.
SET_IF_EMPTY(TEST_TIME_LIMIT 600)

############################3

ADD_CUSTOM_TARGET(tests)

#
# We need a diff tool, preferably numdiff:
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
  # No TEST_DIFF is set, specify one:
  #

  IF(NOT NUMDIFF_EXECUTABLE MATCHES "-NOTFOUND")
    SET(TEST_DIFF ${NUMDIFF_EXECUTABLE} -a 1e-6 -r 1e-8 -s ' \\t\\n:,')
    IF(DIFF_EXECUTABLE MATCHES "-NOTFOUND")
      SET(DIFF_EXECUTABLE ${NUMDIFF_EXECUTABLE})
    ENDIF()
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
  GET_FILENAME_COMPONENT(_test ${_test} NAME_WE)

  CHECK_LINEAR_ALGEBRA(${CMAKE_CURRENT_SOURCE_DIR}/${_test}.prm)

  IF("${_use_test}" STREQUAL "1")

    MATH(EXPR _n_tests "${_n_tests} + 1")

    IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${_test}.cc)
      ADD_LIBRARY(${_test} SHARED EXCLUDE_FROM_ALL ${_test}.cc)
      IF (ASPECT_USE_PETSC)
        SET_PROPERTY(SOURCE ${_test}.cc
	APPEND PROPERTY COMPILE_DEFINITIONS ASPECT_USE_PETSC="1")
      ENDIF()
      DEAL_II_SETUP_TARGET(${_test})
      SET(_testlib 'set Additional shared libraries = ./lib${_test}.so')
      SET(_testdepends ${_test})
    ELSE()
      SET(_testdepends)
      SET(_testlib)
    ENDIF()

    # Create the output directory and a rule to generate the input file therein.
    # This input file is the one found in the source dir, but it has the
    # correct output directory so that we do not need to specify it by hand
    # (a common mistake when copying a previous test) and do not need to
    # specify any additional plugin shared libraries, should they exist
    # for this test
    #
    # We also replace all occurrences of @SOURCE_DIR@ by
    # ${CMAKE_CURRENT_SOURCE_DIR} so that input files can reference
    # mesh, input files, etc, in the tests/ directory without having
    # to know where exactly they are run.
    FILE(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/output-${_test})
    ADD_CUSTOM_COMMAND(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${_test}.x.prm
      COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/${_test}.prm
		 ${CMAKE_CURRENT_BINARY_DIR}/${_test}.x.prm
      COMMAND echo 'set Output directory = output-${_test}'
		 >> ${CMAKE_CURRENT_BINARY_DIR}/${_test}.x.prm
      COMMAND echo '${_testlib}'
		 >> ${CMAKE_CURRENT_BINARY_DIR}/${_test}.x.prm
      COMMAND ${PERL_EXECUTABLE} -pi
		 -e 's!\\@SOURCE_DIR\\@!${CMAKE_CURRENT_SOURCE_DIR}!g;'
		 ${CMAKE_CURRENT_BINARY_DIR}/${_test}.x.prm
      DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${_test}.prm
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
    GET_MPI_COUNT(${CMAKE_CURRENT_SOURCE_DIR}/${_test}.prm)

    # detect the optional script <testname>.sh that a test can use to process
    # output (currently screen-output only).
    IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${_test}.sh)
        SET(_replacement_script "${CMAKE_CURRENT_SOURCE_DIR}/${_test}.sh")
        SET(_testdepends ${_testdepends} ${CMAKE_CURRENT_SOURCE_DIR}/${_test}.sh)
    ELSE()
        SET(_replacement_script "${CMAKE_CURRENT_SOURCE_DIR}/default.sh")
    ENDIF()

    IF("${_mpi_count}" STREQUAL "1")
      ADD_CUSTOM_COMMAND(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/output-${_test}/screen-output
	COMMAND
	  for i in ${CMAKE_CURRENT_BINARY_DIR}/output-${_test}/* \; do
	    if echo \$i | grep -q -v .cmp.notime \; then
	      rm -f \$i \;
	    fi \;
	  done
	COMMAND aspect ${CMAKE_CURRENT_BINARY_DIR}/${_test}.x.prm
		> ${CMAKE_CURRENT_BINARY_DIR}/output-${_test}/screen-output.tmp
	COMMAND ${_replacement_script} screen-output <${CMAKE_CURRENT_BINARY_DIR}/output-${_test}/screen-output.tmp >${CMAKE_CURRENT_BINARY_DIR}/output-${_test}/screen-output
	COMMAND ${PERL_EXECUTABLE} -pi
		  ${CMAKE_CURRENT_SOURCE_DIR}/normalize.pl
		  ${CMAKE_CURRENT_BINARY_DIR}/output-${_test}/*
	DEPENDS aspect ${CMAKE_CURRENT_BINARY_DIR}/${_test}.x.prm ${_testdepends}
	)
    ELSE()
      ADD_CUSTOM_COMMAND(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/output-${_test}/screen-output
	COMMAND
	  for i in ${CMAKE_CURRENT_BINARY_DIR}/output-${_test}/* \; do
	    if echo \$i | grep -q -v .cmp.notime \; then
	      rm -f \$i \;
	    fi \;
	  done
	COMMAND mpirun -np ${_mpi_count} ${CMAKE_CURRENT_BINARY_DIR}/../aspect ${CMAKE_CURRENT_BINARY_DIR}/${_test}.x.prm
		> ${CMAKE_CURRENT_BINARY_DIR}/output-${_test}/screen-output.tmp
	COMMAND ${_replacement_script} screen-output <${CMAKE_CURRENT_BINARY_DIR}/output-${_test}/screen-output.tmp >${CMAKE_CURRENT_BINARY_DIR}/output-${_test}/screen-output
	COMMAND ${PERL_EXECUTABLE} -pi
		  ${CMAKE_CURRENT_SOURCE_DIR}/normalize.pl
		  ${CMAKE_CURRENT_BINARY_DIR}/output-${_test}/*
	DEPENDS aspect ${CMAKE_CURRENT_BINARY_DIR}/${_test}.x.prm ${_testdepends}
	)
    ENDIF()


    # The final target for this test
    ADD_CUSTOM_TARGET(${_test}.run
      DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/output-${_test}/screen-output
      )


    ADD_CUSTOM_TARGET(tests.${_test})


    # for each of the output files saved in the source directory for
    # this test, create targets that
    # - create a 'output.notime' file
    # - diff this file with the saved one and save the result
    FILE(GLOB _outputs ${CMAKE_CURRENT_SOURCE_DIR}/${_test}/[a-zA-Z0-9]*)
    FOREACH(_output ${_outputs})
      GET_FILENAME_COMPONENT(_output ${_output} NAME)

      get_filename_component(ASPECT_BASE_DIR ${CMAKE_CURRENT_SOURCE_DIR} PATH)

      # create the output.notime targets for each of the output files,
      # and also create a similar file for all of the input files (this
      # avoids the need to strip timing information from files we may
      # copy from tests we have just run to the directory where we save
      # these files
      #
      # do the same with the second egrep command to strip out the status
      # message at the top of the run
      ADD_CUSTOM_COMMAND(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/output-${_test}/${_output}.cmp.notime
	COMMAND
	  cat ${CMAKE_CURRENT_SOURCE_DIR}/${_test}/${_output}
	  | egrep -v '^\\|'
	  | egrep -v '^--'
	  > ${CMAKE_CURRENT_BINARY_DIR}/output-${_test}/${_output}.cmp.notime
	DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${_test}/${_output}
      )

      ADD_CUSTOM_COMMAND(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/output-${_test}/${_output}.notime
	COMMAND
	  cat ${CMAKE_CURRENT_BINARY_DIR}/output-${_test}/${_output}
	  | egrep -v '^\\|.*\\|'
	  | egrep -v '^--'
	  | sed 's\#${ASPECT_BASE_DIR}\#ASPECT_DIR\#g'
	  > ${CMAKE_CURRENT_BINARY_DIR}/output-${_test}/${_output}.notime
	DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/output-${_test}/screen-output
		${CMAKE_CURRENT_BINARY_DIR}/output-${_test}/${_output}.cmp.notime
      )

      # create the target that compares the .notime with the saved file
      ADD_CUSTOM_COMMAND(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/output-${_test}/${_output}.diff
	COMMAND
	  if (${TEST_DIFF} ${CMAKE_CURRENT_BINARY_DIR}/output-${_test}/${_output}.cmp.notime
		${CMAKE_CURRENT_BINARY_DIR}/output-${_test}/${_output}.notime
	      > ${CMAKE_CURRENT_BINARY_DIR}/output-${_test}/${_output}.diff) \; then
	    : \;
	  else
	    mv ${CMAKE_CURRENT_BINARY_DIR}/output-${_test}/${_output}.diff
	       ${CMAKE_CURRENT_BINARY_DIR}/output-${_test}/${_output}.diff.failed \;
	    echo "******* Error during diffing output results for ${_test}/${_output}" \;
	    echo "******* Results are stored in ${CMAKE_CURRENT_BINARY_DIR}/output-${_test}/${_output}.diff.failed" \;
	    echo "******* Check: ${CMAKE_CURRENT_BINARY_DIR}/output-${_test}/${_output} ${CMAKE_CURRENT_SOURCE_DIR}/${_test}/${_output}" \;
	    echo "******* Diffs are:" \;
	    cat ${CMAKE_CURRENT_BINARY_DIR}/output-${_test}/${_output}.diff.failed \;
	    false \;
	  fi
	DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/output-${_test}/${_output}.cmp.notime
		${CMAKE_CURRENT_BINARY_DIR}/output-${_test}/${_output}.notime
      )

      # add the target for this output file to the dependencies of this test
      ADD_CUSTOM_TARGET(${_test}.${_output}.diff
	DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/output-${_test}/${_output}.diff
	)
      ADD_DEPENDENCIES(tests.${_test}
		       ${_test}.${_output}.diff)
    ENDFOREACH()

    # add this test to the dependencies of the overall 'tests' target
    # and declare it to ctest
    ADD_DEPENDENCIES(tests tests.${_test})
    ADD_TEST(NAME ${_test}
      COMMAND
	 ${CMAKE_COMMAND}
	    -DBINARY_DIR=${CMAKE_BINARY_DIR}
	    -DTESTNAME=tests.${_test}
	    -DERROR="Test ${_test} failed"
	    -P ${CMAKE_SOURCE_DIR}/tests/run_test.cmake
	  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      )
    SET_TESTS_PROPERTIES(${_test} PROPERTIES
	  TIMEOUT ${TEST_TIME_LIMIT}
      )

  ENDIF()
ENDFOREACH()

MESSAGE(STATUS "Added ${_n_tests} tests")
