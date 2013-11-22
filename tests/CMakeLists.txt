################### top matter ################
FIND_PACKAGE(Perl)

MACRO(SET_IF_EMPTY _variable _value)
  IF("${${_variable}}" STREQUAL "")
    SET(${_variable} ${_value})
  ENDIF()
ENDMACRO()

# set a time limit of 10 minutes per test. this should be long
# enough even for long-running tests, and short enough to not
# get into trouble should we have an infinite loop.
SET_IF_EMPTY(TEST_TIME_LIMIT 600)

############################3

ADD_CUSTOM_TARGET(tests)

SET_IF_EMPTY(DEAL_II_DIFF "$ENV{DEAL_II_DIFF}")
SET_IF_EMPTY(DEAL_II_DIFF "numdiff -a 1e-6 -q -s \" \\t\\n:\"")


FILE(GLOB _tests *.prm)
LIST(SORT _tests)
FOREACH(_test ${_tests})
  GET_FILENAME_COMPONENT(_test ${_test} NAME_WE)
  MESSAGE(STATUS "Adding test ${_test}")

  IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${_test}.cc)
    ADD_LIBRARY(${_test} SHARED EXCLUDE_FROM_ALL ${_test}.cc)
    DEAL_II_SETUP_TARGET(${_test})
    SET(_testdepends ${_test})
  ELSE()
    SET(_testdepends)
  ENDIF()

  # Create the output directory, then generate a rule that runs aspect and normalizes
  # the output files. Before running aspect, delete prior content of the directory
  # to make sure no dead files are left there (this is one way to trip up the
  # 'terminate_user_request' test which terminates the program when a certain
  # file appears)
  FILE(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/output-${_test})
  ADD_CUSTOM_COMMAND(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/output-${_test}/screen-output
    COMMAND rm ${CMAKE_CURRENT_BINARY_DIR}/output-${_test}/*
    COMMAND aspect ${CMAKE_CURRENT_SOURCE_DIR}/${_test}.prm
            > ${CMAKE_CURRENT_BINARY_DIR}/output-${_test}/screen-output 2>&1
    COMMAND ${PERL_EXECUTABLE} -pi
              ${CMAKE_CURRENT_SOURCE_DIR}/normalize.pl
              ${CMAKE_CURRENT_BINARY_DIR}/output-${_test}/*
    DEPENDS aspect ${CMAKE_CURRENT_SOURCE_DIR}/${_test}.prm ${_testdepends}
    )

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

    # create the output.notime target
    ADD_CUSTOM_COMMAND(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/output-${_test}/${_output}.notime
      COMMAND
	cat ${CMAKE_CURRENT_BINARY_DIR}/output-${_test}/${_output}
        | egrep -v '^\\|'
        > ${CMAKE_CURRENT_BINARY_DIR}/output-${_test}/${_output}.notime
      DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/output-${_test}/screen-output
    )

    # create the target that compares the .notime with the saved file
    ADD_CUSTOM_COMMAND(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/output-${_test}/${_output}.diff
      COMMAND
	if (diff ${CMAKE_CURRENT_SOURCE_DIR}/${_test}/${_output}
              ${CMAKE_CURRENT_BINARY_DIR}/output-${_test}/${_output}.notime
            > ${CMAKE_CURRENT_BINARY_DIR}/output-${_test}/${_output}.diff) \; then
          : \;
        else
	  mv ${CMAKE_CURRENT_BINARY_DIR}/output-${_test}/${_output}.diff
             ${CMAKE_CURRENT_BINARY_DIR}/output-${_test}/${_output}.diff.failed \;
          echo "******* Error during diffing output results for ${_test}/${_output}" \;
          echo "******* Results are stored in ${CMAKE_CURRENT_BINARY_DIR}/output-${_test}/${_output}.diff.failed" \;
          echo "******* Diffs are:" \;
          cat ${CMAKE_CURRENT_BINARY_DIR}/output-${_test}/${_output}.diff.failed \;
	  false \;
        fi
      DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${_test}/${_output}
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
ENDFOREACH()

