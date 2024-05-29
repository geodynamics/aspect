# arguments checking
if( NOT TEST_NAME )
  message( FATAL_ERROR "Require TEST_NAME to be defined." )
endif( NOT TEST_NAME )
if( NOT TEST_PROGRAM )
  message( FATAL_ERROR "Require TEST_PROGRAM to be defined." )
endif( NOT TEST_PROGRAM )
if( NOT TEST_OUTPUT )
  message( FATAL_ERROR "Require TEST_OUTPUT to be defined" )
endif( NOT TEST_OUTPUT )
if( NOT TEST_REFERENCE )
  message( FATAL_ERROR "Require TEST_REFERENCE to be defined" )
endif( NOT TEST_REFERENCE )
if( NOT TEST_DIFF )
message( FATAL_ERROR "Require TEST_DIFF to be defined" )
endif( NOT TEST_DIFF )

# create a directory for the test
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/gwb-dat/${TEST_NAME})

if( NOT ${TEST_DAT} EQUAL "" )
  FILE(STRINGS ${TEST_DAT} _input_lines REGEX "EXPECT FAILURE")
endif()
if("${_input_lines}" STREQUAL "")
  SET(EXPECT 0)
else()
  SET(EXPECT 1)
  MESSAGE(STATUS "Test ${testname} is expected to fail.")
endif()

set(EXECUTE_COMMAND ${TEST_PROGRAM} ${TEST_ARGS})

# run the test program, capture the stdout/stderr and the result var ${TEST_ARGS}
execute_process(
  COMMAND ${TEST_PROGRAM} ${TEST_ARGS} 
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/gwb-dat/ 
  OUTPUT_FILE ${TEST_OUTPUT}
  ERROR_VARIABLE TEST_ERROR_VAR
  RESULT_VARIABLE TEST_RESULT_VAR
  OUTPUT_VARIABLE TEST_OUTPUT_VAR
  )

# if the return value is !=0 bail out
if( TEST_ERROR_VAR EQUAL EXPECT )
	message( FATAL_ERROR "Failed: Test program ${TEST_PROGRAM} exited != ${EXPECT}.\n${TEST_ERROR_VAR}" )
endif()

# process test output
file(READ ${TEST_OUTPUT} TMP_TEST_OUTPUT)
String(REGEX REPLACE "[:/a-zA-Z\\/0-9_\\.-]*gwb-grid[.exe]*" "(..path..)bin/gwb-grid" TEST_OUTPUT_PROCESSED "${TMP_TEST_OUTPUT}")
String(REGEX REPLACE "  -j N                  Specify the number of threads the visualizer is allowed to use. Default: [0-9]*." "  -j N                  Specify the number of threads the visualizer is allowed to use. Default: --." TEST_OUTPUT_PROCESSED "${TEST_OUTPUT_PROCESSED}")
String(REGEX REPLACE "GWB Version: [0-9.a-zA-Z]*" "GWB Version: -.-.-.-" TEST_OUTPUT_PROCESSED "${TEST_OUTPUT_PROCESSED}")
String(REGEX REPLACE "git hash: [0-9.a-zA-Z-]* branch: [-0-9.a-zA-Z_]*" "git hash: (..git hash..) branch: (..git branch..)" TEST_OUTPUT_PROCESSED "${TEST_OUTPUT_PROCESSED}")
String(REGEX REPLACE "Starting the world builder with ([0-9]*) threads..." "Starting the world builder with -- threads..." TEST_OUTPUT_PROCESSED "${TEST_OUTPUT_PROCESSED}")
file(WRITE ${TEST_OUTPUT} ${TEST_OUTPUT_PROCESSED})

if( TEST_RESULT_VAR )
  String(REGEX REPLACE "(.*)AssertThrow `(.*)` failed in (.*)/source/([^ )(]*) at line ([01234567890]*): (.*)" "AssertThrow `(..error type..)` failed in (..path..)/source/\\4 at line (..line..): \\6" TEST_ERROR_VAR_PROCESSED "${TEST_ERROR_VAR}")
  string(FIND "${TEST_ERROR_VAR_PROCESSED}" "AssertThrow" FIRST_BRACKET)
  string(FIND "${TEST_ERROR_VAR_PROCESSED}" "Error not recoverable, aborting program." LAST_BRACKET REVERSE)
  MATH(EXPR LAST_BRACKET ${LAST_BRACKET})
  string(SUBSTRING "${TEST_ERROR_VAR_PROCESSED}" "${FIRST_BRACKET}" "${LAST_BRACKET}" TEST_ERROR_VAR_PROCESSED)
  string(REGEX REPLACE "Could not open file <[:/a-zA-Z\\/0-9_\\.-]*>" "Could not open file <..file..>" TEST_ERROR_VAR_PROCESSED "${TEST_ERROR_VAR_PROCESSED}")
  file(APPEND ${TEST_OUTPUT} "Expected fail with: \n${TEST_ERROR_VAR_PROCESSED}")
endif()

file(TO_NATIVE_PATH "${TEST_OUTPUT}" TEST_NATIVE_OUTPUT)
file(TO_NATIVE_PATH "${TEST_REFERENCE}" TEST_NATIVE_REFERENCE)


IF("${TEST_DIFF}" MATCHES ".*exe")
  # windows
  FIND_PROGRAM(DOS2UNIX_EXECUTABLE
	     NAMES dos2unix
	     HINTS ${DIFF_DIR}
	     PATH_SUFFIXES bin
	     )
     IF(NOT DOS2UNIX_EXECUTABLE MATCHES "-NOTFOUND")
	     SET(TEST_D2U ${DOS2UNIX_EXECUTABLE})
     ELSE()
	     MESSAGE(FATAL_ERROR
		     "Could not find dos2unix. This is required for running the testsuite in windows.\n"
		     "Please specify TEST_D2U by hand."
		     )
     ENDIF()
     execute_process(COMMAND ${TEST_D2U} ${TEST_NATIVE_OUTPUT})
     execute_process(COMMAND ${TEST_D2U} ${TEST_NATIVE_REFERENCE})
ENDIF()

# now compare the output with the reference
execute_process(
	COMMAND ${TEST_DIFF} -q ${TEST_NATIVE_OUTPUT} ${TEST_NATIVE_REFERENCE}
  RESULT_VARIABLE TEST_RESULT
  )

  file(READ ${TEST_NATIVE_OUTPUT} TEST_NATIVE_OUTPUT_CONTENT)

# again, if return value is !=0 scream and shout
if( TEST_RESULT )
	execute_process(COMMAND ${TEST_DIFF} ${TEST_NATIVE_OUTPUT} ${TEST_NATIVE_REFERENCE})
	message( FATAL_ERROR "Failed: The output of ${TEST_NAME} stored in ${TEST_NATIVE_OUTPUT} did not match the reference output stored in ${TEST_NATIVE_REFERENCE}\n\n Full output was:\n ${TEST_NATIVE_OUTPUT_CONTENT}")
endif( TEST_RESULT )
