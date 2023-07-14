# arguments checking
if( NOT TEST_NAME )
  message( FATAL_ERROR "Require TEST_NAME to be defined." )
endif( NOT TEST_NAME )
if( NOT TEST_PROGRAM )
  message( FATAL_ERROR "Require TEST_PROGRAM to be defined." )
endif( NOT TEST_PROGRAM )

# create a directory for the test
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/documentation/${TEST_NAME})

set(EXECUTE_COMMAND ${TEST_PROGRAM} ${TEST_ARGS})

# run the test program, capture the stdout/stderr and the result var ${TEST_ARGS}
execute_process(
  COMMAND ${TEST_PROGRAM} ${TEST_ARGS} 
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/documentation/ 
  OUTPUT_FILE ${TEST_OUTPUT}.log
  ERROR_VARIABLE TEST_ERROR_VAR
  RESULT_VARIABLE TEST_RESULT_VAR
  OUTPUT_VARIABLE TEST_OUTPUT_VAR
  )

# if the return value is !=0 bail out
if( TEST_RESULT_VAR )
	message( FATAL_ERROR "Failed: Test program ${TEST_PROGRAM} exited != 0.\n${TEST_ERROR_VAR}" )
endif( TEST_RESULT_VAR )
