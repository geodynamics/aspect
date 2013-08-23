EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND}
    --build ${BINARY_DIR} --target ${TESTNAME}
    RESULT_VARIABLE _result_code
    OUTPUT_VARIABLE _output
    )

IF("${_result_code}" STREQUAL "0")
    MESSAGE("${TESTNAME}: success.")

ELSE()

    MESSAGE("*** ${ERROR}: ***")
    MESSAGE(${_output})
    MESSAGE(FATAL_ERROR "*** Test aborted.")
ENDIF()
