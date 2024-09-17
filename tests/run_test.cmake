execute_process(COMMAND ${CMAKE_COMMAND}
    --build ${BINARY_DIR} --target ${TESTNAME}
    RESULT_VARIABLE _result_code
    OUTPUT_VARIABLE _output
    )

if("${_result_code}" STREQUAL "0")
    message("${TESTNAME}: success.")

else()

    message("*** ${ERROR}: ***")
    message(${_output})
    message(FATAL_ERROR "*** Test aborted.")
endif()
