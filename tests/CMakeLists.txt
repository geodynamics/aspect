
FIND_PACKAGE(Perl)

ADD_CUSTOM_TARGET(tests)

FILE(GLOB _tests *.prm)
FOREACH(_test ${_tests})

  GET_FILENAME_COMPONENT(_test ${_test} NAME_WE) 
  IF(EXISTS ${_test}.cc)
    ADD_LIBRARY(${_test} SHARED EXCLUDE_FROM_ALL ${_test}.cc)
    DEAL_II_SETUP_TARGET(${_test})
    SET(_testdepends ${_test})
  ELSE()
    SET(_testdepends)
  ENDIF()

  ADD_CUSTOM_COMMAND(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/output-${_test}/screen-output
    COMMAND if test ! -d ${CMAKE_CURRENT_BINARY_DIR}/output-${_test} \; then
            mkdir ${CMAKE_CURRENT_BINARY_DIR}/output-${_test} \;
            fi
    COMMAND aspect ${CMAKE_CURRENT_SOURCE_DIR}/${_test}.prm 
            > ${CMAKE_CURRENT_BINARY_DIR}/output-${_test}/screen-output
    COMMAND ${PERL_EXECUTABLE} -pi
            ${CMAKE_CURRENT_SOURCE_DIR}/normalize.pl
            ${CMAKE_CURRENT_BINARY_DIR}/output-${_test}/*
    DEPENDS aspect ${CMAKE_CURRENT_SOURCE_DIR}/${_test}.prm ${_testdepends}
    )

  ADD_CUSTOM_TARGET(${_test}.run 
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/output-${_test}/screen-output
    )
# numdiff numdiff -a 1e-6 -q -s ' \t\n:'
  ADD_CUSTOM_COMMAND(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/diff-${_test}
    COMMAND
      DEAL_II_DIFF=\"diff\" \;
      for i in ${CMAKE_CURRENT_SOURCE_DIR}/${_test}/* \; do
        basename=`echo $$i | perl -pi -e 's!.*/!!g\;'` \;
	cat ${CMAKE_CURRENT_BINARY_DIR}/output-${_test}/$$basename | egrep -v '^\\|' > ${CMAKE_CURRENT_BINARY_DIR}/output-${_test}/$$basename.notime \;
	$${DEAL_II_DIFF} ${CMAKE_CURRENT_SOURCE_DIR}/${_test}/$$basename ${CMAKE_CURRENT_BINARY_DIR}/output-${_test}/$$basename.notime \;
	if $$_ \; then 
          :
        else 
          echo "Files ${CMAKE_CURRENT_SOURCE_DIR}/${_test}/$$basename ${CMAKE_CURRENT_BINARY_DIR}/output-${_test}/$$basename.notime differ" \; exit 1 \;
        fi \;
      done
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/output-${_test}/screen-output
    )

  ADD_CUSTOM_TARGET(${_test}.diff 
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/diff-${_test}
    )

  ADD_DEPENDENCIES(tests ${_test}.diff)

ENDFOREACH()

