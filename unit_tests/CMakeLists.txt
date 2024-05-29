################### top matter ################
CMAKE_MINIMUM_REQUIRED(VERSION 3.13.4)

PROJECT(unittests CXX)

SET(ASPECT_BINARY "${ASPECT_BINARY}" CACHE STRING "" FORCE)

ENABLE_TESTING()

MESSAGE(STATUS "creating unit_tests ctest in ${CMAKE_CURRENT_BINARY_DIR} binary=${ASPECT_BINARY}")

ADD_TEST(NAME unit_tests
         COMMAND ${ASPECT_BINARY} --test
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      )
