################### top matter ################
cmake_minimum_required(VERSION 3.13.4)

project(unittests CXX)

set(ASPECT_BINARY "${ASPECT_BINARY}" CACHE STRING "" FORCE)

enable_testing()

message(STATUS "creating unit_tests ctest in ${CMAKE_CURRENT_BINARY_DIR} binary=${ASPECT_BINARY}")

add_test(NAME unit_tests
         COMMAND ${ASPECT_BINARY} --test
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      )
