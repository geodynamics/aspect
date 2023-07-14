cmake_minimum_required (VERSION 2.6)
project(WorldBuilderFortranExample Fortran)
include_directories(${CMAKE_BINARY_DIR}/../../build/mod/ ${CMAKE_BINARY_DIR}/../../build/lib/)
add_executable(WorldBuilderFortranExample "${CMAKE_CURRENT_SOURCE_DIR}/example.f90")

# Make sure that the whole library is loaded, so the registration is done correctly.
if(NOT APPLE)
    SET(GWB_LIBRARY_WHOLE -L../../build/lib/ -Wl,--whole-archive -lWorldBuilder -Wl,--no-whole-archive -lstdc++ -I../../build/mod/)  
else()
    SET(GWB_LIBRARY_WHOLE -L../../build/lib/ -Wl,-force_load,../../lib/libWorldBuilder.a -lc++ -I../../build/mod/)
endif()

target_link_libraries(WorldBuilderFortranExample ${GWB_LIBRARY_WHOLE})