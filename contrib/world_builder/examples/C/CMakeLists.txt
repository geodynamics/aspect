cmake_minimum_required (VERSION 2.6)
project(WorldBuilderCExample C)
include_directories("${CMAKE_BINARY_DIR}/../../include/" "${CMAKE_BINARY_DIR}/../../build/lib/")
add_executable(WorldBuilderCExample "${CMAKE_CURRENT_SOURCE_DIR}/example.c")

# Make sure that the whole library is loaded, so the registration is done correctly.
if(NOT APPLE)
    SET(GWB_LIBRARY_WHOLE -L../../build/lib/ -Wl,--whole-archive -lWorldBuilder -Wl,--no-whole-archive -lstdc++ -lm )  
else()
    SET(GWB_LIBRARY_WHOLE -L../../build/lib/ -Wl,-force_load,../../lib/libWorldBuilder.a -lc++ -lm)
endif()

target_link_libraries(WorldBuilderCExample ${GWB_LIBRARY_WHOLE})
