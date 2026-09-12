CMAKE_MINIMUM_REQUIRED(VERSION 3.13.4)

FIND_PACKAGE(Aspect REQUIRED HINTS ${Aspect_DIR} ../ ../../ $ENV{ASPECT_DIR})
DEAL_II_INITIALIZE_CACHED_VARIABLES()

SET(TARGET "2D_subduction_with_two_phase_flow")
PROJECT(${TARGET})

# ADD_LIBRARY(${TARGET} SHARED world_builder_too_high.cc fixed_reactions.cc )
ADD_LIBRARY(${TARGET} SHARED bound_fluid_source.cc world_builder_disequilibrium_fluid.cc)
ASPECT_SETUP_PLUGIN(${TARGET})
