# A script for the stream editor sed to update .cc and .h files from the
# naming scheme used in ASPECT 1.5.0 to ASPECT 2.0.0. This script correctly
# updated all files within the official development version of ASPECT,
# but it is not guaranteed to work for all possible names in user plugins.
# Consequently a backup of files is strongly recommended, and the changes
# created by this script should be investigated to ensure a correct renaming.
#
# Usage for a c++ source or header file named FILENAME (possibly containing
# wildcards such as '*.cc') on Linux:
# sed -i -f update_source_files_to_2.0.0.sed FILENAME
# On MacOS:
# sed -i "" -f update_source_files_to_2.0.0.sed FILENAME

# Rename fluid pressure boundary conditions
s/fluid_pressure_boundary_conditions/boundary_fluid_pressure/g
s/FluidPressureBoundaryConditions/BoundaryFluidPressure/g
s/fluid_pressure_boundary/boundary_fluid_pressure/g
s/Fluid pressure boundary/Boundary fluid pressure/g
s/FLUID_PRESSURE_BOUNDARY_CONDITIONS/BOUNDARY_FLUID_PRESSURE_MODEL/g
s/\<\@ingroup BoundaryFluidPressure\>/\@ingroup BoundaryFluidPressures/g

# Rename traction boundary conditions
s/traction_boundary_conditions_model/boundary_traction/g
s/traction_boundary_conditions/boundary_traction/g
s/(?<!prescribed_)traction_boundary/boundary_traction/g
s/TractionBoundaryConditions/BoundaryTraction/g
s/Traction boundary/Boundary traction/g
s/TRACTION_BOUNDARY_CONDITIONS/BOUNDARY_TRACTION_MODEL/g
s/\@ingroup BoundaryTractionModels/\@ingroup BoundaryTractions/g

# Rename velocity boundary conditions
s/velocity_boundary_conditions_model/boundary_velocity/g
s/velocity_boundary_conditions/boundary_velocity/g
s/velocity-boundary-conditions/boundary-velocity/g
s/VelocityBoundaryConditions/BoundaryVelocity/g
s/Velocity boundary conditions/Boundary velocity/g
s/VELOCITY_BOUNDARY_CONDITIONS/BOUNDARY_VELOCITY_MODEL/g
s/\@ingroup BoundaryVelocityModels/\@ingroup BoundaryVelocities/g

# Rename compositional initial conditions
s/compositional_initial_conditions__model/initial_composition/g
s/compositional_initial_conditions/initial_composition/g
s/compositional-initial-conditions/initial-composition/g
s/CompositionalInitialConditions/InitialComposition/g
s/Compositional initial conditions/Initial composition model/g
s/COMPOSITIONAL_INITIAL_CONDITIONS/INITIAL_COMPOSITION_MODEL/g
s/\@ingroup CompositionalInitialConditionsModels/\@ingroup InitialCompositions/g

# Rename initial (temperature) conditions
s/initial_conditions_model/initial_temperature/g
s/initial_conditions\([^\.][^c]\)/initial_temperature\1/g
s/initial-conditions\//initial-temperature\//g
s/initial_condition[^s]/initial_temperature/g
s/InitialConditionsModels/InitialTemperatures/g
s/InitialConditions/InitialTemperature/g
s/enter_subsection*("Initial conditions")/enter_subsection ("Initial temperature model")/g
s/INITIAL_CONDITIONS/INITIAL_TEMPERATURE_MODEL/g
s/\@ingroup InitialConditionsModels/\@ingroup InitialTemperatures/g

# Rename tracers to particles
s/tracer particle/particle/g
s/tracer/particle/g
s/Tracer/Particle/g
s/TRACER/PARTICLE/g
