# A script for the stream editor sed to update .cc and .h files from the
# naming scheme used in ASPECT 1.5.0 to ASPECT 2.0.0.

# This script is now outdated, because it interferes with valid names in the current
# ASPECT version. It is kept in case it is needed, and you can reactivate the
# script by moving it into the parent folder. Use with care and check
# every change that this script suggests.

# Rename fluid pressure boundary conditions
s/fluid_pressure_boundary_conditions/boundary_fluid_pressure/g
s/FluidPressureBoundaryConditions/BoundaryFluidPressure/g
s/fluid_pressure_boundary/boundary_fluid_pressure/g
s/Fluid pressure boundary/Boundary fluid pressure/g
s/FLUID_PRESSURE_BOUNDARY_CONDITIONS/BOUNDARY_FLUID_PRESSURE_MODEL/g
s/\b\@ingroup BoundaryFluidPressure\b/\@ingroup BoundaryFluidPressures/g

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

# Rename assembler base class
s/internal::Assembly::AssemblerLists<dim>/Assemblers::AssemblerLists<dim>/g
s/internal::Assembly::Assemblers/Assemblers/g
s/struct AssemblerLists/class Manager/g
s/AssemblerLists/Manager/g
s/AssemblerBase/Interface/g
s:assembly.h:simulator/assemblers/interface.h:g

# Rename adiabatic conditions plugin initial profile includes
s:#include <aspect/adiabatic_conditions/initial_profile.h>:#include <aspect/adiabatic_conditions/compute_profile.h>:g

# Rename nonlinear solver schemes
s/NonlinearSolver::iterated_IMPES/NonlinearSolver::iterated_Advection_and_Stokes/g
s/NonlinearSolver::IMPES/NonlinearSolver::single_Advection_single_Stokes/g
s/NonlinearSolver::iterated Stokes/NonlinearSolver::single_Advection_iterated_Stokes/g
s/NonlinearSolver::Stokes_only/NonlinearSolver::no_Advection_iterated_Stokes/g
s/NonlinearSolver::Advection_only/NonlinearSolver::single_Advection_no_Stokes/g
s/NonlinearSolver::Newton_Stokes/NonlinearSolver::iterated_Advection_and_Newton_Stokes/g

# Account for moved parameters
s/this->get_parameters().fixed_composition_boundary_indicators/this->get_boundary_composition_manager().get_fixed_composition_boundary_indicators()/g
s/this->get_parameters().fixed_temperature_boundary_indicators/this->get_boundary_temperature_manager().get_fixed_temperature_boundary_indicators()/g
s/parameters.fixed_composition_boundary_indicators/this->get_fixed_composition_boundary_indicators()/g
s/parameters.fixed_temperature_boundary_indicators/this->get_fixed_temperature_boundary_indicators()/g
