# A script for the stream editor sed to update .prm files from the
# naming scheme used in ASPECT 1.5.0 to ASPECT 2.0.0. This script correctly
# updated all files within the official development version of ASPECT,
# but it is not guaranteed to work for all possible names in user files.
# Consequently a backup of files is strongly recommended, and the changes
# created by this script should be investigated to ensure a correct renaming.
#
# Usage for a parameter file named FILENAME (possibly containing
# wildcards such as '*.prm'):
# sed -i -f update_prm_files_to_2.0.0.sed FILENAME
# On MacOS:
# sed -i "" -f update_source_files_to_2.0.0.sed FILENAME

# Rename compositional initial conditions
s/subsection Compositional initial conditions/subsection Initial composition model/g

# Rename initial (temperature) conditions
s/subsection Initial conditions/subsection Initial temperature model/g

# Rename adiabatic conditions plugin initial profile
s/subsection Initial profile/subsection Compute profile/g
s/set Model name = initial profile/set Model name = compute profile/g

# Replace the 'model name' parameter by 'List of model names'
# in all subsections that now use the new parameter.
# Note that this command only works if the parameter is set
# before the next 'end', which is not necessarily the one that
# belongs to the opening subsection (i.e. if the parameter is set
# after a subsection nested inside the 'Boundary temperature model'
# subsection the following will simply do nothing).
/subsection Boundary temperature model/,/\bend\b/ {
     s/set Model name/set List of model names/g
}

/subsection Boundary composition model/,/\bend\b/ {
     s/set Model name/set List of model names/g
}

# Rename tracers to particles
s/tracer/particle/g
s/Tracer/Particle/g

# Rename velocity boundary conditions
s/velocity-boundary-conditions/boundary-velocity/g

# Rename compositional initial conditions
s/compositional-initial-conditions/initial-composition/g

# Rename initial (temperature) conditions
s/initial-conditions/initial-temperature/g

# Replace removed postprocessors
s/seismic vs/named additional outputs/g
s/seismic vp/named additional outputs/g
s/named additional outputs.*named additional outputs/named additional outputs/g
s/friction heating/heating/g
s/viscous dissipation statistics/heating statistics/g
s/heating statistics.*heating statistics/heating statistics/g


# Make all instances of boussinesq into Boussinesq
s/boussinesq/Boussinesq/g

# Remove all instances of `Subtract mean of dynamic topography`
/Subtract mean of dynamic topography/d

# Rename Dynamic topography subsection
s/subsection Dynamic Topography/subsection Dynamic topography/g

# Remove all instances of `melt transport threshold`
/Melt transport threshold/d
