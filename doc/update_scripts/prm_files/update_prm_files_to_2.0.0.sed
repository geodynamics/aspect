# A script for the stream editor sed to update .prm files from the
# naming scheme used in ASPECT 1.5.0 to ASPECT 2.0.0.

# Rename compositional initial conditions
s/subsection Compositional initial conditions/subsection Initial composition model/g

# Rename initial (temperature) conditions
s/subsection Initial conditions/subsection Initial temperature model/g

# Rename adiabatic conditions plugin initial profile
s/subsection Initial profile/subsection Compute profile/g
s/set Model name = initial profile/set Model name = compute profile/g

# Recover previous behavior of the 'Radial linear' gravity model
# In the following part we read the whole subsection 'Radial linear'
# when we find it, check if 'Magnitude at bottom' is inside, and
# insert it, if not.
:jump_before_gravity
/subsection Radial linear/,/[[:space:]]*end[[:space:]]*/ {
# Print the line 'subsection ...' and read next line
n
# If next line was already 'end' jump out of this part
/^[[:space:]]*end/b jump_before_gravity

:jump_in_gravity
# Add next line to current pattern space until we found 'end'
N
/end/!b jump_in_gravity

# If we get here, we found an end. If we did not find the new parameter,
# it needs to be added.
/set Magnitude at bottom/!{
i\
    set Magnitude at bottom = 0.0
}

# Then jump up to finish this part, and print the whole modified subsection
b jump_before_gravity
}

# Replace the 'model name' parameter by 'List of model names'
# in all subsections that now use the new parameter.
# Note that this command only works if the parameter is set
# before the next 'end', which is not necessarily the one that
# belongs to the opening subsection (i.e. if the parameter is set
# after a subsection nested inside the 'Boundary temperature model'
# subsection the following will simply do nothing).
/subsection Boundary temperature model/,/^[[:space:]]*end[[:space:]]*/ {
     s/set Model name/set List of model names/g
}

/subsection Boundary composition model/,/^[[:space:]]*end[[:space:]]*/ {
     s/set Model name/set List of model names/g
}

/subsection Heating model/,/^[[:space:]]*end[[:space:]]*/ {
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

# Rename solver schemes:
s/iterated IMPES/iterated Advection and Stokes/g
s/IMPES/single Advection, single Stokes/g
s/= *iterated Stokes/= single Advection, iterated Stokes/g
s/Stokes only/no Advection, iterated Stokes/g
s/Advection only/single Advection, no Stokes/g
s/= *Newton Stokes/= iterated Advection and Newton Stokes/g
