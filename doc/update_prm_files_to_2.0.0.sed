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
s/subsection Compositional initial conditions/Initial composition model/g

# Rename initial (temperature) conditions
s/subsection Initial conditions/Initial temperature model/g

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
/subsection Boundary temperature model/,/^ *\bend\b/ {
     s/set Model name/set List of model names/g
}

/subsection Boundary composition model/,/^ *\bend\b/ {
     s/set Model name/set List of model names/g
}

/subsection Heating model/,/^ *\bend\b/ {
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



# This part deals with the old options in the 'Model settings' subsection.
# They have been moved into more appropriate subsections, but sed is
# not good at moving things above their current position. Instead
# we store everything correctly formatted in the hold space and
# add it to the end of the file. Since parameter files do not mind
# having the same subsection multiple times, this should work
# just fine.
:jump
/subsection Model settings/,/^ *\bend\b/ {
# Remove empty lines, they would make the rest complicated.
/^ *$/d

/^ *#/{
# If there are comments add the next line
N

# If a comment was not followed by a 'set' we might have
# picked up the 'end' that was supposed to stay there.
# This is a problem, because we can not notice that the subsection
# 'Model settings' is over if we do not find the 'end'.
# Thus, replace the content of the pattern by a singe 'end' and
# jump before the place, where we look for the 'end'.
/\n *end/{
s/^.*$/end/
b jump
}

# Repeat the additions until we add
# a line that does not start with a '#'.
# The ^\n pattern makes sure we do not search in already added lines.
/ *#[^\n][^\n]*$/ b jump
}

# Add more lines to the current pattern, 
# as long as they end with a '\'.
/\\ *$/{
N
/\\ *$/ b jump
}

# Now we have all related comments and the current option in the 
# pattern space. Remove all options, which are set to default
# values, to reduce
# the number of lines we need to move in the next step.
# If there are inline comments, keep the lines.
/set Fixed temperature boundary indicators *= *$/d
/set Fixed composition boundary indicators *= *$/d
/set Tangential velocity boundary indicators *= *$/d
/set Zero velocity boundary indicators *= *$/d
/set Prescribed velocity boundary indicators *= *$/d
/set Prescribed traction boundary indicators *= *$/d
/set Include melt transport *= *false *$/d
/set Free surface boundary indicators *= *$/d
/set Remove nullspace *= *$/d
/set Enable additional Stokes RHS *= *false *$/d
/set Include adiabatic heating *= *false/d
/set Include shear heating *= *false/d
/set Include latent heat *= *false/d

# Now append all options and the appropriate subsection
# in the hold space one by one, and delete them from the stream afterwards.

/set Fixed temperature boundary indicators/ {
s/^.*$/subsection Boundary temperature model\
&\
end\
/
H
d
}

/set Fixed composition boundary indicators/ {
s/^.*$/subsection Boundary composition model\
&\
end\
/
H
d
}

/set Tangential velocity boundary indicators/ {
s/^.*$/subsection Boundary velocity model\
&\
end\
/
H
d
}

/set Zero velocity boundary indicators/ {
s/^.*$/subsection Boundary velocity model\
&\
end\
/
H
d
}

/set Prescribed velocity boundary indicators/ {
s/^.*$/subsection Boundary velocity model\
&\
end\
/
H
d
}

/set Prescribed traction boundary indicators/ {
s/^.*$/subsection Boundary traction model\
&\
end\
/
H
d
}

/set Include melt transport/ {
s/^.*$/subsection Melt settings\
&\
end\
/
H
d
}

/set Free surface boundary indicators/ {
s/^.*$/subsection Free surface\
&\
end\
/
H
d
}

/set Remove nullspace/ {
s/^.*$/&\
/
s/ *set Remove nullspace/set Remove nullspace/
H
d
}

/set Enable additional Stokes RHS/ {
s/^.*$/subsection Formulation\
&\
end\
/
H
d
}

/set Include adiabatic heating *= *true/ {
s/^.*$/# Warning: Merge the parameter in the following subsection with any\
# other parameter of the same name, otherwise your model results might change\
subsection Heating model\
  set List of model names = adiabatic heating\
end\
/
H
d
}

/set Include shear heating *= *true/ {
s/^.*$/# Warning: Merge the parameter in the following subsection with any\
# other parameter of the same name, otherwise your model results might change\
subsection Formulation\
  set List of model names = shear heating\
end\
/
H
d
}


/set Include latent heat *= *true/ {
s/^.*$/# Warning: Merge the parameter in the following subsection with any\
# other parameter of the same name, otherwise your model results might change\
subsection Formulation\
  set List of model names = latent heat\
end\
/
H
d
}



# Now the Model settings subsection should be empty.
# Remove any remaining lines (e.g. empty lines, comments),
# except for the last one (to not jump over it).
/^ *end/ !d

/^ *end/ {
# At the end of the file, add a comment and the content of the hold space
x

# if the hold space is not empty, print comment
/..*/ {
i # The parameters below this comment were created by the update script\
# as replacement for the old 'Model settings' subsection. They can be\
# safely merged with any existing subsections with the same name.

# remove the empty line at the end of the hold space
s/\n *$//
}

# If the hold space was empty, remove the single newline it would produce
/^$/ d
}
}
