# A script for the stream editor sed to update .prm files from the
# naming scheme used in ASPECT 1.5.0 to ASPECT 2.0.0.
# This part deals with the old options in the 'Model settings' subsection.
# They have been moved into more appropriate subsections, but sed is
# not good at moving things above their current position. Instead
# we store everything correctly formatted in the hold space and
# replace the Model settings subsection with the new subsections and
# parameters. Since parameter files do not mind having the same subsection
# multiple times, this works just fine.

:jump
/subsection Model settings/,/^[[:space:]]*end[[:space:]]*/ {
# Remove empty lines, they would make the rest complicated.
/^ *$/d

/^ *#/{
# If there are comments add the next line
N

# If a comment was not followed by a 'set' we might have
# picked up the 'end' that was supposed to stay there.
# This is a problem, because we can not notice that the subsection
# 'Model settings' is over if we do not find the 'end'.
# Thus, replace the content of the pattern by a single 'end' and
# jump before the place, where we look for the 'end'.
/\n *end/{
s/^.*$/end/
b jump
}

# Repeat the additions until we add
# a line that does not start with a '#'.
# The ^\n pattern makes sure we do not search in already added lines.
/\n *#[^\n][^\n]*$/ b jump
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
s/^.*$/subsection Nullspace removal\
&\
end\
/
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
  set List of model names = PARAMETER NEEDS MANUAL CHECK:adiabatic heating\
end\
/
H
d
}

/set Include shear heating *= *true/ {
s/^.*$/# Warning: Merge the parameter in the following subsection with any\
# other parameter of the same name, otherwise your model results might change\
subsection Formulation\
  set List of model names = PARAMETER NEEDS MANUAL CHECK:shear heating\
end\
/
H
d
}


/set Include latent heat *= *true/ {
s/^.*$/# Warning: Merge the parameter in the following subsection with any\
# other parameter of the same name, otherwise your model results might change\
subsection Formulation\
  set List of model names = PARAMETER NEEDS MANUAL CHECK:latent heat\
end\
/
H
d
}



# Now the Model settings subsection should be empty.
# Remove any remaining lines (e.g. empty lines, comments),
# except for the last one (to not jump over it).
/^[[:space:]]*end/ !d

/^[[:space:]]*end/ {
# At the end of the subsection, add a comment and the content of the hold space
x

# if the hold space is not empty, print comment
/..*/ {
i\
# The parameters below this comment were created by the update script\
# as replacement for the old 'Model settings' subsection. They can be\
# safely merged with any existing subsections with the same name.

# remove the empty line at the end of the hold space
s/\n *$//
}

# If the hold space was empty, remove the single newline it would produce
/^$/ d
}
}
