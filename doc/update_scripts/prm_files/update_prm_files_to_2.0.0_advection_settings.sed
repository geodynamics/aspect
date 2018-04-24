# A script for the stream editor sed to update .prm files from the
# naming scheme used in ASPECT 1.5.0 to ASPECT 2.0.0. 
# This file particularly moves advection solver parameters into the new
# 'Solver parameters' subsection.

# Remove lines with default values.
/ *set Composition solver tolerance *= *1e-12/d
/ *set Temperature solver tolerance *= *1e-12/d

:solver_settings

# Look for solver parameters to be moved to the
# 'Solver parameters' subsection. If found,
# store them in hold space (H), and delete them from their
# current place (d), then branch to :solver_settings to check
# next line.
/subsection Solver parameters/,/end/! {
/Temperature solver tolerance/ {
s/set Temperature solver tolerance *=/  set Temperature solver tolerance =/
H
# Delete line and read next, except at end of file
# At end of file, just replace by empty line
$ !d
$ s/.*//
b solver_settings
}

/Composition solver tolerance/ {
s/set Composition solver tolerance *=/  set Composition solver tolerance =/
H
# Delete line and read next, except at end of file
# At end of file, just replace by empty line
$ !d
$ s/.*//
b solver_settings
}
}

# If we find an existing 'Solver parameter' subsection
# insert any found parameters in here
/subsection Solver parameters/ {
x
# If we have stored parameters, print the header and the parameters
/..*/ {
x
# Print the last line if it is not empty
p
s/.*//
x
s/^ *\n//
a\

}

# If the hold space was empty, just switch back to the last line
/^$/ {
x
}
}

# If we have read the whole file, print all the accumulated parameters
# into a new 'Solver parameters subsection' 
$ {
# If the hold space is not empty, print new subsection
x
/..*/ {
x
# Print the last line if it is not empty
/^ *$/! p
x

# Remove empty lines
s/^$//
s/^\n//
s/\n$//

# Add surrounding subsection and print
s/.*/\
subsection Solver parameters\
&\
end/
p
# Delete the content of hold space (it was printed before the new subsection)
x
d
}

# If the hold space was empty, just switch back to the last line
/^$/ {
x
}
}
