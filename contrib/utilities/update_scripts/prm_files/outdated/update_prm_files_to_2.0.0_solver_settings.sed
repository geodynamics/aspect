# A script for the stream editor sed to update .prm files from the
# naming scheme used in ASPECT 1.5.0 to ASPECT 2.0.0. 
# This file particularly moves old parameters into the new
# 'Solver parameters/Stokes solver parameters' subsection.

:solver_settings

# Look for solver parameters to be moved to the new 
# 'Stokes solver parameters' subsection. If found,
# store them in hold space (H), and delete them from their
# current place (d), then branch to :solver_settings to check
# next line.
/subsection Stokes solver parameters/,/end/! {
/Use direct solver for Stokes system/ {
s/set Use direct solver for Stokes system *=/    set Use direct solver for Stokes system =/
H
# Delete line and read next, except at end of file
# At end of file, just replace by empty line
$ !d
$ s/.*//
b solver_settings
}

/Linear solver tolerance/ {
s/set Linear solver tolerance *=/    set Linear solver tolerance =/
H
# Delete line and read next, except at end of file
# At end of file, just replace by empty line
$ !d
$ s/.*//
b solver_settings
}

/Linear solver A block tolerance/ {
s/set Linear solver A block tolerance *=/    set Linear solver A block tolerance =/
H
# Delete line and read next, except at end of file
# At end of file, just replace by empty line
$ !d
$ s/.*//
b solver_settings
}

/Linear solver S block tolerance/ {
s/set Linear solver S block tolerance *=/    set Linear solver S block tolerance =/
H
# Delete line and read next, except at end of file
# At end of file, just replace by empty line
$ !d
$ s/.*//
b solver_settings
}

/Number of cheap Stokes solver steps/ {
s/set Number of cheap Stokes solver steps *=/    set Number of cheap Stokes solver steps =/
H
# Delete line and read next, except at end of file
# At end of file, just replace by empty line
$ !d
$ s/.*//
b solver_settings
}

/Maximum number of expensive Stokes solver steps/ {
s/set Maximum number of expensive Stokes solver steps *=/    set Maximum number of expensive Stokes solver steps =/
H
# Delete line and read next, except at end of file
# At end of file, just replace by empty line
$ !d
$ s/.*//
b solver_settings
}
}

# If we have read the whole file, print all the accumulated parameters
# into the new 'Stokes solver parameters subsection' 
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
  subsection Stokes solver parameters\
&\
  end\
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
