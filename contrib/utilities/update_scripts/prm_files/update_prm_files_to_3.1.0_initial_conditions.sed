# A script for the stream editor sed to update .prm files from
# using Model name in the Initial temperature and Initial composition
# subsections to using List of model names in ASPECT 3.1.0-pre.

# Replace the 'model name' parameter by 'List of model names'.
# Note that this command only works if the parameter is set
# before the next 'end', which is not necessarily the one that
# belongs to the opening subsection (i.e. if the parameter is set
# after a subsection nested inside the 'Initial temperature model'
# subsection the following will simply do nothing).
/subsection Initial temperature model/,/^[[:space:]]*end[[:space:]]*/ {
     s/set Model name/set List of model names/g
}

/subsection Initial composition model/,/^[[:space:]]*end[[:space:]]*/ {
     s/set Model name/set List of model names/g
}
