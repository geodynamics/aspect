#!/bin/sh

# This is the default test processing script, which is used to pre-process the
# screen output unless you create a <testname>.sh next to your test. The
# script gets one argument ($1), which is the name of the file (currently this
# will always be "screen-output") and the input is piped into the script.
# This implementation here just prints the input unmodified.

cat
