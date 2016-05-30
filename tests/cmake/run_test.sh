#!/bin/bash

# This script is executed to run ASPECT with a given .prm test and process the
# return value of ASPECT correctly.

# usage: EXPECT ./aspect test.prm output-name
# where EXPECT is the expected return value

EXPECT="$1"
BINARY="$2"
PRM="$3"
OUTPUT="$4"

$BINARY $PRM >${OUTPUT} 2>&1
ret=$?
if [[ "$EXPECT" == "$ret" ]]
then
  exit 0
fi

exit 1
