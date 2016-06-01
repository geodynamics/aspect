#!/bin/bash

# This script is executed to run ASPECT with a given .prm test and process the
# return value of ASPECT correctly.

# usage: EXPECT_FAIL ./aspect test.prm output-name
#
# ASPECT is expected to return a non-zero return code if EXPECT_FAIL=1 and
# 0 (=success) if EXPECT_FAIL=0

EXPECT_FAIL="$1"
BINARY="$2"
PRM="$3"
OUTPUT="$4"

$BINARY $PRM >${OUTPUT} 2>&1
ret=$?
if [[ ( "$EXPECT_FAIL" == "0" && "$ret" == "0" ) 
	    ||
      ( "$EXPECT_FAIL" == "1" && "$ret" != "0" ) ]]
then
  exit 0
fi

exit 1
