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
  # success: ASPECT succeeded (or failed if we expect it to fail)
  exit 0
fi

# ASPECT failed (or succeeded even though it should not have), so report
# failure. Before we do so, dump the screen-output to the screen:
cat ${OUTPUT}
exit 1
