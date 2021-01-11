#!/bin/bash

# This script is used by the test framework of ASPECT to provide the GNU bash
# "timeout" routine to abort hanging tests. Sadly, it is not available on OSX
# by default. For OSX, we ignore the timeout and run the command directly.  On
# Linux machines where "timeout" exists, this script behaves identically.

if command -v timeout &> /dev/null
then
    # timeout exists, so just pass along all arguments
    timeout "$@"
else
    # Work-around: remove the specified timeout time and run
    shift
    "$@"
fi
