#!/bin/bash

# This script is executed by ctest when ASPECT_RUN_ALL_TESTS is disabled. Due
# to a bug in cmake we can not print this info directly inside
# CTestCustom.ctest if you are running cmake 2.8.11.x so that is why we do it
# in a separate shell script.

echo ""
echo ""
echo "---------------------------------------------------------------"
echo "Completed running the minimal testsuite. Please check the test"
echo "results above. If any test failed, run"
echo "    ctest -V"
echo "to get detailed error messages. Contact the ASPECT mailing list"
echo "with that information if you need help."
echo ""
echo "To run the complete testsuite, reconfigure your project with"
echo "the option ASPECT_RUN_ALL_TESTS=ON or by running"
echo "    make setup_tests"
echo "Be aware that installing the 'numdiff' tools is strongly"
echo "recommended before running the full testsuite. Even then, the"
echo "output is likely different on your machine due to configuration"
echo "differences, so do not expect all tests to succeed on your"
echo "machine."
echo "---------------------------------------------------------------"
