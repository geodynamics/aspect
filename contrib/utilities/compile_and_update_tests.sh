#!/bin/bash

# This script can modify the test results in the tests/
# directory to the official tester results, if docker is installed
# on the system. This script has to be started from the main ASPECT
# folder (the one containing include/, cmake/, source/, and tests/).
# Navigate to that folder and run the following command in a terminal:
#
# docker run -v $PWD:/home/dealii/aspect --name=aspect-tester --rm -it geodynamics/aspect-tester:focal-dealii-9.4-v3 bash /home/dealii/aspect/contrib/utilities/compile_and_update_tests.sh

# This command executes the following shell script *inside* the docker container
# that contains the official ASPECT test system. Note that by mounting your
# ASPECT folder into the container you are actually changing reference test
# results on the host system (i.e. your computer) instead of inside the
# container, and you also keep the build folder outside of the container in
# the new directory 'tester-build'.

echo "Starting official tester and compiling source ..."

SRC_PATH=`dirname $0`
SRC_PATH=`cd $SRC_PATH/../..;pwd`

mkdir -p ${SRC_PATH}/tester-build
cd ${SRC_PATH}/tester-build
cmake -G "Ninja" -D ASPECT_TEST_GENERATOR=Ninja -D ASPECT_RUN_ALL_TESTS=ON -D ASPECT_PRECOMPILE_HEADERS=ON $SRC_PATH
ninja
ASPECT_TESTS_VERBOSE=1 bash ${SRC_PATH}/cmake/generate_reference_output.sh
