#!/bin/bash

# Script for diffing a certain output file of a test and printing nice output
#
# While we could put this code into CMakeLists.txt, cmake prints the whole
# command that failed, which will be really long and ugly. This way, only the
# line executing this script will show up.

# Call this script with 4 parameters:

# 1. name/path of diff executable
DIFF_EXE=$1

# 2. short name of the reference file: "testname/filename"
PRETTY_TEST_AND_FILENAME=$2

# 3. absolute cmake source directory
CMAKE_CURRENT_SOURCE_DIR=$3

# 4. absolute cmake binary directory
CMAKE_CURRENT_BINARY_DIR=$4


#
# now generate filenames with full paths
#

# the reference file to compare with: "filename.cmp.notime"
REF_FILE=${CMAKE_CURRENT_BINARY_DIR}/output-${PRETTY_TEST_AND_FILENAME}.cmp.notime

# filename to compare: "filename.notime"
GEN_FILE=${CMAKE_CURRENT_BINARY_DIR}/output-${PRETTY_TEST_AND_FILENAME}.notime

# filename to write the diff output into: "filename.diff"
DIFF_OUTPUT=${CMAKE_CURRENT_BINARY_DIR}/output-${PRETTY_TEST_AND_FILENAME}.diff

# full path of the original (unprocessed reference file)
ORIGINAL_REF_FULL_PATH=${CMAKE_CURRENT_SOURCE_DIR}/${PRETTY_TEST_AND_FILENAME}

# full path of the generated filename (without .notime)
ORIGINAL_GEN_FULL_PATH=${CMAKE_CURRENT_BINARY_DIR}/output-${PRETTY_TEST_AND_FILENAME}


#
# Finally do the work.
#

rm -f ${DIFF_OUTPUT}.failed ${DIFF_OUTPUT}

case ${DIFF_EXE} in
    *numdiff)
	${DIFF_EXE} -a 1e-6 -r 1e-8 -s ' \t\n:<>=,;' \
	    ${REF_FILE} ${GEN_FILE} > ${DIFF_OUTPUT}.tmp
	;;
    *)
	"${DIFF_EXE}" \
	    ${REF_FILE} ${GEN_FILE} > ${DIFF_OUTPUT}.tmp
esac

if [ $? -ne 0 ]; then
  mv ${DIFF_OUTPUT}.tmp ${DIFF_OUTPUT}.failed
  echo "******* Error during diffing output results for ${PRETTY_TEST_AND_FILENAME}"
  echo "******* Results are stored in ${DIFF_OUTPUT}.failed"
  echo "******* Check ${ORIGINAL_GEN_FULL_PATH} ${ORIGINAL_REF_FULL_PATH}"
  echo "******* 50 lines of diffs are:"
  head -n 50 ${DIFF_OUTPUT}.failed
  exit 1
else
  mv ${DIFF_OUTPUT}.tmp ${DIFF_OUTPUT}
fi

exit 0
