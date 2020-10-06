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


# Grab ASPECT_GENERATE_REFERENCE_OUTPUT from the environment. If set to
# something (not ""), do not run tests normally but generate reference output
# instead. Also see the generate_reference_output make target and the file
# ./cmake/generate_reference_output.sh
GENERATE_REFERENCE_OUTPUT=0
if [ -n "$ASPECT_GENERATE_REFERENCE_OUTPUT" ]; then
  GENERATE_REFERENCE_OUTPUT=1
  echo "generating reference output" 1>&2
fi

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
	${DIFF_EXE} -V -a 1e-10 -r 1e-8 -s ' \t\r\n:<>=,;' \
	    ${REF_FILE} ${GEN_FILE} > ${DIFF_OUTPUT}.tmp
	;;
    *)
	"${DIFF_EXE}" \
	    ${REF_FILE} ${GEN_FILE} > ${DIFF_OUTPUT}.tmp
esac

if [ $? -ne 0 ]; then

  if [ "$GENERATE_REFERENCE_OUTPUT" -eq 1 ]; then
    echo "modifying $ORIGINAL_REF_FULL_PATH"
    cp ${GEN_FILE} $ORIGINAL_REF_FULL_PATH
    exit 0
  fi

  mv ${DIFF_OUTPUT}.tmp ${DIFF_OUTPUT}.failed
  echo "******* Error during diffing output results for ${PRETTY_TEST_AND_FILENAME}"
  echo "******* Results are stored in ${DIFF_OUTPUT}.failed"
  echo "******* Check ${ORIGINAL_GEN_FULL_PATH} ${ORIGINAL_REF_FULL_PATH}"
  nlines="`cat ${DIFF_OUTPUT}.failed | wc -l`"
  if test "$nlines" -ge 50 ; then
    echo "******* First 50 of $nlines lines of diffs are:"
  else
    echo "******* All of $nlines lines of diffs are:"
  fi
  head -n 50 ${DIFF_OUTPUT}.failed

  exit 1

else
  mv ${DIFF_OUTPUT}.tmp ${DIFF_OUTPUT}
fi

exit 0
