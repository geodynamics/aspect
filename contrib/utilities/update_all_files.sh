#!/bin/bash

# This script executes 'update_source_files.sh' and 'update_prm_files.sh' for
# all applicable files in the ASPECT repository. It should usually not have any
# effect, because all files should be on the current version, but if you change
# the code, or add new files it can be handy to execute this script.

UTIL_DIR=`dirname $0`
BASE_DIR=`cd ${UTIL_DIR}/../..;pwd`
echo "Scanning ${BASE_DIR} for changes..."

# Update source files
SOURCE_FILES=`find $BASE_DIR -type f \( -name *.cc -or -name *.h \) -and -not -name *.bak | grep -v doc`
bash ${UTIL_DIR}/update_source_files.sh $SOURCE_FILES

# Update prm files
PRM_FILES=`find $BASE_DIR -type f -name *.prm* -not -name *update_script* -not -name *.bak`
bash ${UTIL_DIR}/update_prm_files.sh $PRM_FILES

# To remove the backup files that are created you will likely want to use the
# following command. It is commented out by default, because you should think
# carefully before automatically removing files.
# find ${BASE_DIR} -name *.bak -delete
