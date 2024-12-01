#!/bin/bash

# A script to update .cc and .h files to the current ASPECT naming scheme. This
# script correctly updated all files within the official development version of
# ASPECT, but it is not guaranteed to work for all possible names in user
# files. Consequently a backup of files is created, and the changes by this
# script should be investigated to ensure a correct renaming.
#
# Usage for a source file named FILENAME
# (possibly containing wildcards such as '*.cc'):
#
# bash update_source_files.sh FILENAME

FOLDER=`dirname $0`
SCRIPT_FOLDER=${FOLDER}/update_scripts/source_files

# Create the backup
for source_file in "$@"; do
  cp $source_file ${source_file}.bak
done

# Run all update scripts in ${SCRIPT_FOLDER} on all files.
# In order to make the -i command portable between Linux and MacOS,
# we need to provide a backup file-ending (.tmp) and then later
# remove the backup file. We can not use this file
# instead of the .bak file, because it is overwritten by every script,
# and so it is only a backup of the last execution.
for script in `find ${SCRIPT_FOLDER} -maxdepth 1 -name "*.sed"`; do
  sed -i.tmp -f $script "$@"
done

for script in `find ${SCRIPT_FOLDER} -maxdepth 1 -name "*.pl"`; do
  for file in "$@" ; do
    cat "$file" | perl $script > "$file.tmp"
    mv "$file.tmp" "$file"
  done
done
