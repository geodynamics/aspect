#!/bin/bash

# A script to update .prm files to the current ASPECT naming scheme. This
# script correctly updated all files within the official development version of
# ASPECT, but it is not guaranteed to work for all possible names in user
# files. Consequently a backup of files is created, and the changes by this
# script should be investigated to ensure a correct renaming.
#
# Usage for a parameter file named FILENAME
# (possibly containing wildcards such as '*.prm'):
#
# bash update_prm_files.sh FILENAME

FOLDER=`dirname $0`
SCRIPT_FOLDER=${FOLDER}/update_scripts/prm_files

# Create the backup
for parameter_file in "$@"; do
  cp $parameter_file ${parameter_file}.bak
done

# Run all update scripts in ${SCRIPT_FOLDER} on all files.
# In order to make the -i command portable between Linux and MacOS,
# we need to provide a backup file-ending (.tmp) and then later
# remove the backup file. We can not use this file
# instead of the .bak file, because it is overwritten by every script,
# and so it is only a backup of the last execution.
for script in `find ${SCRIPT_FOLDER} -maxdepth 1 -name \*.sed`; do
  sed -i.tmp -f $script "$@"
done

for script in `find ${SCRIPT_FOLDER} -maxdepth 1 -name \*.pl`; do
  for file in "$@" ; do
    cat "$file" | perl $script > "$file.tmp"
    mv "$file.tmp" "$file"
  done
done

for script in `find ${SCRIPT_FOLDER} -maxdepth 1 -name \*.py`; do
  for file in "$@" ; do
    python3 $script $file "$file.tmp"
    mv "$file.tmp" "$file"
  done
done
