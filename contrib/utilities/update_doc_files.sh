#!/bin/bash

# This script fixes some artifacts in the automatically converted documentation

UTIL_DIR=`dirname $0`
BASE_DIR=`cd ${UTIL_DIR}/../..;pwd`
echo "Replacing files in ${BASE_DIR} ..."

cd $BASE_DIR

PROCESS_DIR="benchmarks cookbooks"
SOURCE_FILES=`find cookbooks benchmarks -name '*.md' $(printf "! -name %s " $(cat contrib/utilities/excluded_doc_files))`
#echo $SOURCE_FILES

for file in $SOURCE_FILES; do
  echo $file
  perl -0777 $BASE_DIR/contrib/utilities/update_doc_files.pl $file > "$file.tmp"
  mv "$file.tmp" "$file"
done